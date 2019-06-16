/************************************************************************
convect for greens - depends on solute number isp
TWS August 2010
Set up convective fluxes and alpha matrix
Version 3.0 May 17, 2011.
Version 4.0, March 1, 2018.
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void blood(float c, float hem, float *p, float *pp);
float bloodconc(float p, float h);
float bloodconcp(float p, float h);

void convect(int isp)
{
	extern int nnodbc, nseg, nnodfl, nodsegm, nsp, nnv, nsegfl;
	extern int *bcnod, *nodtyp, *nodrank, *nodout, *segtyp, *permsolute, *oxygen, *nspoint, *istart;
	extern int *segname, *nodname, **nodseg, *mainseg;
	extern float *bifpar, *hd, *qq, *q, *bchd, *diam, **pv, **bcp, *ds, *segc, **cv, **qv, flowfac, *solutefac;
	extern float **al;

	int i, j, k, ii, jj, inod, iseg, jseg, in, isegk, nodt, nin, nout, ineg, ihigh;
	float fluxsumin, pb, pp;
	float sumin, sumout, hdsumin, hdsumout;

	int *isegkk; //added June 2013 to check for errors in segment sequence
	isegkk = ivector(1, nseg);
	for (iseg = 1; iseg <= nseg; iseg++) isegkk[iseg] = 0;

	isegk = 0;	//number of segments processed
	for (i = 1; i <= nnv; i++) for (j = 1; j <= nnv; j++) al[i][j] = 0.;
	for (iseg = 1; iseg <= nseg; iseg++) segc[iseg] = 0.;

	//set convective fluxes in segments connected to inflow boundary nodes, segc is the convective flux
	for (j = 1; j <= nnodbc; j++) {
		inod = bcnod[j];
		if (nodout[inod] == 1) {
			iseg = nodseg[1][inod];
			//values modified according to VaryParams.dat
			if (oxygen[isp] == 1) segc[iseg] = bloodconc(bcp[j][isp] * solutefac[isp], hd[iseg])*qq[iseg] * flowfac;
			else segc[iseg] = bcp[j][isp] * solutefac[isp] * qq[iseg] * flowfac;
			isegkk[iseg] = 1;
		}
	}
	ineg = 0;
	ihigh = 0;

	for (in = 1; in <= nnodfl; in++) {	//scan all nodes in downstream order
		inod = nodrank[in];
		nodt = nodtyp[inod];
		nout = nodout[inod];
		nin = nodt - nout;
		if (nodt > 1) {	//don't do this part for network boundary nodes
			sumin = 0.;
			hdsumin = 0.;
			fluxsumin = 0.;
			for (ii = nout + 1; ii <= nodt; ii++) { //inflows
				iseg = nodseg[ii][inod];
				if (isegkk[iseg] == 0) printf("*** Error: wrong segment sequence in convect, segment %i ***\n", iseg);
				sumin += qq[iseg] * flowfac;
				hdsumin += qq[iseg] * flowfac*hd[iseg];
				fluxsumin += segc[iseg];
			}
			//calculate solute level going into node
			if (oxygen[isp] == 1) blood(fluxsumin / sumin, hdsumin / sumin, &pb, &pp);
			else pb = fluxsumin / sumin;
			//assign solute levels going out of node
			sumout = 0.;
			hdsumout = 0.;
			for (ii = 1; ii <= nout; ii++) {		//outflows	
				iseg = nodseg[ii][inod];
				isegkk[iseg] = 1;
				sumout += qq[iseg] * flowfac;	//check conservation of flow and hematocrit
				hdsumout += qq[iseg] * flowfac*hd[iseg];
				if (oxygen[isp] == 1) segc[iseg] = bloodconc(pb, hd[iseg])*qq[iseg] * flowfac;
				else segc[iseg] = pb * qq[iseg] * flowfac;
				if (q[iseg] >= 0.) i = istart[iseg];
				else i = istart[iseg] + nspoint[iseg] - 1;
				for (jj = nout + 1; jj <= nodt; jj++) {	//inflows
					jseg = nodseg[jj][inod];
					if (q[jseg] >= 0.) j = istart[jseg] + nspoint[jseg] - 1;
					else j = istart[jseg];
					al[i][j] = qq[iseg] * flowfac / sumin;	//calculate alpha values across node
					if (oxygen[isp] == 1 && nout > 1)	//if nout=1, hd[iseg] = hdsumin/sumin
						al[i][j] *= bloodconcp(pb, hd[iseg]) / bloodconcp(pb, hdsumin / sumin);
					for (k = 1; k <= nnv; k++) al[i][k] += al[i][j] * al[j][k];//calculate other alpha values
				}
			}
			if (sumin + sumout != 0.) if (fabs(sumin - sumout) / (sumin + sumout) > 0.01)
				printf("*** Error: Flow conservation violation at node %i\n", nodname[inod]);
			if (hdsumin + hdsumout != 0.) if (fabs(hdsumin - hdsumout) / (hdsumin + hdsumout) > 0.01)
				printf("*** Error: Hematocrit conservation violation at node %i\n", nodname[inod]);
		}
		//subsegments of outflow segments - convective fluxes and alpha matrix, including network boundary nodes
		for (ii = 1; ii <= nout; ii++) {				//outflows
			iseg = nodseg[ii][inod];
			for (jj = 0; jj < nspoint[iseg]; jj++) {	//convective fluxes
				if (q[iseg] >= 0.) i = istart[iseg] + jj;
				else i = istart[iseg] + nspoint[iseg] - jj - 1;
				segc[iseg] -= qv[i][isp] / 2.;
				cv[i][isp] = segc[iseg] / qq[iseg] / flowfac;
				segc[iseg] -= qv[i][isp] / 2.;
			}
			for (jj = 1; jj < nspoint[iseg]; jj++) {	//alpha matrix
				if (q[iseg] >= 0.) {
					j = istart[iseg] + jj - 1;
					i = istart[iseg] + jj;
				}
				else {
					j = istart[iseg] + nspoint[iseg] - jj;
					i = istart[iseg] + nspoint[iseg] - jj - 1;
				}
				al[i][j] = 1.;
				for (k = 1; k <= nnv; k++) al[i][k] += al[i][j] * al[j][k];
			}
		}
		isegk += nout;
	}
	if (isegk != nsegfl) printf("*** Error in convect, %i of %i segments processed\n", isegk, nseg);
	for (i = 1; i <= nnv; i++) al[i][i] = 0.5;
	free_ivector(isegkk, 1, nseg);
}