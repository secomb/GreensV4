/*****************************************************
eval - Evaluate solute field from source strengths.  TWS November 07.
Modified to include non-diffusible solutes.  May 2010.
Version 3.0, May 17, 2011.
Version 4.0, March 1, 2018.
******************************************************/
#include <math.h>
#include "nrutil.h"
#include <stdio.h>

float *eval(int slsegdiv, float req, float *x)
{
	extern int mxx, myy, mzz, nnt, nnv, nseg, nnod, nsp, is2d;
	extern int *mainseg, **tisspoints, *permsolute, *diffsolute, ***nbou;
	extern float fac, w2d, r2d;
	extern float *axt, *ayt, *azt, *ds, *diff, *g0, *y, *p;
	extern float **qt, **start, **scos, **qv, **ax, **pt;

	float dist2, gtt, gtv, lamx, lamy, lamz, r2d2 = SQR(r2d), req2 = SQR(req);
	int i, j, k, ii, jj, kk, iseg, itp, isp;

	for (isp = 1; isp <= nsp; isp++) p[isp] = g0[isp];	//initialize to g0
	for (itp = 1; itp <= nnt; itp++) {	//add contributions from tissue sources
		dist2 = SQR(x[1] - axt[tisspoints[1][itp]])
			+ SQR(x[2] - ayt[tisspoints[2][itp]])
			+ SQR(x[3] - azt[tisspoints[3][itp]]);
		if (dist2 <= req2) {
			if (is2d) gtt = fac / w2d * (log(r2d2 / req2) + 1. - dist2 / req2);
			else gtt = fac * (1.5 - 0.5*dist2 / req2) / req;
		}
		else {
			if (is2d) gtt = fac / w2d * log(r2d2 / dist2);
			else gtt = fac / sqrt(dist2);
		}
		for (isp = 1; isp <= nsp; isp++)	if (diffsolute[isp]) p[isp] += gtt / diff[isp] * qt[itp][isp];
	}
	
	for (i = 1; i <= nnv; i++) {	//add contributions from vessel sources.  Subdivide subsegments.
		iseg = mainseg[i];			//Note that vessel point is at midpoint of subsegment
		for (k = 1; k <= slsegdiv; k++) {
			for (j = 1; j <= 3; j++)	y[j] = ax[j][i] + scos[j][iseg] * ds[iseg] * (-0.5 + (k - 0.5) / slsegdiv);
			dist2 = SQR(x[1] - y[1]) + SQR(x[2] - y[2]) + SQR(x[3] - y[3]);
			if (dist2 <= req2) {
				if (is2d) gtv = fac / w2d * (log(r2d2 / req2) + 1. - dist2 / req2);
				else gtv = fac * (1.5 - 0.5*dist2 / req2) / req;
			}
			else {
				if (is2d) gtv = fac / w2d * log(r2d2 / dist2);
				else gtv = fac / sqrt(dist2);
			}
			for (isp = 1; isp <= nsp; isp++)	if (permsolute[isp]) p[isp] += gtv / diff[isp] * qv[i][isp] / slsegdiv;
		}
	}
	//for non-diffusible solute, calculate by interpolating values at tissue points.  May 2010.
	for (isp = 1; isp <= nsp; isp++)	if (diffsolute[isp] == 0) {
		i = 0;
		while (x[1] > axt[i + 1] && i < mxx) i++;
		if (i == 0) lamx = 1.;
		else if (i == mxx) lamx = 0.;
		else lamx = (x[1] - axt[i]) / (axt[i + 1] - axt[i]);
		j = 0;
		while (x[2] > ayt[j + 1] && j < myy) j++;
		if (j == 0) lamy = 1.;
		else if (j == myy) lamy = 0.;
		else lamy = (x[2] - ayt[j]) / (ayt[j + 1] - ayt[j]);
		k = 0;
		while (x[3] > azt[k + 1] && k < mzz) k++;
		if (k == 0) lamz = 1.;
		else if (k == mzz) lamz = 0.;
		else lamz = (x[3] - azt[k]) / (azt[k + 1] - azt[k]);
		p[isp] = 0.;
		for (ii = 0; ii <= 1; ii++) for (jj = 0; jj <= 1; jj++) for (kk = 0; kk <= 1; kk++)
			if (i + ii >= 1 && i + ii <= mxx && j + jj >= 1 && j + jj <= myy && k + kk >= 1 && k + kk <= mzz) {
				itp = nbou[i + ii][j + jj][k + kk];
				if (itp != 0) p[isp] += (1. - lamx + ii * (2.*lamx - 1))*(1. - lamy + jj * (2.*lamy - 1))
					*(1. - lamz + kk * (2.*lamz - 1))*pt[itp][isp];
			}
	}
	return p;
}