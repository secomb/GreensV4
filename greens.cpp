/************************************************************************
Green's function approach for multiple reacting species.
T.W. Secomb July 2007 - based on Greens.f by R. Hsu.
See http://www.physiology.arizona.edu/people/secomb/greens.html

Variable array sizes and bounds, using Numerical Recipes utilities.
Tissue-vessel and tissue-tissue matrices computed on the fly to save memory.
No nondimensionalization.  Lengths, diameters in microns, times in s.
Flows in nanoliters/min
Oxygen concentrations in cm^3/cm^3
Consumption rates in cm^3/cm^3/s - changed in V2
Partial pressures in mmHg
Diffusivities in cm2/s, converted to micron^2/s

Special parameters for oxygen
p50, fn --- parameters in the Hill equation
cs --- red blood cell oxygen binding capacity in cm^3 O2/cm^3
alphab --- average solubility in blood in cm^3 O2/cm^3/mmHg
gamma1 --- intravascular resistance, varies with vessel diameter, in mmHg.cm.s/cm^3 O2

Main variables:
  gvv --- Green's function matrix for vessels
  mat --- matrix for vessel strengths
  al --- matrix giving dependence of vessel convective fluxes on source strengths
  lseg --- segment length
  ds --- subsegment length
  qv --- oxygen efflux from subsegment
  pv --- PO2 in the subsegment
  cv --- oxygen concentration in the subsegment
  qt --- tissue source strength
  pvt --- PO2 on vessels due to source terms in tissue
  ptv --- PO2 in tissue due to source terms on vessels
  q --- flow rate, qq = abs(q)

Version 2.0, May 2010.
With 9/08 updates.  New vessel-vesel interaction coefficient. January 2009
With alternate terms for 2D version.  May 2009
With g0 computed as part of linear system, for permeable solutes.  March 2010
  g0method = 1:  include g0 in linear system to be solved - fastest method *****
  g0method = 2:  theoretical estimates of dqsum/dg0 - was used in Version 1
For impermeable solutes, method 2 is always used
With choice of Gauss-Jordan, LU or biconjugate gradient linear solvers.  March 2010
  linmethod = 1:  Gaussian elimination - was used in Version 1
  linmethod = 2:  LU decomposition
  linmethod = 3:  biconjugate gradient (iterative) - fastest method *****
Does not require that species 1 is oxygen. Allows for non-diffusing solutes.  April 2010
Creates log file. April 2010.
During tissue loop, scales vessel sources so that qvsum = qtsum, for faster convergence.  April 2010.
Modified for compatibility with Mac XCode compiler.  April 2010.
Includes intravascular resistance for all solutes.  May 2010.
Includes non-diffusible solutes.  May 2010.

Version 3.0, May 17, 2011.
Uses convect instead of genalpha and genalphahd.
This gives improved convergence if hematocrit is non-uniform in network

Version 4.0, March 1, 2018.
Includes cmgui to generate files needed for visualization using CMGUI
In soluteparams.dat, the parameter for total inflow is replaced by a flow factor
Includes capability to run multiple cases, controlled by VaryParams.dat
Includes capability to compute dependent variables (e.g. cell survival) after
  greens has run, as specified by PostGreensParams.dat
Problem-specific c code (used as #include files) is needed as follows:
  tissrate.cpp.dat specifies dependence of reaction rates on concentrations
  postgreens.cpp.dat specifies computation of dependent variables
Produces a summary output file from multiple runs (summary.out)
Output files (TissueLevels.dat, TissueSources.dat, VesselLevels.dat, VesselSources.dat)
  are formatted differently, with one line per tissue point or vessel point
Previous combined output file (GreensRes.out) is no longer generated
Capability to restart a run has been removed (was not functional)
All output files and a copy of input files are placed in a folder called "Current"
Algorithm for non-diffusible solutes is improved by using reaction rates based on
  updated levels of other solutes at each iteration of Newton method

**************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

void putrank(void);
void initgreens();
void blood(float c, float hem, float *p, float *pp);
float bloodconc(float p, float h);
void tissrate(int nsp, float *p, float *mtiss, float *mptiss);
void convect(int isp);
void testconvect(int isp);
float *eval(int slsegdiv, float req, float *x);

void gaussj(double **a, int n, double **b, int m);
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double *b);
double bicgstab(double **a, double *b, double *x, int n, double eps, int itmax);

void greens(void)
{
	extern int nmaxvessel, nmaxtissue, nmax, g0method, linmethod, kmain, imain;
	extern int mxx, myy, mzz, nnt, nnv, nseg, nsp, nnodbc;
	extern int is2d; //needed for 2d version
	extern int *mainseg, **tisspoints, *permsolute, *segtyp;
	extern int *segname, *nspoint, *istart, *bcnod, *lowflow;
	extern int *errvesselcount, *errtissuecount;
	extern int *imaxerrvessel, *imaxerrtissue, *nresis, *indx;
	extern int *oxygen, *diffsolute, **nodseg, ***nbou;

	extern float p50, cs, req, totalq, q0fac, fac, flowfac, lowflowcrit, errfac, v, vol, vdom, tlength, pi1;
	extern float tlengthq, tlengthqhd;
	extern float alx, aly, alz, w2d, r2d; //needed for 2d version
	extern float *axt, *ayt, *azt, *ds, *diff, *pmin, *pmax, *pmeant, *pmeanv, *psdt, *psdv, *g0, *g0fac, *g0facnew, *pref, *dtmin;
	extern float *diam,*rseg,*qdata,*q,*qq,*hd,*bchd,*qvtemp,*qvfac;
	extern float *x, *lseg, *mtiss, *mptiss, *dqvsumdg0, *dqtsumdg0;
	extern float *epsvessel, *epstissue, *eps, *errvessel, *errtissue, *p;
	extern float *rhs, *rhstest, *g0old, *ptt, *ptpt, *qtsum, *qvsum, *intravascfac;
	extern float **tissparam, **start,**scos,**ax,**bcp;
	extern float **qv, **qt, **pv, **pev, **pt, **resisdiam, **resis;
	extern float **qvseg,**pvseg,**pevseg, **pvt,**pvprev,**qvprev,**cv,**dcdp;
	extern float **ptprev,**ptv,**gamma1,**cv0,**conv0, **gvv,**end,**al;
	extern double **mat, **rhsg, *rhsl, *matx;
	extern float ***dtt;


	int i, j, k, ix, iy, iz, jx, jy, jz, iseg, nt, ineg, ihigh, isp, imaxerr;
	int ixdiff, iydiff, izdiff, isn, jseg, ktissue, kvessel, itp, jtp, convflag, convflagt, convflagv;
	int greensverbose = 1;
	int bicgstabit = 5000; //parameter for biconjugate gradient method
	double dd,bicgstaberr = 0.0001; //parameter for biconjugate gradient method

	float x11, x22, x33, duration, rsegmax, dsmax, gvarfac;
	float gtt, gtv, gvt, disp2, ds2, dist, d, de2, dtave, den, dqsumdg0;
	float dif, err, qhdcrit, pathlength, tlengthdiam;
	float lam, lam3d, lam2d;
	char fname[80];
	FILE *ofp,*ofp1;
	clock_t tstart, tfinish, tstart1, tfinish1;

	w2d = alz;  //needed for 2d
	r2d = sqrt(SQR(alx) + SQR(aly) + SQR(alz));
	lam3d = 0.2;//	lam3d = 0.2;	//Underrelax iteration of tissue levels
	lam2d = 0.1;	//Use this one for 2D permeable solutes only
	if (is2d) lam = lam2d;
	else lam = lam3d;

	//setup mainseg (must be done after setuparrays2)
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5)
		for (i = 0; i < nspoint[iseg]; i++) mainseg[istart[iseg] + i] = iseg;
	//identify vessel points
	for (i = 1; i <= nnv; i++) {
		iseg = mainseg[i];
		isn = i - istart[iseg];
		for (j = 1; j <= 3; j++) ax[j][i] = start[j][iseg] + scos[j][iseg] * ds[iseg] * (isn + 0.5);
	}
	//index tissue points
	for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) for (k = 1; k <= mzz; k++) {
		nt = nbou[i][j][k];
		if (nt > 0) {
			tisspoints[1][nt] = i;
			tisspoints[2][nt] = j;
			tisspoints[3][nt] = k;
		}
	}
	//calculate the distance of tissue points to the nearest vessel
	dtave = 0.;
	for (itp = 1; itp <= nnt; itp++) {
		i = tisspoints[1][itp];
		j = tisspoints[2][itp];
		k = tisspoints[3][itp];
		dtmin[itp] = 1.e6;
		for (jseg = 1; jseg <= nseg; jseg++) if (segtyp[jseg] == 4 || segtyp[jseg] == 5) {
			x11 = (axt[i] - start[1][jseg])*scos[2][jseg] - (ayt[j] - start[2][jseg])*scos[1][jseg];
			x22 = (ayt[j] - start[2][jseg])*scos[3][jseg] - (azt[k] - start[3][jseg])*scos[2][jseg];
			x33 = (azt[k] - start[3][jseg])*scos[1][jseg] - (axt[i] - start[1][jseg])*scos[3][jseg];
			disp2 = SQR(x11) + SQR(x22) + SQR(x33);
			ds2 = SQR(axt[i] - start[1][jseg]) + SQR(ayt[j] - start[2][jseg]) + SQR(azt[k] - start[3][jseg]);
			de2 = SQR(axt[i] - end[1][jseg]) + SQR(ayt[j] - end[2][jseg]) + SQR(azt[k] - end[3][jseg]);
			if (FMAX(ds2, de2) - disp2 > SQR(lseg[jseg])) d = FMAX(sqrt(FMIN(ds2, de2)) - rseg[jseg], 0.);
			else d = FMAX(sqrt(disp2) - rseg[jseg], 0.);
			if (d < dtmin[itp]) dtmin[itp] = d;
		}
		dtave += dtmin[itp];
	}
	dtave = dtave / nnt;
	vdom = nnt * vol;
	tlength = 0.;
	tlengthq = 0.;
	tlengthqhd = 0.;
	tlengthdiam = 0.;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
		q[iseg] = qdata[iseg] * q0fac;		//scaling by q0fac
		qq[iseg] = fabs(q[iseg]);
		tlength += lseg[iseg];
		tlengthq += lseg[iseg] * qq[iseg];
		tlengthqhd += lseg[iseg] * qq[iseg] * hd[iseg];
		tlengthdiam += lseg[iseg] * diam[iseg];
	}
	den = sqrt(vdom / tlength);
	if (greensverbose) {
		printf("Average distance from tissue node to the nearest vessel = %f\n", dtave);
		printf("Vessel length) = %f\n", tlength);
		printf("Sqrt(Tissue Volume/vessel length) = %f\n", den);
		printf("Capillary density = %8.1f /mm2\n", tlength / vdom * 1.e6);
		printf("Total inflow to network = %f nl/min (for flow values in network.dat)\n", totalq);
		printf("Perfusion = %f cm3/cm3/min (uncorrected for path length effect)\n", totalq*q0fac / vdom * 1.e6);
		pathlength = 0.;
		if (q0fac > 0.) for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5)
			pathlength += fabs(q[iseg])*lseg[iseg] / totalq / q0fac;
		printf("Flow-weighted path length = %f micron\n", pathlength);
		printf("Length-weighted mean diameter = %f micron\n", tlengthdiam/tlength);
	}
	//Calculate intravascular or wall transport resistance.  Zero unless specified in intravascfac.dat.
	//If not oxygen, assume value from data is 1/(wall permeability in um/s)
	for (isp = 1; isp <= nsp; isp++) {
		for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
			gamma1[iseg][isp] = 0.;
			if (nresis[isp] != 0) {
				gamma1[iseg][isp] = resis[1][isp];
				for (j = 2; j <= nresis[isp]; j++) if (diam[iseg] <= resisdiam[j][isp] && diam[iseg] > resisdiam[j - 1][isp])
					gamma1[iseg][isp] = resis[j - 1][isp] + (resis[j][isp] - resis[j - 1][isp])
					*(diam[iseg] - resisdiam[j - 1][isp]) / (resisdiam[j][isp] - resisdiam[j - 1][isp]);
				if (diam[iseg] > resisdiam[nresis[isp]][isp]) gamma1[iseg][isp] = resis[nresis[isp]][isp];
				if (oxygen[isp] != 1) gamma1[iseg][isp] = gamma1[iseg][isp] / pi1 / diam[iseg];	
				gamma1[iseg][isp] = gamma1[iseg][isp] * intravascfac[isp];		//scaling from VaryParams.dat
			}
		}
	}
	//vessel ~ vessel matrix elements gvv
	//Uses empirical fit to results from elliptical integral form for diagonal elements, updated 2009
	//if center of one segment lies within the other segment, calculate Gvv as for self-interaction term
	//based on larger radius and larger length
	for (i = 1; i <= nnv; i++) {
		iseg = mainseg[i];
		for (j = 1; j <= nnv; j++) {
			jseg = mainseg[j];
			//this section modified to give better behavior for 'hairpin' structures
			dist = sqrt(SQR(ax[1][j] - ax[1][i]) + SQR(ax[2][j] - ax[2][i]) + SQR(ax[3][j] - ax[3][i]));
			if (dist < FMAX(sqrt(ds[iseg] * rseg[iseg]), sqrt(ds[jseg] * rseg[jseg]))) {
				dsmax = FMAX(ds[iseg], ds[jseg]);
				rsegmax = FMAX(rseg[iseg], rseg[jseg]);
				//Version 2.0 of 3-D interaction coefficients for close or coincident segments.  See Sep. 2009 notes.
				gvarfac = 0.6*exp(-0.45*dsmax / rsegmax);
				//for distinct vessels close together, make distance rsegmax in following calculation, to improve convergence
				if (iseg != jseg) dist = rsegmax;
				gvv[i][j] = (1.298 / (1. + 0.297*powf(dsmax / rsegmax, 0.838)) - gvarfac * SQR(dist / rsegmax))*fac / rsegmax;
				//for 2D version, additional terms give effect of boundaries (reflection method)
				if (is2d) {
					gvv[i][j] -= fac * 2. / w2d * 0.926*SQR(1. - 1. / (1. + 0.36*dsmax / w2d))*powf(1. + dsmax / w2d, 0.27);
					gvv[i][j] += fac * 2. / w2d * (log(r2d / w2d + 0.5 + 0.27 / r2d * w2d) - 0.117);
				}
			}
			else {
				if (is2d) gvv[i][j] = fac * 2. / w2d * log(r2d / dist);
				else gvv[i][j] = fac / dist;
			}
		}
	}
	// tissue ~ vessel, vessel ~ tissue: compute matrix elements gvt, gtv on fly as needed
	// tissue ~ tissue: construct matrix of distances from a corner node
	//diagonal elements of tissue ~ tissue matrix gtt
	if (is2d) gtt = fac / w2d * (2.*log(r2d / req) + 0.5);
	else gtt = 1.2*fac / req;
	for (jx = 1; jx <= mxx; jx++) for (jy = 1; jy <= myy; jy++) for (jz = 1; jz <= mzz; jz++) {
		dist = sqrt(SQR(axt[1] - axt[jx]) + SQR(ayt[1] - ayt[jy]) + SQR(azt[1] - azt[jz]));
		if (jx*jy*jz != 1) {
			if (is2d) dtt[jx][jy][jz] = fac * 2. / w2d * log(r2d / dist);
			else dtt[jx][jy][jz] = fac / dist;
		}
		else dtt[jx][jy][jz] = gtt;
	}
	//detect and label vessels with very low q*hd - test if oxygen is one of the solutes
	qhdcrit = 0.;
	for (isp = 1; isp <= nsp; isp++) if (oxygen[isp] == 1) qhdcrit = lowflowcrit * tissparam[1][isp];
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
		lowflow[iseg] = 0;
		if (qq[iseg] * (hd[iseg] + 0.01) < qhdcrit) lowflow[iseg] = 1;//Added 0.01 to allow for high flow, zero hematocrit channels
	}
	initgreens();
	putrank();
	//	for(isp=1; isp<=nsp; isp++){//for testing purposes only
	//		convect(isp);
	//		testconvect(isp);	
	//	}

	//create log file
	ofp1 = fopen("Current/GreensLog.txt", "w");
	fprintf(ofp1, "GreensLog.txt\n");
	fclose(ofp1);
	tstart = clock();
	//********************** start of main loop *****************************
	for (kmain = 1; kmain <= nmax; kmain++) {
		tstart1 = clock();
		if (greensverbose) printf("\n----- kmain = %i -----\n", kmain);
		else printf(" %i", kmain);
		for (isp = 1; isp <= nsp; isp++) {
			for (itp = 1; itp <= nnt; itp++)	ptprev[itp][isp] = pt[itp][isp];
			if (permsolute[isp] == 1) for (i = 1; i <= nnv; i++) pvprev[i][isp] = pv[i][isp];
			g0old[isp] = g0[isp];
		}
		//********************** start of vessel loop *****************************
		//compute contribution pvt from tissue source strengths qt
		for (i = 1; i <= nnv; i++) {
			for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1) {
				if (g0method == 1 && permsolute[isp] == 1) pvt[i][isp] = 0.;
				else pvt[i][isp] = g0[isp];
			}
			for (itp = 1; itp <= nnt; itp++) {
				dist = sqrt(SQR(ax[1][i] - axt[tisspoints[1][itp]])
					+ SQR(ax[2][i] - ayt[tisspoints[2][itp]])
					+ SQR(ax[3][i] - azt[tisspoints[3][itp]]));
				if (dist <= req) {
					if (is2d) gvt = fac / w2d * (2.*log(r2d / req) + 1. - SQR(dist / req));
					else gvt = fac * (1.5 - 0.5*SQR(dist / req)) / req;
				}
				else {
					if (is2d) gvt = fac * 2. / w2d * log(r2d / dist);
					else gvt = fac / dist;
				}
				for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1) pvt[i][isp] += gvt / diff[isp] * qt[itp][isp];//permsolute
			}
		}
		//compute blood solute levels and PO2
		for (kvessel = 1; kvessel <= nmaxvessel; kvessel++) {
			convflagv = 1;
			for (isp = 1; isp <= nsp; isp++) {
				qvsum[isp] = 0.;
				dqvsumdg0[isp] = 0.;
				if (permsolute[isp] == 1) {
					ineg = 0;
					ihigh = 0;
					convect(isp);
					for (i = 1; i <= nnv; i++) {
						iseg = mainseg[i];
						qvprev[i][isp] = qv[i][isp];
						if (oxygen[isp] == 1) {
							if (lowflow[iseg] != 1) {//only do this if not a lowflow segment. 
								if (cv[i][isp] < 0.) {
									ineg++;
									if (ineg == 1 && greensverbose) printf("*** Warning: cblood is negative -%i", segname[iseg]);
									if (ineg > 1 && greensverbose) printf("-%i", segname[iseg]);
								}
								if (cv[i][isp] > bloodconc(150., hd[iseg])) {
									ihigh++;
									if (ihigh == 1 && greensverbose) printf("*** Warning: cblood is high +%i", segname[iseg]);
									if (ihigh > 1 && greensverbose) printf("+%i", segname[iseg]);
								}
								blood(cv[i][isp], hd[iseg], &pv[i][isp], &dcdp[i][isp]);
							}
						}
						else {
							pv[i][isp] = cv[i][isp];
							dcdp[i][isp] = 1.;
						}
					}
					if (ineg > 0 || ihigh > 0) if (greensverbose) printf("\n");
					for (i = 1; i <= nnv; i++) {	//generate linear system to be solved
						iseg = mainseg[i];
						rhs[i] = pv[i][isp] - pvt[i][isp];
						for (j = 1; j <= nnv; j++) {
							jseg = mainseg[j];
							mat[i][j] = gvv[i][j] / diff[isp] + al[i][j] / dcdp[i][isp] / qq[iseg] / flowfac;
							if (i == j) mat[i][j] += gamma1[iseg][isp] / ds[iseg];
							rhs[i] += al[i][j] * qv[j][isp] / dcdp[i][isp] / qq[iseg] / flowfac;
							if (oxygen[isp] == 1 && lowflow[mainseg[i]] == 1) {  //low q*hd
								if (i == j) mat[i][j] = 1.;
								else mat[i][j] = 0.;
							}
						}
						if (oxygen[isp] == 1 && lowflow[iseg] == 1) rhs[i] = qvprev[i][isp];  //low q*hd
					}
					//solve system of linear algebraic equations: Sum mat[i][j]*qv[j]= rhs[i]
					if (g0method == 1) {
						for (i = 1; i <= nnv; i++) {
							if (oxygen[isp] == 1 && lowflow[mainseg[i]] == 1) mat[i][nnv + 1] = 0.;  //low q*hd
							else mat[i][nnv + 1] = 1.;
							mat[nnv + 1][i] = 1.;
							mat[nnv + 1][nnv + 1] = 0.;
						}
						if (linmethod == 1) {
							for (i = 1; i <= nnv; i++) rhsg[i][1] = rhs[i];
							rhsg[nnv + 1][1] = -qtsum[isp];
							gaussj(mat, nnv + 1, rhsg, 1);
							for (i = 1; i <= nnv; i++) {
								qv[i][isp] = rhsg[i][1];
								qvsum[isp] += qv[i][isp];
							}
							g0[isp] = rhsg[nnv + 1][1];
						}
						if (linmethod == 2) {
							ludcmp(mat, nnv + 1, indx, &dd);
							for (i = 1; i <= nnv; i++) rhsl[i] = rhs[i];
							rhsl[nnv + 1] = -qtsum[isp];
							lubksb(mat, nnv + 1, indx, rhsl);
							for (i = 1; i <= nnv; i++) {
								qv[i][isp] = rhsl[i];
								qvsum[isp] += qv[i][isp];
							}
							g0[isp] = rhsl[nnv + 1];
						}
						if (linmethod == 3) {
							for (i = 1; i <= nnv; i++) rhsl[i] = rhs[i];
							rhsl[nnv + 1] = -qtsum[isp];
							for (i = 1; i <= nnv; i++) matx[i] = qv[i][isp];
							matx[nnv + 1] = g0[isp];
							bicgstab(mat, rhsl, matx, nnv + 1, bicgstaberr, bicgstabit);
							for (i = 1; i <= nnv; i++) {
								qv[i][isp] = matx[i];
								qvsum[isp] += qv[i][isp];
							}
							g0[isp] = matx[nnv + 1];
						}
					}
					if (g0method == 2) {
						if (linmethod == 1) {
							for (i = 1; i <= nnv; i++) {
								rhsg[i][1] = rhs[i];
								rhsg[i][2] = -1.;
							}
							gaussj(mat, nnv, rhsg, 2);
							for (i = 1; i <= nnv; i++) {
								qv[i][isp] = rhsg[i][1];
								qvsum[isp] += qv[i][isp];
								if (oxygen[isp] != 1 || lowflow[mainseg[i]] != 1) dqvsumdg0[isp] += rhsg[i][2];
							}
						}
						if (linmethod == 2) {
							ludcmp(mat, nnv, indx, &dd);
							for (i = 1; i <= nnv; i++) rhsl[i] = rhs[i];
							lubksb(mat, nnv, indx, rhsl);
							for (i = 1; i <= nnv; i++) {
								qv[i][isp] = rhsl[i];
								qvsum[isp] += qv[i][isp];
							}
							for (i = 1; i <= nnv; i++) rhsl[i] = -1.;
							lubksb(mat, nnv, indx, rhsl);
							for (i = 1; i <= nnv; i++)
								if (oxygen[isp] != 1 || lowflow[mainseg[i]] != 1) dqvsumdg0[isp] += rhsl[i];
						}
						if (linmethod == 3) {
							for (i = 1; i <= nnv; i++) {
								rhsl[i] = rhs[i];
								matx[i] = qv[i][isp];
							}
							bicgstab(mat, rhsl, matx, nnv, bicgstaberr, bicgstabit);
							for (i = 1; i <= nnv; i++) {
								qv[i][isp] = matx[i];
								qvsum[isp] += qv[i][isp];
							}
							for (i = 1; i <= nnv; i++) {
								rhsl[i] = -1.;
								matx[i] = 0.;
							}
							bicgstab(mat, rhsl, matx, nnv, bicgstaberr, bicgstabit);
							for (i = 1; i <= nnv; i++)
								if (oxygen[isp] != 1 || lowflow[mainseg[i]] != 1) dqvsumdg0[isp] += matx[i];
						}
					}
					//for low q*hd segments, calculate efflux based on change in extravascular oxygen level
					//save values in qvtemp to avoid influence on eval, update qv, underrelax
					for (i = 1; i <= nnv; i++) {
						iseg = mainseg[i];
						if (oxygen[isp] == 1 && lowflow[iseg] == 1) {
							for (j = 1; j <= 3; j++) x[j] = ax[j][i] - 0.5*scos[j][iseg] * ds[iseg];
							p = eval(1, req, x);
							p[isp] = FMAX(p[isp], 0.);
							pv[i][isp] = p[isp] / 2.;
							qvtemp[i] = q[iseg]*flowfac*bloodconc(p[isp],hd[iseg]);
							for (j = 1; j <= 3; j++) x[j] = ax[j][i] + 0.5*scos[j][iseg] * ds[iseg];
							p = eval(1, req, x);
							p[isp] = FMAX(p[isp], 0.);
							pv[i][isp] += p[isp] / 2.;
							qvtemp[i] -= q[iseg] * flowfac*bloodconc(p[isp], hd[iseg]);
							cv[i][isp] = bloodconc(pv[i][isp], hd[iseg]);
						}
					}
					for (i = 1; i <= nnv; i++) if (oxygen[isp] == 1 && lowflow[mainseg[i]] == 1)
						qv[i][isp] = 0.5*qvtemp[i] + 0.5*qvprev[i][isp];	//underrelax
					errvessel[isp] = 0.;
					imaxerr = 0;
					errvesselcount[isp] = 0;
					for (i = 1; i <= nnv; i++) {
						dif = qv[i][isp] - qvprev[i][isp];
						//If qv is large, use relative rather than absolute error 
						if (qv[i][isp] != 0.) dif = dif * FMIN(1., epsvessel[isp] / errfac / fabs(qv[i][isp]));
						if (fabs(dif) >= errvessel[isp]) {
							imaxerrvessel[isp] = mainseg[i];
							errvessel[isp] = fabs(dif);
						}
						if (fabs(dif) > epsvessel[isp]) errvesselcount[isp]++;
					}
					if (greensverbose) printf("Solute %i: qtsum = %f, qvsum = %f\n", isp, qtsum[isp], qvsum[isp]);
					if (greensverbose) printf("Solute %i: kvessel = %i, errvessel_q = %f, imaxerr = %i, g0 = %f\n",
						isp, kvessel, errvessel[isp], imaxerrvessel[isp], g0[isp]);
					if (errvesselcount[isp] > 0) convflagv = 0;
				}
			}
			if (convflagv) goto vesselconv;
		}
		for (isp = 1; isp <= nsp; isp++) if (errvesselcount[isp] > 0)
			if (greensverbose) printf("*** Warning: solute %i, %i vessel source strengths not converged\n",
				isp, errvesselcount[isp]);
	vesselconv:;
		//********************** end of vessel loop *****************************	
		//********************** start of tissue loop *****************************
		//Compute tissue source strengths iteratively by successive relaxation: updated qt values are immediately used.
		//Continually scales up qv values so that their sum equals updated sum of qt values.
		for (itp = 1; itp <= nnt; itp++) {
			for (isp = 1; isp <= nsp; isp++) ptv[itp][isp] = 0.;
			for (i = 1; i <= nnv; i++) {
				dist = sqrt(SQR(ax[1][i] - axt[tisspoints[1][itp]])
					+ SQR(ax[2][i] - ayt[tisspoints[2][itp]]) + SQR(ax[3][i] - azt[tisspoints[3][itp]]));
				if (dist <= req) {
					if (is2d) gtv = fac / w2d * (2.*log(r2d / req) + 1. - SQR(dist / req));
					else gtv = fac * (1.5 - 0.5*SQR(dist / req)) / req;
				}
				else {
					if (is2d) gtv = fac * 2. / w2d * log(r2d / dist);
					else gtv = fac / dist;
				}
				for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1) ptv[itp][isp] += gtv / diff[isp] * qv[i][isp];
			}
		}
		for (isp = 1; isp <= nsp; isp++) qvfac[isp] = 1.;
		for (ktissue = 1; ktissue <= nmaxtissue; ktissue++) {	//Scale all av, qvsum and ptv values so that qvsum = qtsum
			for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1 && g0method == 1) {
				qvfac[isp] = -qtsum[isp] / qvsum[isp];
				if (fabs(qvfac[isp]) > 2.) qvfac[isp] = 1.;  //avoid extreme values
				if (fabs(qvfac[isp]) < 0.5) qvfac[isp] = 1.;  //avoid extreme values
			}
			convflagt = 1;
			for (isp = 1; isp <= nsp; isp++) {
				qtsum[isp] = 0;
				errtissue[isp] = 0.;
				dqtsumdg0[isp] = 0.;
				errtissuecount[isp] = 0;
			}
			//----------------------------------------------------
			for (itp = 1; itp <= nnt; itp++) {	//contribution ptt from tissue source strengths qt
				ix = tisspoints[1][itp];
				iy = tisspoints[2][itp];
				iz = tisspoints[3][itp];
				for (isp = 1; isp <= nsp; isp++)	ptt[isp] = 0.;//all solutes
				for (jtp = 1; jtp <= nnt; jtp++) {
					jx = tisspoints[1][jtp];
					jy = tisspoints[2][jtp];
					jz = tisspoints[3][jtp];
					ixdiff = abs(ix - jx) + 1;
					iydiff = abs(iy - jy) + 1;
					izdiff = abs(iz - jz) + 1;
					for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1) ptt[isp] += dtt[ixdiff][iydiff][izdiff] * qt[jtp][isp];
				}
				for (isp = 1; isp <= nsp; isp++) ptpt[isp] = pt[itp][isp];
				for (isp = 1; isp <= nsp; isp++) {
						if(diffsolute[isp])		//diffusible solutes, underrelax 
							pt[itp][isp] = (1.-lam)*pt[itp][isp] + lam*(ptv[itp][isp]*qvfac[isp] + g0[isp] + ptt[isp]/diff[isp]);
					else {	//non-diffusible - use Newton method to solve for pt.
							tissrate(nsp,ptpt,mtiss,mptiss);	//update tissrates
						if (mptiss[isp] == 0.) printf("*** Error: mptiss[%i] = 0 at tissue point %i\n", isp, itp);
						else pt[itp][isp] -= mtiss[isp] / mptiss[isp];
					}
					ptpt[isp] = pt[itp][isp];
				}
				tissrate(nsp, ptpt, mtiss, mptiss);	//update tissrates
					for(isp=1; isp<=nsp; isp++){		//replace qt with value based on updated pt
					dif = mtiss[isp] * vol - qt[itp][isp];
					qt[itp][isp] = mtiss[isp] * vol;
					qtsum[isp] += qt[itp][isp];
					if (diffsolute[isp] == 1) dqtsumdg0[isp] += mptiss[isp] * vol;
					if (fabs(dif) > errtissue[isp]) {
						errtissue[isp] = fabs(dif);
						imaxerrtissue[isp] = itp;
					}
					if (fabs(dif) > epstissue[isp]) errtissuecount[isp]++;
				}
			}
			for (isp = 1; isp <= nsp; isp++) {
				if(greensverbose) printf("Solute %i: qtsum = %f, qvsum = %f\n",isp,qtsum[isp],qvsum[isp]*qvfac[isp]);
				if (greensverbose) printf("Solute %i: ktissue = %i, errtissue_q = %f, imaxerr = %i, g0 = %f\n",
					isp, ktissue, errtissue[isp], imaxerrtissue[isp], g0[isp]);
				if (errtissuecount[isp] > 0) convflagt = 0;
			}
			if (kmain > 1 && convflagt) goto tissueconv;  //force full number of iterations when kmain = 1.
		}
		for (isp = 1; isp <= nsp; isp++) if (errtissuecount[isp] > 0)
			if (greensverbose) printf("*** Warning: solute %i, %i tissue source strengths not converged\n", isp, errtissuecount[isp]);
	tissueconv:;

		ofp1 = fopen("Current/GreensLog.txt", "a");		//Print log file
		kvessel = IMIN(kvessel, nmaxvessel);
		ktissue = IMIN(ktissue, nmaxtissue);
		fprintf(ofp1, "\n----- kmain = %i, kvessel = %i, ktissue = %i -----\n", kmain, kvessel, ktissue);
		for (isp = 1; isp <= nsp; isp++) {
			if (diffsolute[isp] == 1) fprintf(ofp1, "Solute %i: qtsum = %f, qvsum = %f, g0 = %f\n",
				isp, qtsum[isp], qvsum[isp] * qvfac[isp], g0[isp]);
			if (permsolute[isp] == 1) fprintf(ofp1, "Solute %i: errvessel_q = %f, imaxerr = %i\n",
				isp, errvessel[isp], segname[imaxerrvessel[isp]]);
			if (diffsolute[isp] == 1) fprintf(ofp1, "Solute %i: errtissue_q = %f, imaxerr = %i\n",
				isp, errtissue[isp], imaxerrtissue[isp]);
		}
		fclose(ofp1);
		//********************** end of tissue loop *****************************
		//Update g0.  If permsolute[isp] != 1, always use method 2.
		//Method 2 is based on derivative wrt g0, automatic estimation of g0fac
		for (isp = 1; isp <= nsp; isp++) g0facnew[isp] = 0.;
		for (itp = 1; itp <= nnt; itp++) {
			for (isp = 1; isp <= nsp; isp++) ptpt[isp] = pt[itp][isp];
			tissrate(nsp, ptpt, mtiss, mptiss);
			ix = tisspoints[1][itp];
			iy = tisspoints[2][itp];
			iz = tisspoints[3][itp];
			for (jtp = 1; jtp <= nnt; jtp++) {
				jx = tisspoints[1][jtp];
				jy = tisspoints[2][jtp];
				jz = tisspoints[3][jtp];
				ixdiff = abs(ix - jx) + 1;
				iydiff = abs(iy - jy) + 1;
				izdiff = abs(iz - jz) + 1;
				for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1) g0facnew[isp] += dtt[ixdiff][iydiff][izdiff] / diff[isp] * mptiss[isp] * vol;
			}
		}
		for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1 && (g0method == 2 || permsolute[isp] == 0)) {
			g0facnew[isp] = 1. / (1. - g0facnew[isp] / nnt);
			dqsumdg0 = FMIN(dqvsumdg0[isp], 0.) + FMIN(dqtsumdg0[isp], 0.)*g0facnew[isp];
			if (fabs(dqsumdg0) > 1.e-6) {
				dif = (qvsum[isp] + qtsum[isp]) / dqsumdg0 * g0fac[isp];//This g0fac should normally be 1.0.
				g0[isp] -= dif;
			}
		}
		convflag = 1;		//Convergence based on changes in pv, pt and g0.  Express relative to eps[isp].
		if (greensverbose) printf("\n");
		for (isp = 1; isp <= nsp; isp++) {
			err = 0.;
			imaxerr = 0;
			if (permsolute[isp] == 1) for (i = 1; i <= nnv; i++) {
				dif = fabs(pv[i][isp] - pvprev[i][isp]) / eps[isp];
				if (dif > err) {
					imaxerr = mainseg[i];
					err = dif;
				}
			}
			errvessel[isp] = err;
			imaxerrvessel[isp] = imaxerr;
			err = 0.;
			imaxerr = 0;
			for (itp = 1; itp <= nnt; itp++) {
				dif = fabs(pt[itp][isp] - ptprev[itp][isp]) / eps[isp];
				if (dif > err) {
					imaxerr = itp;
					err = dif;
				}
			}
			errtissue[isp] = err;
			imaxerrtissue[isp] = imaxerr;
			if (errvessel[isp] > err) {
				imaxerr = imaxerrvessel[isp];
				err = errvessel[isp];
			}
			else imaxerr = -imaxerr;
			dif = fabs(g0[isp] - g0old[isp]) / eps[isp];
			if (dif > err) {
				imaxerr = 0;
				err = dif;
			}
			if (greensverbose) printf("Solute %i: err = %f, imaxerr = %i (- for tissue point)\n", isp, err, imaxerr);
			if (greensverbose && imaxerr > 0) if (lowflow[imaxerr]) printf("Solute %i: max error is at a low-flow segment\n", isp);
			if (err > 1.) convflag = 0;
		}
		ofp1 = fopen("Current/GreensLog.txt", "a");		//Print log file
		for (isp = 1; isp <= nsp; isp++) {
			if (permsolute[isp] == 1) fprintf(ofp1, "Solute %i: errvessel_p = %f, imaxerr = %i\n",
				isp, errvessel[isp], segname[imaxerrvessel[isp]]);
			fprintf(ofp1, "Solute %i: errtissue_p = %f, imaxerr = %i\n", isp, errtissue[isp], imaxerrtissue[isp]);
		}
		fclose(ofp1);
		tfinish1 = clock();
		duration = (float)(tfinish1 - tstart1) / CLOCKS_PER_SEC;
		if (greensverbose) printf("\nkmain = %i, %2.3f seconds for step\n", kmain, duration);
		if (convflag && convflagv && convflagt) goto mainconv;
	}
	printf("\n*** Warning: tissue or vessel solute levels not converged");
mainconv:;
	//********************** end of main loop *****************************
	tfinish = clock();
	duration = (float)(tfinish - tstart) / CLOCKS_PER_SEC;
	printf("\n%i iterations, %2.1f seconds for main loop\n", kmain, duration);
	ofp1 = fopen("Current/GreensLog.txt", "a");
	fprintf(ofp1, "\n%i iterations, %2.1f seconds for main loop\n", kmain, duration);
	fclose(ofp1);

	for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1 && g0method == 1) {
		qvsum[isp] *= qvfac[isp];		//Scale all qv values so that qvsum = qtsum
		for (i = 1; i <= nnv; i++) qv[i][isp] *= qvfac[isp];
	}

	sprintf(fname, "Current/GreensRes%03i.out", imain);		//general output file
	ofp = fopen(fname, "w");
	fprintf(ofp, "%i %i %i %i %i %i\n", nnv, nseg, mxx, myy, mzz, nnt);
	fprintf(ofp, "Scaling factor for flows q0fac = %f\n", q0fac);
	for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "g0[%i] = %f\n", isp, g0[isp]);
	fprintf(ofp, "\n");
	for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1) {	//extravascular solute levels
		fprintf(ofp, "\nSolute %i\n", isp);
		fprintf(ofp, "Segment");
		if (permsolute[isp] == 1) fprintf(ofp, "Efflux Pvessel Ptissue Cvessel");
		fprintf(ofp, "\n");
		for (i = 1; i <= nnv; i++) {
			pev[i][isp] = pvt[i][isp];
			if (g0method == 1 && permsolute[isp] == 1) pev[i][isp] += g0[isp];
			if (permsolute[isp] == 1) for (j = 1; j <= nnv; j++) pev[i][isp] += gvv[i][j] * qv[j][isp] / diff[isp];
		}
		for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
			qvseg[iseg][isp] = 0.;
			pevseg[iseg][isp] = 0.;
			pvseg[iseg][isp] = 0.;
		}
		for (i = 1; i <= nnv; i++) {
			iseg = mainseg[i];
			if (permsolute[isp] == 1) fprintf(ofp, "%4i %4i %10.4f %10.4f %10.4f %10.4f\n",
				i, iseg, qv[i][isp], pv[i][isp], pev[i][isp], cv[i][isp]);
			qvseg[iseg][isp] += qv[i][isp];
			pevseg[iseg][isp] += pev[i][isp] / nspoint[iseg];
			pvseg[iseg][isp] += pv[i][isp] / nspoint[iseg];
		}
		fprintf(ofp, "Solute %i: qtsum = %f, qvsum = %f\n", isp, qtsum[isp], qvsum[isp]);
	}
	for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1) {
		fprintf(ofp, "Solute %i: segment length pvseg pevseg qvseg gamma\n", isp);
		for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5)
			fprintf(ofp, "%4i %10.4f %10.4f %10.4f %10.4f %10.4f\n",
				segname[iseg], lseg[iseg], pvseg[iseg][isp], pevseg[iseg][isp], qvseg[iseg][isp], gamma1[iseg][isp]);
	}
	fclose(ofp);

	sprintf(fname, "Current/VesselLevels%03i.out", imain);	//Vessel levels for all vessel points
	ofp = fopen(fname, "w");
	fprintf(ofp, "Vessel levels\n");
	for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) {
		pmax[isp] = -1.e8;
		pmeanv[isp] = 0.;
		psdv[isp] = 0.;
		pmin[isp] = 1.e8;
		for (i = 1; i <= nnv; i++) {
			pmeanv[isp] += pv[i][isp];
			psdv[isp] += SQR(pv[i][isp]);
			pmax[isp] = FMAX(pv[i][isp], pmax[isp]);
			pmin[isp] = FMIN(pv[i][isp], pmin[isp]);
		}
		pmeanv[isp] = pmeanv[isp]/nnv;
		psdv[isp] = sqrt(psdv[isp]/nnv - SQR(pmeanv[isp]));
		fprintf(ofp, "   Solute %i  ", isp);
	}
	fprintf(ofp,"\npmeanv\n");
	for(isp=1; isp<=nsp; isp++) if(permsolute[isp])	fprintf(ofp,"%12f ", pmeanv[isp]);
	fprintf(ofp,"\npsdv\n");
	for(isp=1; isp<=nsp; isp++) if(permsolute[isp])	fprintf(ofp,"%12f ", psdv[isp]);
	fprintf(ofp, "\npmin\n");
	for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) fprintf(ofp, "%12f ", pmin[isp]);
	fprintf(ofp, "\npmax\n");
	for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) fprintf(ofp, "%12f ", pmax[isp]);
	fprintf(ofp, "\nvalues\n");
	for (i = 1; i <= nnv; i++) {
		for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) fprintf(ofp, "%12f ", pv[i][isp]);
		fprintf(ofp, "\n");
	}
	fclose(ofp);

	sprintf(fname, "Current/VesselSources%03i.out", imain);	//Vessel source strengths for all vessel points
	ofp = fopen(fname, "w");
	fprintf(ofp, "Vessel sources\n");
	for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) fprintf(ofp, "   Solute %i  ", isp);
	fprintf(ofp, "\n");
	for (i = 1; i <= nnv; i++) {
		for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) fprintf(ofp, "%12f ", qv[i][isp]);
		fprintf(ofp, "\n");
	}
	fclose(ofp);

	sprintf(fname, "Current/TissueLevels%03i.out", imain);		//Tissue levels for all tissue points
	ofp = fopen(fname, "w");
	fprintf(ofp, "Tissue levels\n");
	for (isp = 1; isp <= nsp; isp++) {
		pmax[isp] = -1.e8;
		pmeant[isp] = 0.;
		psdt[isp] = 0.;
		pmin[isp] = 1.e8;
		for (itp = 1; itp <= nnt; itp++) {
			pmeant[isp] += pt[itp][isp];
			psdt[isp] += SQR(pt[itp][isp]);
			pmax[isp] = FMAX(pt[itp][isp], pmax[isp]);
			pmin[isp] = FMIN(pt[itp][isp], pmin[isp]);
		}
		pmeant[isp] = pmeant[isp]/nnt;
		psdt[isp] = sqrt(psdt[isp]/nnt - SQR(pmeant[isp]));
		fprintf(ofp, "   Solute %i  ", isp);
	}
	fprintf(ofp,"\npmeant\n");
	for(isp=1; isp<=nsp; isp++)	fprintf(ofp,"%12f ", pmeant[isp]);
	fprintf(ofp,"\npsdt\n");
	for(isp=1; isp<=nsp; isp++)	fprintf(ofp,"%12f ", psdt[isp]);
	fprintf(ofp, "\npmin\n");
	for (isp = 1; isp <= nsp; isp++)	fprintf(ofp, "%12f ", pmin[isp]);
	fprintf(ofp, "\npmax\n");
	for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "%12f ", pmax[isp]);
	fprintf(ofp, "\nvalues\n");
	for (itp = 1; itp <= nnt; itp++) {
		for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "%12f ", pt[itp][isp]);
		fprintf(ofp, "\n");
	}
	fclose(ofp);

	sprintf(fname, "Current/TissueSources%03i.out", imain);		//Tissue source strengths for all tissue points
	ofp = fopen(fname, "w");
	fprintf(ofp, "Tissue sources\n");
	for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "   Solute %i  ", isp);
	fprintf(ofp, "\n");
	for (itp = 1; itp <= nnt; itp++) {
		for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "%12f ", qt[itp][isp]);
		fprintf(ofp, "\n");
	}
	fclose(ofp);
}