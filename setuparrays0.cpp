/************************************************************************
setuparrays0 - for Greens.  TWS January 08
Set up arrays with fixed dimensions
Version 2.0, May 1, 2010.
Version 3.0, May 17, 2011.
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setuparrays0()
{
	extern int nsp, mxx, myy, mzz;
	extern int *errvesselcount, *errtissuecount;
	extern int *imaxerrvessel, *imaxerrtissue;  //added April 2010
	extern int ***nbou;
	extern float *pmin, *pmax, *pmean, *mtiss, *mptiss, *g0old, *ptt, *ptpt, *qtsum, *qvsum, *errvessel;
	extern float *errtissue, *dqvsumdg0, *dqtsumdg0, *pinit, *p, *epsvessel, *epstissue, *eps, *g0, *g0facnew;
	extern float *qvfac; //May 2010
	extern float *x, *y, *ss, *axt, *ayt, *azt;
	extern float ***dtt;
	extern float *mptissref;//September 2010

	errvesselcount = ivector(1, nsp);
	errtissuecount = ivector(1, nsp);
	imaxerrvessel = ivector(1, nsp);
	imaxerrtissue = ivector(1, nsp);
	nbou = i3tensor(1, mxx, 1, myy, 1, mzz);

	pmin = vector(1, nsp);
	pmax = vector(1, nsp);
	pmean = vector(1, nsp);
	mtiss = vector(1, nsp);
	mptiss = vector(1, nsp);
	g0old = vector(1, nsp);
	ptt = vector(1, nsp);
	ptpt = vector(1, nsp);
	qtsum = vector(1, nsp);
	qvsum = vector(1, nsp);
	errvessel = vector(1, nsp);
	errtissue = vector(1, nsp);
	dqvsumdg0 = vector(1, nsp);
	dqtsumdg0 = vector(1, nsp);
	g0facnew = vector(1, nsp);
	pinit = vector(1, nsp);
	p = vector(1, nsp);
	epsvessel = vector(1, nsp);
	epstissue = vector(1, nsp);
	eps = vector(1, nsp);
	qvfac = vector(1, nsp);
	mptissref = vector(1, nsp);

	x = vector(1, 3);
	y = vector(1, 3);
	ss = vector(1, 3);

	axt = vector(1, mxx);
	ayt = vector(1, myy);
	azt = vector(1, mzz);

	dtt = f3tensor(1, mxx, 1, myy, 1, mzz);
}