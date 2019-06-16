/************************************************************************
setuparrays2 - for Greens.  TWS January 08
Set up arrays with dimensions nnv and nnt
Version 2.0, May 1, 2010.
Version 3.0, May 17, 2011.
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setuparrays2(int nnv, int nnt)
{
	extern int nsp, *indx, *mainseg, **tisspoints;
	extern float *rhs, *qvtemp, *sumal, *dtmin;
	extern float **qv, **pv, **pev, **pvt, **pvprev, **qvprev, **cv, **dcdp, **cv0, **conv0, **gvv;
	extern float **qt, **pt, **ptprev, **ptv, **ax, **al;
	extern double **mat, **rhsg, *rhsl, *matx;

	mainseg = ivector(1, nnv);
	indx = ivector(1, nnv + 1);
	tisspoints = imatrix(1, 3, 1, nnt);

	dtmin = vector(1, nnt);
	rhs = vector(1, nnv);
	qvtemp = vector(1, nnv);
	sumal = vector(1, nnv);

	qv = matrix(1, nnv, 1, nsp);
	pv = matrix(1, nnv, 1, nsp);
	pev = matrix(1, nnv, 1, nsp);
	pvt = matrix(1, nnv, 1, nsp);
	pvprev = matrix(1, nnv, 1, nsp);
	qvprev = matrix(1, nnv, 1, nsp);
	cv = matrix(1, nnv, 1, nsp);
	dcdp = matrix(1, nnv, 1, nsp);
	cv0 = matrix(1, nnv, 1, nsp);
	conv0 = matrix(1, nnv, 1, nsp);
	gvv = matrix(1, nnv, 1, nnv);	//this is the large matrix
	qt = matrix(1, nnt, 1, nsp);
	pt = matrix(1, nnt, 1, nsp);
	ptprev = matrix(1, nnt, 1, nsp);
	ptv = matrix(1, nnt, 1, nsp);
	ax = matrix(1, 3, 1, nnv);
	al = matrix(1, nnv, 1, nnv);
	mat = dmatrix(1, nnv + 1, 1, nnv + 1);
	rhsg = dmatrix(1, nnv + 1, 1, 2);
	rhsl = dvector(1, nnv + 1);
	matx = dvector(1, nnv + 1);
}