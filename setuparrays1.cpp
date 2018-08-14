/************************************************************************
setuparrays1 - for Greens.  TWS January 08
Set up arrays with dimensions nnod and nseg
Version 2.0, May 1, 2010.
Version 3.0, May 17, 2011.
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setuparrays1(int nseg, int nnod)
{
	extern int nsp, nodsegm;
	extern int *nspoint, *istart, *lowflow, *nk, *nodrank, *nodout, *nodtyp, *ista, *iend;
	extern int **nodnod, **nodseg;
	extern float *segvar, *q, *qq, *oxflux, *rseg, *cbar, *lseg, *ds, *nodvar, *segc;//added November 2016
	extern float **gamma1, **qvseg, **pvseg, **pevseg, **start, **end, **scos;
	extern float ***rsta, ***rend;

	nspoint = ivector(1, nseg);
	istart = ivector(1, nseg);
	lowflow = ivector(1, nseg);
	nk = ivector(1, nnod); //not nseg - error fixed 20 April 2010
	ista = ivector(1, nseg);
	iend = ivector(1, nseg);
	nodrank = ivector(1, nnod);
	nodout = ivector(1, nnod);
	nodtyp = ivector(1, nnod);

	nodnod = imatrix(1, nodsegm, 1, nnod);
	nodseg = imatrix(1, nodsegm, 1, nnod);

	segvar = vector(1, nseg);
	q = vector(1, nseg);//added Novamber 2016
	qq = vector(1, nseg);
	oxflux = vector(1, nodsegm);
	rseg = vector(1, nseg);
	cbar = vector(1, nseg);
	lseg = vector(1, nseg);
	ds = vector(1, nseg);
	nodvar = vector(1, nnod);
	segc = vector(1, nseg);

	gamma1 = matrix(1, nseg, 1, nsp);
	qvseg = matrix(1, nseg, 1, nsp);
	pvseg = matrix(1, nseg, 1, nsp);
	pevseg = matrix(1, nseg, 1, nsp);
	start = matrix(1, 3, 1, nseg);
	scos = matrix(1, 3, 1, nseg);
	end = matrix(1, 3, 1, nseg);
	rsta = f3tensor(1, 3, 1, 16, 1, nseg);
	rend = f3tensor(1, 3, 1, 16, 1, nseg);
}