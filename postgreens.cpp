/************************************************************************
postgreens.cpp - analyzes results from greens
Uses parameters from PostGreensParams.dat
Includes problem-specific code from postgreens.cpp.dat
Example of usage: compute survival fraction of cells from drug concentration
TWS, May 2015
**************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

void postgreens(void)
{
	extern int max, nsp, nnt, npostgreensparams, npostgreensout;
	extern float **pt, *dtmin, *postgreensparams, *postgreensout;
	extern char numstr[6];
	char fname[80];
	FILE *ofp;

	int	i, isp, itp;

	strcpy(fname, "Current/PostGreens");
	strcat(fname, numstr);
	strcat(fname, ".out");
	ofp = fopen(fname, "w");
	//**************************************************************
#include "postgreens.cpp.dat"
//**************************************************************
	fclose(ofp);
}
