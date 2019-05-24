/************************************************************************
Main program to call greens
Version 2.0, May 1, 2010.
Version 3.0, May 17, 2011.
Version 4.0, March 1, 2018.
See greens.cpp for description of changes.
***********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#if defined(__linux__)
	// Requires c++17 support, should be included in all current linux releases
	#include <experimental/filesystem> 
	namespace fs = std::experimental::filesystem::v1;
#elif defined(__APPLE__)
	// Requires removal of the -lstdc++fs flag from makefile
	#include <filesystem>
	namespace fs = std::filesystem;
#elif defined(_WIN32)    //Windows version
	#include <Windows.h>
#endif

void input(void);
void analyzenet(void);
void picturenetwork(float *nodvar, float *segvar, const char fname[]);
void greens(void);
void contour(const char fname[]);
void histogram(const char fname[]);
void setuparrays0();
void setuparrays1(int nseg, int nnod);
void setuparrays2(int nnv, int nnt);
void cmgui(float *segvar);
void postgreens(void);

int max = 100, nmaxvessel, nmaxtissue, nmax, rungreens, initgreens, g0method, linmethod;
int mxx, myy, mzz, nnt, nseg, nnod, nnodfl, nnv, nsp, nnodbc, nodsegm, nsegfl, kmain;
int slsegdiv, nsl1, nsl2;
int is2d; //needed for 2d version
int nvaryparams, nruns, ntissparams, npostgreensparams, npostgreensout;	//needed for varying parameters, postgreens

float *dtmin;//added July 2011
int *mainseg, *permsolute, *nodrank, *nodtyp, *nodout, *bcnodname, *bcnod, *bctyp, *lowflow;
int *nodname, *segname, *segtyp, *nspoint, *istart, *nl, *nk, *indx, *ista, *iend;
int *errvesselcount, *errtissuecount;
int *imaxerrvessel, *imaxerrtissue, *nresis;  //added April 2010
int *oxygen, *diffsolute; //added April 2010
int **segnodname, **nodseg, **tisspoints, **nodnod;
int ***nbou;

int **tissfix;	//added September 2010
float **tisserr, **dmtissdp, *mptissref;//September 2010;
int **ivaryparams;	//added April 2015

float gtt;	//added September 2010
float fn, c, alphab, p50, cs, cext, hext, req, q0fac, totalq, flowfac = 1.e6 / 60.;
float plow, phigh, clowfac, chighfac, pphighfac;//added January 2012
float pi1 = atan(1.)*4., fac = 1. / 4. / pi1;
float lb, maxl, v, vol, vdom, errfac, tlength, alx, aly, alz, lowflowcrit;
float tlengthq, tlengthqhd;//added 8/09
float xmax, ymax, scalefac;
float w2d, r2d; //needed for 2d version

float *axt, *ayt, *azt, *ds, *diff, *pmin, *pmax, *pmean, *pref, *g0, *g0fac, *g0facnew, *sumal;
float *diam, *rseg, *q, *qdata, *qq, *hd, *oxflux, *segc, *bcprfl, *bchd, *nodvar, *segvar, *qvtemp, *qvfac;//added qdata November 2016
float **start, **scos, **ax, **cnode, **resisdiam, **resis, **bcp; //added April 2010
float **qv, **qt, **pv, **pev, **pt;
float **qvseg, **pvseg, **pevseg;
float **paramvalue, *solutefac, *intravascfac, *postgreensparams, *postgreensout;	//added April 2015

float *x, *y, *lseg, *ss, *cbar, *mtiss, *mptiss, *dqvsumdg0, *dqtsumdg0;
float *epsvessel, *epstissue, *eps, *errvessel, *errtissue, *pinit, *p;
float *rhs, *rhstest, *g0old, *ptt, *ptpt, *qtsum, *qvsum;
float **pvt, **pvprev, **qvprev, **cv, **dcdp, **tissparam;
float **ptprev, **ptv, **gamma1, **qcoeff1, **cv0, **conv0;
float **gvv, **end, **al;
float ***rsta, ***rend, ***dtt;
float *xsl0, *xsl1, *xsl2, *clmin, *clint, *cl, **zv, ***psl;
float **qtp;
double *rhstiss, *matxtiss;
double **mat, **rhsg, *rhsl, *matx;

char numstr[6];

int main(int argc, char *argv[])
{
	int iseg, inod, imain, j, isp;
	char fname[80];

	FILE *ofp;

	// Create current subdirectory if doesn't exist; copy .dat files there
	// Updated May 2019 for cross-platform support
	#if defined(__unix__)
		if (!fs::exists("Current")) fs::create_directory("Current");	

		fs::copy_file("ContourParams.dat", fs::path("Current/ContourParams.dat"), fs::copy_options::overwrite_existing);
		fs::copy_file("SoluteParams.dat", fs::path("Current/SoluteParams.dat"), fs::copy_options::overwrite_existing);
		fs::copy_file("Network.dat", fs::path("Current/Network.dat"), fs::copy_options::overwrite_existing);
		fs::copy_file("IntravascRes.dat", fs::path("Current/IntravascRes.dat"), fs::copy_options::overwrite_existing);
		if (fs::exists("Varyparams.dat")) 
			fs::copy_file("VaryParams.dat", fs::path("Current/VaryParams.dat"), fs::copy_options::overwrite_existing);
		fs::copy_file("tissrate.cpp.dat", fs::path("Current/tissrate.cpp.dat"), fs::copy_options::overwrite_existing);
		
	#elif defined(_WIN32)
		BOOL NoOverwrite = FALSE;
		FILE *ofp;

		CopyFile("SoluteParams.dat","Current/SoluteParams.dat",NoOverwrite);
		CopyFile("IntravascRes.dat","Current/IntravascRes.dat",NoOverwrite);
		CopyFile("ContourParams.dat","Current/ContourParams.dat",NoOverwrite);
		CopyFile("VaryParams.dat","Current/VaryParams.dat",NoOverwrite);
		CopyFile("network.dat","Current/Network.dat",NoOverwrite);
		CopyFile("tissrate.cpp.dat","Current/tissrate.cpp.dat",NoOverwrite);
	#endif

	input();

	is2d = 0; //set to 1 for 2d version, 0 otherwise
	if (mzz == 1) is2d = 1; //assumes 2d version if all tissue points lie in one z-plane

	setuparrays0();

	setuparrays1(nseg, nnod);

	analyzenet();

	setuparrays2(nnv, nnt);

	for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = segname[iseg];
	for (inod = 1; inod <= nnod; inod++) nodvar[inod] = nodname[inod];
	picturenetwork(nodvar, segvar, "Current/NetNodesSegs.ps");
	//for(iseg=1; iseg<=nseg; iseg++) segvar[iseg] = fabs(diam[iseg]);
	for (iseg = 1; iseg <= nseg; iseg++)
		segvar[iseg] = log(fabs(qdata[iseg]));
	cmgui(segvar);

	ofp = fopen("Current/summary.out", "w");
	//print headings for summary output file
	fprintf(ofp, "imain kmain ");
	for (j = 1; j <= nvaryparams; j++) {
		switch (ivaryparams[j][1]) {
		case 1:
		{
			fprintf(ofp, "   q0fac    ");
			break;
		}
		case 2:
		{
			fprintf(ofp, " solutefac[%i]", ivaryparams[j][2]);
			break;
		}
		case 3:
		{
			fprintf(ofp, " diff[%i]     ", ivaryparams[j][2]);
			break;
		}
		case 4:
		{
			fprintf(ofp, " intravascfac[%i]", ivaryparams[j][2]);
			break;
		}
		case 5:
		{
			fprintf(ofp, " tissparam[%i][%i]", ivaryparams[j][2], ivaryparams[j][3]);
			break;
		}
		case 6:
		{
			fprintf(ofp, "   p50     ");
			break;
		}
		}
	}
	for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "  pmean[%i]  ", isp);
	for (j = 1; j <= npostgreensout; j++) fprintf(ofp, " postgreens[%i]", j);
	fprintf(ofp, "\n");

	//The following loop allows running a series of cases with varying parameters
	for (imain = 1; imain <= nruns; imain++) {
		sprintf(numstr, "%03i", imain);	//need 3-digit frame number for file name. November 2016
		for (j = 1; j <= nvaryparams; j++) {
			switch (ivaryparams[j][1])
			{
			case 1:
			{
				q0fac = paramvalue[imain][j];
				break;
			}
			case 2:
			{
				isp = ivaryparams[j][2];	//updated November 2016
				if (isp <= nsp) solutefac[isp] = paramvalue[imain][j];
				break;
			}
			case 3:
			{
				isp = ivaryparams[j][2];
				if (isp <= nsp) diff[isp] = paramvalue[imain][j];
				break;
			}
			case 4:
			{
				isp = ivaryparams[j][2];
				if (isp <= nsp) intravascfac[isp] = paramvalue[imain][j];
				break;
			}
			case 5:
			{
				isp = ivaryparams[j][3];
				if (isp <= nsp) tissparam[ivaryparams[j][2]][isp] = paramvalue[imain][j];
				break;
			}
			case 6:
			{
				p50 = paramvalue[imain][j];
				break;
			}
			}
		}

		greens();

		fprintf(ofp, "%4i  %4i  ", imain, kmain);
		for (j = 1; j <= nvaryparams; j++) fprintf(ofp, "%12f ", paramvalue[imain][j]);
		for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "%12f ", pmean[isp]);

		if (npostgreensparams) postgreens();

		if (npostgreensout) for (j = 1; j <= npostgreensout; j++) fprintf(ofp, "%12f ", postgreensout[j]);
		fprintf(ofp, "\n");

		for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = pvseg[iseg][1];
		for (inod = 1; inod <= nnod; inod++) nodvar[inod] = nodname[inod];

		strcpy(fname, "Current/NetNodesOxygen");
		strcat(fname, numstr);
		strcat(fname, ".ps");
		picturenetwork(nodvar, segvar, fname);

		cmgui(segvar);

		strcpy(fname, "Current/Contour");
		strcat(fname, numstr);
		strcat(fname, ".ps");
		contour(fname);

		strcpy(fname, "Current/Histogram");
		strcat(fname, numstr);
		strcat(fname, ".out");
		histogram(fname);

	}
	fclose(ofp);
	return 0;
}