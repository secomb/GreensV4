/**********************************************************
contour.cpp - generate data for contour plot.  TWS Dec. 07
Version 3.0, May 17, 2011.
Produces a single postscript file with a page for each solute
***********************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void contr_lines(FILE *ofp, int m, int n, float scalefac, int nl,
	float xmin, float xmax, float ymin, float ymax, float *cl, float **zv);
void contr_shade(FILE *ofp, int m, int n, float scalefac, int nl,
	float xmin, float xmax, float ymin, float ymax, float *cl, float **zv,
	int showscale, int lowcolor, int hatch, int plotcontour);
float *eval(int slsegdiv, float req, float *x);

void contour(const char fname[])
{
	extern int max, nsp, nseg, *segtyp, *nl, *ista, *iend;
	extern int slsegdiv, nsl1, nsl2;
	extern float pi1, req, scalefac;
	extern float *x, *p, *diam, **cnode, **pvseg;
	extern float *xsl0, *xsl1, *xsl2, *clmin, *clint, *cl, **zv, ***psl;

	int i, iseg, isp, isl1, isl2, ilevel, nlevel = 100;
	float xmin, ymin, xmax, ymax, xs1, ys1, xs2, ys2, **cos;
	float red, green, blue, xz, xzmin, xzmax;
	float diamfac = 1., zcoord, zbottom, ztop, zmin, zmax;


	FILE *ofp;

	printf("Generating data for contour plots...");

	xmin = 0.;
	xmax = sqrt(SQR(xsl1[1] - xsl0[1]) + SQR(xsl1[2] - xsl0[2]) + SQR(xsl1[3] - xsl0[3]));
	ymin = 0.;
	ymax = sqrt(SQR(xsl2[1] - xsl0[1]) + SQR(xsl2[2] - xsl0[2]) + SQR(xsl2[3] - xsl0[3]));
	cos = matrix(1, 3, 1, 3);
	for (i = 1; i <= 3; i++) {	//set up matrix of direction cosines
		cos[1][i] = (xsl1[i] - xsl0[i]) / xmax;
		cos[2][i] = (xsl2[i] - xsl0[i]) / ymax;
	}
	cos[3][1] = cos[1][2] * cos[2][3] - cos[1][3] * cos[2][2];
	cos[3][2] = cos[1][3] * cos[2][1] - cos[1][1] * cos[2][3];
	cos[3][3] = cos[1][1] * cos[2][2] - cos[1][2] * cos[2][1];

	//Determine range of z values
	zmin = 1.e6;
	zmax = -1.e6;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
		zcoord = 0.;
		for (i = 1; i <= 3; i++)	zcoord += (cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2.*cos[3][i];
		zmin = FMIN(zmin, zcoord - 1.);
		zmax = FMAX(zmax, zcoord + 1.);
	}

	//Calculate P on a planar slice through the region, specified by three corners and number of points along each edge
	//Subdivide vessel segments by slsegdiv for smoother plots
	for (isl1 = 1; isl1 <= nsl1; isl1++)	for (isl2 = 1; isl2 <= nsl2; isl2++) {
		for (i = 1; i <= 3; i++) x[i] = xsl0[i] + (isl1 - 1)*(xsl1[i] - xsl0[i]) / (nsl1 - 1) + (isl2 - 1)*(xsl2[i] - xsl0[i]) / (nsl2 - 1);
		p = eval(slsegdiv, req, x);
		for (isp = 1; isp <= nsp; isp++) psl[isl1][isl2][isp] = p[isp];
	}
	xmin = 0.;
	ymin = 0.;
	ofp = fopen(fname, "w");
	fprintf(ofp, "%%!PS-Adobe-2.0\n");
	fprintf(ofp, "%%%%Pages: %i\n", nsp);
	fprintf(ofp, "%%%%EndComments\n");
	for (isp = 1; isp <= nsp; isp++) {
		fprintf(ofp, "%%%%Page: %i %i\n", isp, isp);
		for (isl1 = 1; isl1 <= nsl1; isl1++)	for (isl2 = 1; isl2 <= nsl2; isl2++) zv[isl1][isl2] = psl[isl1][isl2][isp];
		for (i = 1; i <= nl[isp]; i++) cl[i] = clmin[isp] + (i - 1)*clint[isp];

		if (isp == 1) {
			cl[1] = 1.;		//override contour levels to give contours at p = 1, 2 and 5
			cl[2] = 2.;
			cl[3] = 5.;
		}


		contr_shade(ofp, nsl1, nsl2, scalefac, nl[isp], xmin, xmax, ymin, ymax, cl, zv, 1, 1, 0, 0);
		//		contr_lines(ofp,nsl1,nsl2,scalefac,nl[isp],xmin,xmax,ymin,ymax,cl,zv);

		fprintf(ofp, "/sl {setlinewidth} def\n");
		fprintf(ofp, "/sc {setrgbcolor} def\n");
		fprintf(ofp, "/s {stroke} def\n");
		fprintf(ofp, "1 setlinecap\n");

		//Plot projection of network in contour plane
		//plot vessels according to pvseg in order from bottom to top according to z-coordinate
		xzmin = clmin[isp];
		xzmax = clmin[isp] + (nl[isp] - 1)*clint[isp];
		for (ilevel = 1; ilevel <= nlevel; ilevel++) {
			zbottom = zmin + (ilevel - 1)*(zmax - zmin) / nlevel;
			ztop = zmin + ilevel * (zmax - zmin) / nlevel;
			for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
				zcoord = 0.;
				for (i = 1; i <= 3; i++)	zcoord += (cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2.*cos[3][i];
				if (zcoord >= zbottom && zcoord < ztop) {
					if (xzmin != xzmax) xz = (pvseg[iseg][isp] - xzmin) / (xzmax - xzmin);
					else xz = 0.75;
					blue = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.25), 0.), 1.);//Set up colors using Matlab 'jet' scheme
					green = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.5), 0.), 1.);
					red = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.75), 0.), 1.);

					xs1 = 0.;
					ys1 = 0.;
					xs2 = 0.;
					ys2 = 0.;
					for (i = 1; i <= 3; i++) {
						xs1 += (cnode[i][ista[iseg]] - xsl0[i])*cos[1][i];
						ys1 += (cnode[i][ista[iseg]] - xsl0[i])*cos[2][i];
						xs2 += (cnode[i][iend[iseg]] - xsl0[i])*cos[1][i];
						ys2 += (cnode[i][iend[iseg]] - xsl0[i])*cos[2][i];
					}
					fprintf(ofp, "0 0 0 sc\n");	//Plot vessels slightly larger in black to outline
					fprintf(ofp, "%g sl\n", scalefac*diam[iseg] * diamfac + 2.);
					fprintf(ofp, "%g mx %g my m %g mx %g my l s \n", xs1, ys1, xs2, ys2);
					fprintf(ofp, "%f %f %f sc\n", red, green, blue);
					fprintf(ofp, "%g sl\n", scalefac*diam[iseg] * diamfac);//line widths scaled up by diamfac
					fprintf(ofp, "%g mx %g my m %g mx %g my l s \n", xs1, ys1, xs2, ys2);
				}
			}
		}
		fprintf(ofp, "0 0 0 setrgbcolor\n");		//black
		fprintf(ofp, "/Times-Roman findfont 12 scalefont setfont\n");
		fprintf(ofp, "50 30 moveto\n");
		fprintf(ofp, "(Solute %i) show\n", isp);  //show solute number
//create a scale bar
		float barlength = 100;
		if (xmax > 500.) barlength = 200.;
		if (xmax > 1500.) barlength = 500.;
		fprintf(ofp, "%g sl\n", 2.);
		fprintf(ofp, "%g %g m %g %g l stroke\n", 120., 25., 120. + scalefac * barlength, 25.);
		fprintf(ofp, "%g %g m (%g mm) show\n", 120., 30., barlength / 1000.);
		fprintf(ofp, "showpage\n");
	}
	fclose(ofp);
	printf("done\n");
}