/*****************************************************
Tissue uptake rates of solutes as a function of levels
Must be provided for each application.  TWS November 07.
Version 1.0, May 1, 2008.
Version 2.0, May 1, 2010.
Version 3.0, May 17, 2011.
***  Must rebuild project after switching tissrate.cpp.dat files ***
******************************************************/
#include <math.h>
#include "nrutil.h"
#include <stdio.h>
void tissrate(int nsp, float *p, float *mtiss, float *mptiss)
{
	extern int *oxygen;
	extern float **tissparam;

	int isp;

#include "tissrate.cpp.dat"
}



/*

Some sample code for tissrate.cpp.dat


	float pcr,m0,gf0,gf1,gf2,aterm;//not external, April 2010
	int isp;
	for(isp=1; isp<=nsp; isp++){
		switch (isp)
		{
		case 1: //oxygen
			m0 = tissparam[1][isp];
			pcr = tissparam[2][isp];
			aterm = 0.;
			if(p[isp] >= 0.){
				mtiss[isp] = -m0*p[isp]/(p[isp] + pcr);
				mptiss[isp] = -m0*pcr/SQR(p[isp] + pcr);
			}
			else{
				mtiss[isp] = 0.;
				mptiss[isp] = 0.;
			}
			break;
/*
		case 2: //VEGF: non-permeable diffusible solute, based on Mac Gabhann and Popel - 2010
			if(p[1] <= 1.) mtiss[2] = 6.*tissparam[1][2];
			else if(p[1] <= 20.) mtiss[2] = (1. + 5.*pow((20. - p[1])/19.,3.))*tissparam[1][2];
			else mtiss[2] = tissparam[1][2];
			mtiss[2] -= tissparam[2][2]*p[2];
			mptiss[2] = -tissparam[2][2];

		case 2: //non-permeable diffusible solute produced in hypoxic regions - old version
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			if(p[1] >= 0.) mtiss[isp] = gf0*pcr/(p[1] + pcr) - gf1*p[isp];
			else mtiss[isp] = gf0 - gf1*p[isp];
			mptiss[isp] = -gf1;
			break;
		case 2: //permeable solute delivered in blood with linear consumption in hypoxic regions
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			gf2 = tissparam[3][isp];
			if(p[1] >= 0.) mtiss[isp] = - gf1*p[isp]*pcr/(p[1] + pcr);
			else mtiss[isp] =  - gf1*p[isp];
			mptiss[isp] = -gf1*pcr/(p[1] + pcr);
			break;

		case 2: //permeable solute delivered in blood with linear consumption in hypoxic regions - KOH, updated TWS July 2011
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			gf2 = tissparam[3][isp];
			if(p[1] >= 0.){
				aterm = gf1*p[isp]*pcr/(p[1] + pcr);
				mtiss[isp] = -aterm;
				mptiss[isp] = -gf1*pcr/(p[1] + pcr);
			}
			else{
				aterm = gf1*p[isp];
				mtiss[isp] = -aterm;
				mptiss[isp] = -gf1;
			}

			break;
		case 3: //permeable diffusiblbe solute produced in hypoxic regions by solute 2
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			gf2 = tissparam[3][isp];
			if(p[1] >= 0.) mtiss[isp] =  aterm;
			else mtiss[isp] =  aterm;
			mptiss[isp] = 0;
			break;
/*
		case 3: //permeable solute delivered in blood with linear consumption in hypoxic regions
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			gf2 = tissparam[3][isp];
			if(p[1] >= 0.) mtiss[isp] = - gf1*p[isp]*pcr/(p[1] + pcr);
			else mtiss[isp] =  - gf1*p[isp];
			mptiss[isp] = -gf1*pcr/(p[1] + pcr);
			break;
		case 3: //permeable solute delivered in blood with linear consumption
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			mtiss[isp] = -gf1*p[isp];
			mptiss[isp] = -gf1;
			break;

		case 4: //non-diffusible solute produced in hypoxic regions
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			if(p[1] >= 0.) mtiss[isp] = gf0*pcr/(p[1] + pcr) - gf1*p[isp];
			else mtiss[isp] = gf0 - gf1*p[isp];
			mptiss[isp] = -gf1;
			break;

		case 4: //non-diffusible non permeable solute produced from case 2 in hypoxic regions
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			gf2 = tissparam[3][isp];
			if(p[1] >= 0.) mtiss[isp] =  aterm - gf0*p[isp];
			else mtiss[isp] =  aterm - gf0*p[isp];
			mptiss[isp] = -gf0;
			break;
		default:
			printf("Error: nsp is too large in tissrate\n");
		}
	}
*/