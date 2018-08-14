/*****************************************************
testconvect - test that convect gives correct results for alpha matrix
Compare matrix values with values obtained by numerical differentiation
TWS August 2010
******************************************************/
#include <math.h>
#include "nrutil.h"
#include <stdio.h>
void convect(int isp);

void testconvect(int isp)
{
	extern int nnv;
	extern int *mainseg, *istart, *nspoint, *nodout, **nodseg, *segname;
	extern float *q, *qq, **al, **cv, **qv, flowfac;

	int i, j;
	float *conv, altest, delq;
	conv = vector(1, nnv);

	delq = qv[1][isp] / 5.;	//small change used for numerical derivatives
	for (j = 1; j <= nnv; j++) {
		qv[j][isp] -= delq;
		convect(isp);
		for (i = 1; i <= nnv; i++) conv[i] = cv[i][isp];
		qv[j][isp] += 2.*delq;
		convect(isp);
		for (i = 1; i <= nnv; i++) conv[i] -= cv[i][isp];
		qv[j][isp] -= delq;	//restore original values
		convect(isp);	//restore original values
		for (i = 1; i <= nnv; i++) {
			altest = conv[i] * qq[mainseg[i]] * flowfac / 2. / delq;
			if (fabs(altest - al[i][j]) > 0.02)
				printf("*** Warning:  Possible error in al matrix %i %i %f %f %i %i\n",
					i, j, al[i][j], altest, segname[mainseg[i]], segname[mainseg[j]]);
		}
	}
	free_vector(conv, 1, nnv);
}