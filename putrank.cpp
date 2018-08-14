/************************************************************************
putrank - generate list of nodes in order of flow direction
nodrank --- if nodrank[i] < nodrank[j], node j is not upstream of node i
Version 2.0, May 1, 2010.
Version 3.0, May 17, 2011.
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void putrank(void)
{
	extern int nseg, nnod, nnodfl, nnodbc, nsegfl;
	extern int *nodrank, *nodtyp, *nodout, *segtyp, *nk, *ista, *iend;
	extern int **nodseg, **nodnod;
	extern float *q;

	int inod, j, jseg, nod1, nod2, iseg, flag;

	//construct node table; count outputs from node; output nodes precede input nodes
	for (inod = 1; inod <= nnod; inod++) {
		nodtyp[inod] = 0;
		nodout[inod] = 0;
	}
	nsegfl = 0;	//added TWS 2010
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
		if (q[iseg] >= 0) {
			nod1 = ista[iseg];
			nod2 = iend[iseg];
		}
		else {
			nod1 = iend[iseg];
			nod2 = ista[iseg];
		}
		nodtyp[nod1]++;
		nodseg[nodtyp[nod1]][nod1] = iseg;
		nodnod[nodtyp[nod1]][nod1] = nod2;
		nodout[nod1]++;
		nsegfl++;
	}
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
		if (q[iseg] >= 0) {
			nod1 = ista[iseg];
			nod2 = iend[iseg];
		}
		else {
			nod1 = iend[iseg];
			nod2 = ista[iseg];
		}
		nodtyp[nod2]++;
		nodseg[nodtyp[nod2]][nod2] = iseg;
		nodnod[nodtyp[nod2]][nod2] = nod1;
	}
	for (inod = 1; inod <= nnod; inod++)	if (nodtyp[inod] == 0) printf("***Warning: Node %i is not related to any segment\n", inod);
	//assign low ranks to inflow nodes
	nnodfl = 0;
	for (inod = 1; inod <= nnod; inod++) {
		nk[inod] = 0;
		if (nodtyp[inod] == 1 && nodout[inod] == 1) {
			nnodfl++;
			nk[inod] = 1;
			nodrank[nnodfl] = inod;
		}
	}
	//assign increasing ranks to downstream connected nodes
	flag = 1;
	while (flag == 1) {
		flag = 0;
		for (inod = 1; inod <= nnod; inod++)	if (nk[inod] == 0 && nodtyp[inod] > 0) {
			for (j = nodout[inod] + 1; j <= nodtyp[inod]; j++) {
				jseg = nodseg[j][inod];
				if (inod == iend[jseg] && (nk[ista[jseg]] == 0 || q[jseg] <= 0.)) goto skipnode;
				if (inod == ista[jseg] && (nk[iend[jseg]] == 0 || q[jseg] >= 0.)) goto skipnode;
			}
			nnodfl++;
			nk[inod] = 1;
			nodrank[nnodfl] = inod;
			flag = 1;
		skipnode:;
		}
	}
}