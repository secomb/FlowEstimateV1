/************************************************************************
putrank - generate list of nodes in order of flow direction
needed for dishem
Considers only type 4 and 5 segments
nodrank --- if i < j, node nodrank[j] is not upstream of node nodrank[i]
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void putrank(void)
{
	extern int nseg, nnod, nnodfl, nsegfl;
	extern int *nodrank, *nodtyp, *nodout, *nodname, *segtyp, *nk, *ista, *iend;
	extern int **nodseg, **nodnod;
	extern float *q;

	int inod, j, iseg, nod1, nod2, flag;

	//reconstruct node table; count outputs from node; output nodes precede input nodes
	for (inod = 1; inod <= nnod; inod++) {
		nodtyp[inod] = 0;
		nodout[inod] = 0;
	}
	nsegfl = 0;	//added TWS 2010
	for (iseg = 1; iseg <= nseg; iseg++) {//segtypes 4 and 5 are flowing segments; all segments in this network are flowing
		if (q[iseg] >= 0) {
			nod1 = ista[iseg];
			nod2 = iend[iseg];
		}
		//If q>=0, flow goes from ista to iend. If q<0, then flow goes from iend to ista. This if-else statement ensures that nod1 is always the supplier of blood flow. -B
		else {
			nod1 = iend[iseg];
			nod2 = ista[iseg];
		}
		nodtyp[nod1]++;
		nodseg[nodtyp[nod1]][nod1] = iseg;
		nodnod[nodtyp[nod1]][nod1] = nod2;
		nodout[nod1]++;//nodout is the number of segments for which blood flows from the node to the segment
		nsegfl++;
	}
	for (iseg = 1; iseg <= nseg; iseg++) { //two for loops guarantees that the outflow nodes come first
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
	//assign low ranks to inflow nodes
	nnodfl = 0;
	for (inod = 1; inod <= nnod; inod++) {
		nk[inod] = 0;//If node has been ranked, nk=0. Otherwise, nk=1.
		if (nodtyp[inod] == 1 && nodout[inod] == 1) {//node is connected and provides flow to only one segment
			nnodfl++;
			nk[inod] = 1;
			nodrank[nnodfl] = inod;//indices of nodrank are the ranks, elements are nod indices
		}
	}
	//assign increasing ranks to downstream connected nodes
	flag = 1;
	while (flag == 1) {
		flag = 0;
		for (inod = 1; inod <= nnod; inod++)	if (nk[inod] == 0 && nodtyp[inod] > 0) {
			for (j = nodout[inod] + 1; j <= nodtyp[inod]; j++) {
				iseg = nodseg[j][inod];		//an inflow segment
				if (inod == iend[iseg] && (nk[ista[iseg]] == 0 || q[iseg] < 0.)) goto skipnode;
				if (inod == ista[iseg] && (nk[iend[iseg]] == 0 || q[iseg] >= 0.)) goto skipnode;
			}
			nnodfl++;
			nk[inod] = 1;
			nodrank[nnodfl] = inod;
			flag = 1;
		skipnode:;
		}
	}
	//check for unprocessed nodes--should be none.
	for (inod = 1; inod <= nnod; inod++) if (nodtyp[inod] != 0 && nk[inod] == 0)
		printf("*** Error: unprocessed node %i in putrank\n", inod);
}
