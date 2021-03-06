/************************************************************************
analyzenet - for flowest.  TWS August 09
Modified for FlowEstimateV1, June 2019
Set up nodtyp, nodseg, nodnod arrays based on flowing segments
ista[i]: the node where the ith segment starts
iend[i]: the node where the ith segment ends
nodtyp[i]: the number of segments connected to the ith node
************************************************************************
knowntyp	0 for known pressures - Nkp
			1 for interior nodes
			2 for known flow boundary conditions
			3 for unknown boundary conditions - Nup
***********************************************************************/ 
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void analyzenet()
{
	extern int nseg, nnod, nnodbc, nnodbck, nodsegm, solvetyp;
	extern int numberknownpress, numberunknown, matrixdim;
	extern int *bcnodname, *bcnod, *bctyp, *nodtyp, *nodname, *ista, *iend, *knowntyp, *nodelambda, *idx;
	extern int *known_flow_direction, *segtyp;
	extern int **segnodname, **nodnod, **nodseg;
	extern float pi1, totallength, consthd;
	extern float *diam, *q, *qq, *hd, *bcprfl, *bchd, **cnode;
	extern double presstarget1, *lseg, *length_weight;
	extern double *nodpress, *lambda, *nodeinflow, *dd;
	extern double **hmat, **kmat;
	extern double **fullmatrix, *bvector, *xvector, *precond;

	int i, k, iseg, inod, inod1, inod2, inodbc, counter, nnodtyp0;
	double totallength1;

	for (iseg = 1; iseg <= nseg; iseg++) {	//Search for nodes corresponding to this segment
		for (inod = 1; inod <= nnod; inod++) if (nodname[inod] == segnodname[1][iseg]) {
			ista[iseg] = inod;
			goto foundit1;
		}
		printf("*** Error: No matching node found for nodname %i\n", segnodname[1][iseg]);
	foundit1:;
		for (inod = 1; inod <= nnod; inod++) if (nodname[inod] == segnodname[2][iseg]) {
			iend[iseg] = inod;
			goto foundit2;
		}
		printf("*** Error: No matching node found for nodname %i\n", segnodname[2][iseg]);
	foundit2:;
	}
	//Setup nodtyp, nodseg and nodnod
	for (inod = 1; inod <= nnod; inod++) nodtyp[inod] = 0;
	for (iseg = 1; iseg <= nseg; iseg++) {
		inod1 = ista[iseg]; 
		inod2 = iend[iseg]; 
		nodtyp[inod1]++; 
		nodtyp[inod2]++; 
		if (nodtyp[inod1] > nodsegm) printf("*** Error: Too many segments connected to node %i\n", inod1);
		if (nodtyp[inod2] > nodsegm) printf("*** Error: Too many segments connected to node %i\n", inod2);
		nodseg[nodtyp[inod1]][inod1] = iseg; 
		nodseg[nodtyp[inod2]][inod2] = iseg; 
		nodnod[nodtyp[inod1]][inod1] = inod2; 
		nodnod[nodtyp[inod2]][inod2] = inod1; 
	}
	//segment lengths
	totallength = 0.;
	for (iseg = 1; iseg <= nseg; iseg++)  {	// calculate the length of each segment by Pythagoras
		lseg[iseg] = 0;
		for (k = 1; k <= 3; k++) lseg[iseg] += SQR(cnode[k][iend[iseg]] - cnode[k][ista[iseg]]);
		lseg[iseg] = sqrt(lseg[iseg]);
		totallength += lseg[iseg];
		if (lseg[iseg] < 1.e-3) printf("*** Error: Segment %i has length %f\n", iseg, lseg[iseg]);
	}
	//node weighting factors: 1/2 the sum of the lengths of the segments connected to the node
	totallength1 = 0.;
	for (inod = 1; inod <= nnod; inod++) {
		length_weight[inod] = 0;
		for (i = 1; i <= nodtyp[inod]; i++) {
			iseg = nodseg[i][inod];
			length_weight[inod] += 0.5*lseg[iseg];
		}
		totallength1 += length_weight[inod];
	}
	if (fabs(totallength - totallength1) / totallength > 1.e-5) printf("*** Error in nodal length weights\n");
	// Classify node types
	numberknownpress = 0;
	numberunknown = 0;
	nnodbc = 0;
	nnodtyp0 = 0;
	for (inod = 1; inod <= nnod; inod++) {	//initialize all nodes
		if (nodtyp[inod] == 0) {
			nnodtyp0++;
			if(nnodtyp0 == 1) printf("*** Warning: Nodes with no attached flowing segment\n");
			printf("%i ", inod);
			if(nnodtyp0 % 20 == 0) printf("\n");
		}
		if (nodtyp[inod] == 1) {			//identify boundary nodes
			nnodbc++;
			bcnod[nnodbc] = inod;			//set the inodbc-th boundary node to node inod
			bchd[nnodbc] = consthd;			//boundary hematocrits are unknown
			knowntyp[inod] = 3;				//initially assume all boundary nodes are unknown
			numberunknown++;
		}
		else knowntyp[inod] = 1;			//internal nodes are type 1
		nodeinflow[inod] = 0.;
		lambda[inod] = 0.;
		//set pressures on first iteration with some randomness to avoid zero flow segments
		nodpress[inod] = presstarget1 + rand()*10. / RAND_MAX - 5.;		
	}
	if (nnodtyp0 % 20 != 0) printf("\n");
	for (iseg = 1; iseg <= nseg; iseg++) hd[iseg] = consthd;		//initialize hematocrits
	//boundary nodes with known boundary conditions or flow
	for (iseg = 1; iseg <= nseg; iseg++) known_flow_direction[iseg] = 0;
	for (inodbc = 1; inodbc <= nnodbck; inodbc++) { 		//Search for node corresponding to this node name
		for (inod = 1; inod <= nnod; inod++) if (nodname[inod] == bcnodname[inodbc]) {
			if (nodtyp[inod] != 1) printf("*** Error: Known boundary node %i is not a 1-segment node, nodtyp = %i\n", inod, nodtyp[inod]);
			if (bctyp[inodbc] == 0) {
				knowntyp[inod] = 0;
				nodpress[inod] = bcprfl[inodbc];
				numberknownpress++;
				numberunknown--;
			}
			if (bctyp[inodbc] == 2) {
				knowntyp[inod] = 2;
				nodeinflow[inod] = bcprfl[inodbc];
				numberunknown--;
			}
			if (bctyp[inodbc] == -2) {		//known flow direction only
				knowntyp[inod] = 3;
				iseg = nodseg[1][inod];
				if(bcprfl[inodbc] > 0. && inod == ista[iseg]) known_flow_direction[iseg] = 1;
				else if (bcprfl[inodbc] < 0. && inod == iend[iseg]) known_flow_direction[iseg] = 1;
				else known_flow_direction[iseg] = -1;
			}
			goto foundit;
		}
		printf("*** Error: No matching node found for nodname %i\n", bcnodname[inodbc]);
	foundit:;
	}
	matrixdim = 2 * nnod - numberunknown - numberknownpress;		//size of matrix
	counter = nnod;
	for (inod = 1; inod <= nnod; inod++) if (knowntyp[inod] == 1 || knowntyp[inod] == 2) {	//has an associated constraint
		counter++;
		nodelambda[inod] = counter;			//position of lambda corresponding corresponding to inod in vector of unknowns
	}
	if (counter != 2 * nnod - numberunknown - numberknownpress) printf("*** Error: incorrect number of constraints\n");

	//set up vectors and matrices
	bvector = dvector(1, matrixdim);
	xvector = dvector(1, matrixdim);
	precond = dvector(1, matrixdim);
	if (solvetyp == 1) {	//not needed for sparse version
		idx = ivector(1, matrixdim);
		dd = dvector(1, matrixdim);
		fullmatrix = dmatrix(1, matrixdim, 1, matrixdim);
	}
}