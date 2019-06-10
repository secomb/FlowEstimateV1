/************************************************************************
bvector.cpp
Calculate RHS of linear system for FlowEstimateV1
Create hmat and kmat, which are needed for Amatrix and Amultiply
Create symmetric preconditioner to put 1 on diagonals of submatrices
Used for Methods 1 and 2
TWS, August 2018
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void bvectorc(double *output)
{
	extern int nnod, matrixdim;
	extern int *nodtyp, *ista, *iend, *knowntyp, *nodelambda;
	extern int **nodseg, **nodnod;
	extern double ktau, kpress, presstarget1;
	extern double *cond, *shearfac, *sheartarget, *nodpress, *length_weight, *nodeinflow, *precond;
	extern double *hfactor1, *hfactor2;
	extern double **hmat, **kmat;

	int i, inod, iseg, currentnode;

	for (inod = 1; inod <= nnod; inod++) {
		hmat[0][inod] = kpress * length_weight[inod];	//hmat includes ktau and kpress contributions
		kmat[0][inod] = 0.;
		for (i = 1; i <= nodtyp[inod]; i++) {
			iseg = nodseg[i][inod];
			hmat[0][inod] += hfactor2[iseg];
			hmat[i][inod] = -hfactor2[iseg];
			kmat[0][inod] += cond[iseg] *FMAX(1., 1000. * ktau);	//helps keep lambda values O(1), for error control
			kmat[i][inod] = -cond[iseg] *FMAX(1., 1000. * ktau);	//flip sign of kmat for better convergence of bicgstab
		}
	}
	for (inod = 1; inod <= nnod; inod++) {	//symmetric preconditioners, chosen to put 1 on diagonals of submatrices
		if (knowntyp[inod] != 0) {
			precond[inod] = 1. / sqrt(hmat[0][inod]);
			if (knowntyp[inod] != 3) precond[nodelambda[inod]] = sqrt(hmat[0][inod]) / kmat[0][inod];
		}
		else precond[inod] = 1.;
	}
	for (inod = 1; inod <= nnod; inod++) {		//row in d/dp equations and known p equations (knowntyp = 0,1,2,3)
		if (knowntyp[inod] == 0) output[inod] = nodpress[inod];		//1 in diagonal element of K, element of b is the pressure
		else {
			output[inod] = kpress * length_weight[inod] * presstarget1;
			if (knowntyp[inod] != 3) output[nodelambda[inod]] = nodeinflow[inod] *FMAX(1., 1000. * ktau);
			for (i = 1; i <= nodtyp[inod]; i++) {
				iseg = nodseg[i][inod];
				if (ista[iseg] == inod) output[inod] += hfactor1[iseg];
				if (iend[iseg] == inod) output[inod] -= hfactor1[iseg];
				currentnode = nodnod[i][inod];
				if (knowntyp[currentnode] == 0) {	//contributions to RHS from known pressure nodes (moved out of A matrix to make it symmetric)
					output[inod] -= hmat[i][inod] * nodpress[currentnode];
					if (knowntyp[inod] != 3) output[nodelambda[inod]] -= kmat[i][inod] * nodpress[currentnode];
				}
			}
		}
	}
	for (i = 1; i <= matrixdim; i++) output[i] = output[i] * precond[i];
}