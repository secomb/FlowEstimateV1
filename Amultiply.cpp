/************************************************************************
Amultiply.cpp
Multiply given vector by symmetric matrix for FlowEst2018
Sparse version, with preconditioning
TWS, August 2018
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void Amultiply(double *input, double *output) {
	extern int nnod, matrixdim;
	extern int *nodtyp, *knowntyp, *nodelambda;
	extern int **nodnod;
	extern double *length_weight, *precond;
	extern double **hmat, **kmat;

	int inod, i, currentnode;

	for (inod = 1; inod <= nnod; inod++) {						//row in d/dp equations and known p equations (knowntyp = 0,1,2,3)
		if (knowntyp[inod] == 0) output[inod] = input[inod];	//known pressure node
		else {
			output[inod] = input[inod];// *hmat[0][inod] * precond[inod] * precond[inod];	//preconditioner chosen to put 1 on diagonal
			if (knowntyp[inod] != 3) {
				output[inod] += input[nodelambda[inod]];// *kmat[0][inod] * precond[inod] * precond[nodelambda[inod]];	//preconditioner chosen to put 1 on diagonal
				output[nodelambda[inod]] = input[inod];// *kmat[0][inod] * precond[nodelambda[inod]] * precond[inod];	//preconditioner chosen to put 1 on diagonal
			}
			for (i = 1; i <= nodtyp[inod]; i++) {
				currentnode = nodnod[i][inod];
				if (knowntyp[currentnode] != 0) {
					output[inod] += input[currentnode] * hmat[i][inod] * precond[inod] * precond[currentnode];
					if (knowntyp[currentnode] != 3) output[inod] += input[nodelambda[currentnode]] * kmat[i][inod] * precond[inod] * precond[nodelambda[currentnode]];
					if (knowntyp[inod] != 3) output[nodelambda[inod]] += input[currentnode] * kmat[i][inod] * precond[nodelambda[inod]] * precond[currentnode];
				}
			}
		}
	}
}
