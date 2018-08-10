/************************************************************************
Amatrix.cpp
Create full symmetric matrix for FlowEst10, with preconditioning
TWS, August 2018
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void Amatrix(double **matrix) {
	extern int nnod, matrixdim;
	extern int *nodtyp, *knowntyp, *nodelambda;
	extern int **nodnod;
	extern double *length_weight, *precond;
	extern double **hmat, **kmat;

	int inod, i, j, currentnode;

	for (i = 1; i <= matrixdim; i++) for (j = 1; j <= matrixdim; j++) matrix[i][j] = 0.;
	for (inod = 1; inod <= nnod; inod++) {
		if (knowntyp[inod] == 0) matrix[inod][inod] = 1.;	//known pressure node
		else {
			matrix[inod][inod] = hmat[0][inod] * precond[inod] * precond[inod];
			if (knowntyp[inod] != 3) {
				matrix[inod][nodelambda[inod]] = kmat[0][inod] * precond[inod] * precond[nodelambda[inod]];
				matrix[nodelambda[inod]][inod] = kmat[0][inod] * precond[nodelambda[inod]] * precond[inod];
			}
			for (i = 1; i <= nodtyp[inod]; i++) {
				currentnode = nodnod[i][inod];
				if (knowntyp[currentnode] != 0) {
					matrix[inod][currentnode] = hmat[i][inod] * precond[inod] * precond[currentnode];
					if (knowntyp[currentnode] != 3) matrix[inod][nodelambda[currentnode]] = kmat[i][inod] * precond[inod] * precond[nodelambda[currentnode]];
					if (knowntyp[inod] != 3) matrix[nodelambda[inod]][currentnode] = kmat[i][inod] * precond[nodelambda[inod]] * precond[currentnode];
				}
			}
		}
	}
}
