/************************************************************************
flow - for FlowEst2018
Flows in nl/min, pressures in mmHg, viscosity in cP, shear stress in dyn/cm2
Lengths, diameters in microns, times in s
Use double precision pressures for networks with many segments
System setup, based on Fry, Lee, Smith, Secomb. 2012
Version with square submatrices, TWS June 2018
With strictly symmetric matrix and preconditioning, using conjugate gradient solver
TWS August 2018

				N				IUB'
     __                                    __ __  __          __               __
     |                       |              | | p1 |          |                 |
     |                                      | | p2 |          |                 |
     |       1 + ktH/kpw     |     -KT      | | p3 |          |p+kt*ScLMtau/kpw |
  N  |         or 1                         | | .  |          |   or fixed p    |
     |                       |              | | .  |          |                 |
     |                                      | | pn |          |                 |
     |-----------------------+--------------| |----|  _____   |-----------------|
     |                                      | | l1 |  _____   |      q01        |
     |                       |              | | l2 |          |      q02        |
     |                                      | | l3 |          |      q03        |
     |          -K           |        0     | | .  |          |  also includes  |
 IUB'|                                      | | .  |          |   terms from    |
     |                       |              | | .  |          | fixed p values  |
     |                                      | | .  |          |        .        |
     |_                      |             _| |_  _|          |_               _|

N = nnod = total number of nodes
Nkp = numberknownpress = number of boundary nodes with known pressures
Nu  = numberunknown = number of boundary nodes with unknown pressures and flows
IUB' = number of internal nodes + number of boundary nodes with known flow
      = N - Nu - Nkp
IUB' + N = matrixdim = 2*N - Nkp - Nu
***********************************************************************
knowntyp	0 for known pressures - Nkp
			1 for interior nodes
			2 for known flow boundary conditions
			3 for unknown boundary conditions - Nu
***********************************************************************
solvetyp	1 lu decomposition
			2 sparse conjugate gradient
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void bvectorc(double *output);
void Amatrix(double **matrix);
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double *b);
double sparse_cgsymm(double *b, double *x, int n, double eps, int itmax);
void putrank(void);
void dishem(void);
float viscor(float d, float h);

void Amultiply(double *input, double *output);

void flow()
{
	extern int numberknownpress, numberunknown, matrixdim, solvetyp;
	extern int nnodbc, nseg, nnod, nodsegm, nitmax, phaseseparation, varyviscosity;
	extern int *bcnodname, *nodname, *flow_direction, *knowntyp, *bcnod, *bctyp, *nodtyp, *ista, *iend, *segname, *index;
	extern int **nodseg, **nodnod;
	extern float pi1, consthd, constvisc, totallength, tol;
	extern float *diam, *q, *qq, *diam, *segpress, *tau, *bchd, *hd, *bcprfl;
	extern double ktau, kpress, sheartarget1, presstarget1;
	extern double *length_weight, *lseg, *nodpress, *cond, *lambda, *sheartarget, *shearfac, *nodeinflow, *precond;
	extern double **hmat, **kmat, **fullmatrix, *bvector, *xvector, *dd;
	extern FILE *ofp1;

	int notconservedcount = 0, minpressnod = 0, maxpressnod = 0;
	int iseg, inod, i;
	float total_dev, shear_rms, press_rms, currentflow = 0.;
	float viscosity, facfp = pi1 * 1333. / 128. / 0.01*60. / 1.e6;
	float shearconstant = 32. / pi1 * 1.e4 / 60.; //converts shearfac (c_j) from cP/micron^3 to (dyn/cm^2)/(nL/min)
	double bcgerror = 0.;

	kpress = 1.;	//was 0.1 in previous versions

	for (iseg = 1; iseg <= nseg; iseg++) {		// conductances g_j
		if (varyviscosity == 1) viscosity = viscor(diam[iseg], hd[iseg]);
		else viscosity = constvisc;
		cond[iseg] = facfp * pow(diam[iseg], 4) / lseg[iseg] / viscosity;	// Poiseuille's Law: C = 1/R = Q/P = (pi1/128)*D^4/(mu*L)
		shearfac[iseg] = shearconstant * viscosity / pow(diam[iseg], 3);	// c_j
		sheartarget[iseg] = flow_direction[iseg] * sheartarget1;			// constant tau0j
	}

	for (inod = 1; inod <= matrixdim; inod++) {			//initialize unknown vector
		if (inod <= nnod) xvector[inod] = nodpress[inod];
		else xvector[inod] = lambda[inod - nnod];
	}
	//set up RHS of system
	bvectorc(bvector);	//also calculates hmat and  kmat
	//////////////////////
	// LU DECOMPOSITION //
	//////////////////////
	if (solvetyp == 1) {
		for (i = 1; i <= matrixdim; i++) {		// Initialize LU vectors
			idx[i] = 1;
			dd[i] = 1.;
		}
		Amatrix(fullmatrix);
		ludcmp(fullmatrix, matrixdim, idx, dd);
		lubksb(fullmatrix, matrixdim, idx, bvector);
		printf("LU decomposition complete\n");
		for (i = 1; i <= matrixdim; i++) xvector[i] = bvector[i];
	}
	///////////////////////////////////
	// CONJUGATE GRADIENT - SPARSE   //
	// Requires symmetric matrix     //
	///////////////////////////////////
	if (solvetyp == 2) sparse_cgsymm(bvector, xvector, matrixdim, tol, nitmax);
	/////////////////////////////////
	//Process results for both methods
	for (inod = 1; inod <= matrixdim; inod++) {
		if (inod <= nnod) nodpress[inod] = xvector[inod] * precond[inod];
		else lambda[inod - nnod] = xvector[inod] * precond[inod];
	}
	for (iseg = 1; iseg <= nseg; iseg++) {
		q[iseg] = (nodpress[ista[iseg]] - nodpress[iend[iseg]])*cond[iseg];			// qj = (p1-p2)*gj
		segpress[iseg] = (nodpress[ista[iseg]] + nodpress[iend[iseg]]) / 2.;		// pressure in segment j = (p1+p2)/2
		tau[iseg] = (nodpress[ista[iseg]] - nodpress[iend[iseg]])*1333.*diam[iseg] / lseg[iseg] / 4.; // tau = r*|p1-p2|/(2L)
		qq[iseg] = fabs(q[iseg]);
	}
	// check whether flow is conserved
	notconservedcount = 0;
	for (inod = 1; inod <= nnod; inod++) if (knowntyp[inod] == 1) {
		currentflow = 0.;
		for (i = 1; i <= nodtyp[inod]; i++) {
			iseg = nodseg[i][inod];
			if (ista[iseg] == inod) currentflow += q[iseg];
			else currentflow -= q[iseg];
		}
		if (fabs(currentflow) > 0.001) { 		//flow is NOT conserved
			notconservedcount++;
			if (notconservedcount == 1) printf("*** Warning: Flow conservation error node: net flow");
			if (notconservedcount % 4 == 1) printf("\n");
			printf("%i: %e ", inod, currentflow);
		}
	}
	if (notconservedcount > 0) printf("\n");
	press_rms = 0.;
	shear_rms = 0.;	// Calculate quantities to be minimized	
	for (inod = 1; inod <= nnod; inod++) press_rms += length_weight[inod] * SQR(nodpress[inod] - presstarget1);
	for (iseg = 1; iseg <= nseg; iseg++) shear_rms += lseg[iseg] * SQR(tau[iseg] - sheartarget[iseg]);
	press_rms = press_rms / totallength;
	shear_rms = shear_rms / totallength;
	total_dev = kpress * press_rms + ktau * shear_rms;
	press_rms = sqrt(press_rms);
	shear_rms = sqrt(shear_rms);

	ofp1 = fopen("Run_summary.out", "a");
	fprintf(ofp1, "%f %f %f %f\n", ktau, press_rms, shear_rms, total_dev);
	fclose(ofp1);

	if (phaseseparation == 1 && ktau > 0.) {	// update hematocrit, but not for ktau = 0 because this gives some zero flows
		putrank();
		dishem();
	}
	else for (i = 1; i <= iseg; i++) hd[iseg] = consthd;
}