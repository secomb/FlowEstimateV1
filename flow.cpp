/************************************************************************
flow - for FlowEstV1
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
			3 successive over-relaxation
************************************************************************
Factor kappa multiplying target shear stress introduced to compensate
for bias in least squares estimation. 
kappa = mean(tau^2)/(mean(tau))^2
Otherwise, estimated tau values are biased to be small, because this 
reduces the variance in tau. TWS, May 2019.
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

void bvectorc(double *output);
void Amatrix(double **matrix);
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double *b);
double sparse_cgsymm(double *b, double *x, int n, double eps, int itmax);
void relax_method(double *nodpress, double *lambda);
void putrank(void);
void dishem(void);
float viscor(float d, float h);
void Amultiply(double *input, double *output);

void flow()
{
	extern int numberknownpress, numberunknown, matrixdim, solvetyp, varytargetshear, nitmax2;
	extern int insideit, nnodbc, nseg, nnod, nodsegm, phaseseparation, varyviscosity, ktausteps;
	extern int *bcnodname, *nodname, *flow_direction, *segtyp, *knowntyp, *bcnod, *bctyp, *nodtyp, *ista, *iend, *segname, *idx;
	extern int **nodseg, **nodnod, *known_flow_direction;
	extern float pi1, consthd, constvisc, totallength, known_flow_weight, kappa;
	extern float *diam, *q, *diam, *tau, *bchd, *hd, *bcprfl;
    extern float *tmpq0;
	extern double ktau, kpress, sheartarget1, presstarget1, eps, omega1, ktaustart;
	extern double *length_weight, *lseg, *nodpress, *cond, *lambda, *sheartarget, *shearfac, *nodeinflow, *precond;
	extern double *hfactor1, *hfactor2, *condsum, *hfactor2sum;
	extern double **hmat, **kmat, **fullmatrix, *bvector, *xvector, *dd;
	extern FILE *ofp1;

	int notconservedcount = 0, minpressnod = 0, maxpressnod = 0;
	int iseg, inod, i;
	float total_dev, shear_rms, press_rms, currentflow, sumflow;
	float viscosity, facfp = pi1 * 1333. / 128. / 0.01*60. / 1.e6;
	float shearconstant = 32. / pi1 * 1.e4 / 60.; //converts shearfac (c_j) from cP/micron^3 to (dyn/cm^2)/(nL/min)
	float kappasum1, kappasum2;
	double vesspress, vess_shear_target;
	double bcgerror = 0.;
	float duration;
	clock_t tstart, tfinish;

    float known_scale_fac = 10.;
    float known_velocity_weight = FMIN(known_scale_fac * ktau, pow(2.,(float)ktausteps - 1.) * ktaustart);
    
    // Hardcoded solution for forcing higher scale fac at low ktau values
    //if (ktau < 0.128) known_velocity_weight = known_scale_fac * ktau;
    //else known_velocity_weight = ktau;

	kpress = 1.;	//was 0.1 in previous versions

	for (iseg = 1; iseg <= nseg; iseg++) {		// conductances g_j
		if (varyviscosity == 1) viscosity = viscor(diam[iseg], hd[iseg]);
		else viscosity = constvisc;
		cond[iseg] = facfp * pow(diam[iseg], 4) / lseg[iseg] / viscosity;	// Poiseuille's Law: C = 1/R = Q/P = (pi1/128)*D^4/(mu*L)
		shearfac[iseg] = shearconstant * viscosity / pow(diam[iseg], 3);	// c_j
		if (varytargetshear) {
			vesspress = (nodpress[ista[iseg]] + nodpress[iend[iseg]]) / 2;		// compute vessel pressure as average of nodes
			vesspress = FMAX(vesspress, 10.);									// set minimum pressure to 10 for shear target estimation
            if (segtyp[iseg] == 6) sheartarget[iseg] = 80 * viscosity * tmpq0[iseg] / diam[iseg];
            else {
                vess_shear_target = 100 - 86 * exp(-5000 * pow(log10(log10(vesspress)), 5.4));	//pressure-shear relationship
                sheartarget[iseg] = flow_direction[iseg] * vess_shear_target;
            }

		}
		else sheartarget[iseg] = flow_direction[iseg] * sheartarget1;  // constant tau0j

        if (segtyp[iseg] == 6) {
            //treat these the same as known flow directions???  I don't know but we'll certainly try
            //apparently this doesn't work... but treating them as normal segments does.
            //FMIN(pow(2.,ktausteps) * ktaustart, known_scale_fac * ktau) * 
			hfactor1[iseg] = known_velocity_weight * lseg[iseg] * shearfac[iseg] * kappa * cond[iseg] * sheartarget[iseg];	//extra weight for known flow direction segments
			hfactor2[iseg] = known_velocity_weight * lseg[iseg] * SQR(shearfac[iseg] * cond[iseg]);

        }
        else if (known_flow_direction[iseg] == 0) {
			hfactor1[iseg] = ktau * lseg[iseg] * shearfac[iseg] * kappa * cond[iseg] * sheartarget[iseg];	//needed for hmatrix terms
			hfactor2[iseg] = ktau * lseg[iseg] * SQR(shearfac[iseg] * cond[iseg]);
		}
		else{
			hfactor1[iseg] = known_flow_weight * shearfac[iseg] * kappa * cond[iseg] * sheartarget[iseg];	//extra weight for known flow direction segments
			hfactor2[iseg] = known_flow_weight * SQR(shearfac[iseg] * cond[iseg]);
		}
	}
	tstart = clock();
	//////////////////////
	// LU DECOMPOSITION //
	//////////////////////
	if (solvetyp == 1) {
		bvectorc(bvector);	//set up RHS of system, also calculates hmat and  kmat
		for (i = 1; i <= matrixdim; i++) {		// Initialize LU vectors
			idx[i] = 1;
			dd[i] = 1.;
		}
		Amatrix(fullmatrix);
		ludcmp(fullmatrix, matrixdim, idx, dd);
		lubksb(fullmatrix, matrixdim, idx, bvector);
		for (inod = 1; inod <= matrixdim; inod++) {
			if (inod <= nnod) nodpress[inod] = bvector[inod] * precond[inod];
			else lambda[inod - nnod] = bvector[inod] * precond[inod];
		}
	}
	///////////////////////////////////
	// CONJUGATE GRADIENT (CG)       //
	// Requires symmetric matrix     //
	///////////////////////////////////
	if (solvetyp == 2) {
		bvectorc(bvector);	//set up RHS of system, also calculates hmat and  kmat
		for (inod = 1; inod <= matrixdim; inod++) {			//initialize unknown vector
			if (inod <= nnod) xvector[inod] = nodpress[inod];
			else xvector[inod] = lambda[inod - nnod];
		}
		sparse_cgsymm(bvector, xvector, matrixdim, eps, nitmax2);
		for (inod = 1; inod <= matrixdim; inod++) {
			if (inod <= nnod) nodpress[inod] = xvector[inod] * precond[inod];
			else lambda[inod - nnod] = xvector[inod] * precond[inod];
		}
	}
	////////////////////////////////////////
	// SUCCESSIVE OVER-RELAXATION (SOR)   //
	////////////////////////////////////////
	if (solvetyp == 3) {
		for (inod = 1; inod <= nnod; inod++) {
			condsum[inod] = 0.;
			hfactor2sum[inod] = 0.;
			for (i = 1; i <= nodtyp[inod]; i++) {
				iseg = nodseg[i][inod];
				condsum[inod] += cond[iseg];
				hfactor2sum[inod] += hfactor2[iseg];
			}
		}
		relax_method(nodpress, lambda);
	}
	/////////////////////////////////////////////////
	tfinish = clock();
	duration = (float)(tfinish - tstart) / CLOCKS_PER_SEC;
	printf("Iteration %i: %2.1f seconds for solution by method %i\n", insideit, duration, solvetyp);
	///////////////// Process results ///////////////
	kappasum1 = 0.;	//update kappa factor according to results from previous iteration
	kappasum2 = 0.;
	for (iseg = 1; iseg <= nseg; iseg++) {
            q[iseg] = (nodpress[ista[iseg]] - nodpress[iend[iseg]])*cond[iseg];			// qj = (p1-p2)*gj
            tau[iseg] = (nodpress[ista[iseg]] - nodpress[iend[iseg]])*1333.*diam[iseg] / lseg[iseg] / 4.; // tau = r*|p1-p2|/(2L)
            kappasum2 += lseg[iseg] * SQR(tau[iseg]);
            kappasum1 += lseg[iseg] * fabs(tau[iseg]);
	}
	kappasum1 = kappasum1 / totallength;
	kappasum2 = kappasum2 / totallength;
	kappa = kappasum2 / SQR(kappasum1);

	notconservedcount = 0;		// check whether flow is conserved
	for (inod = 1; inod <= nnod; inod++) if (knowntyp[inod] == 1) {
		currentflow = 0.;
		sumflow = 0.;
		for (i = 1; i <= nodtyp[inod]; i++) {
			iseg = nodseg[i][inod];
			if (ista[iseg] == inod) currentflow += q[iseg];
			else currentflow -= q[iseg];
			sumflow += fabs(q[iseg]);
		}
		if (fabs(currentflow)/(1. + sumflow) > 10.*eps) { 		//flow is not conserved
			notconservedcount++;
			if (notconservedcount == 1) printf("*** Warning: Flow conservation error node: net flow");
			if (notconservedcount % 8 == 1) printf("\n");
			printf("%i: %8.1e ", inod, currentflow);
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

	ofp1 = fopen("Current/Run_summary.out", "a");
	fprintf(ofp1, "%f %f %f %f %f\n", ktau, press_rms, shear_rms, total_dev, kappa);
	fclose(ofp1);

	if (phaseseparation == 1 && ktau > 0.) {	// update hematocrit, but not for ktau = 0 because this gives some zero flows
		putrank();
		dishem();
	}
}
