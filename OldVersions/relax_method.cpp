/************************************************************************
relax_method - for FlowEstimateV1 - TWS May 2019
Use relaxation approach to solve linear system. See notes.
Original version with inner and outer iteration
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void relax_method(double *nodpress, double *lambda)
{
	extern int nnod, nseg, ninner, nouter;
	extern int *nodtyp, *knowntyp, *nodelambda, *ista, *iend;
	extern int **nodseg, **nodnod;	
	extern double ktau, kpress, presstarget1, omega1, omega2, eps;
	extern double *cond, *condsum, *shearfac, *sheartarget, *length_weight, *lseg, *nodeinflow, *hfactor1, *hfactor2;
	extern double *nodpressprev, *lambdaprev, *hfactor2sum;
	
	int i, inod, inod1, iseg, iter, iter1, iter2, iter3, errnode_p, errnode_l;
	double maxerr1, maxerr2, maxerr3, resid, lhscoeff, dnodpressmax, dlambdamax;
	double dnodpress, dlambda;

	for (iter = 1; iter <= nouter; iter++) {
		for (inod = 1; inod <= nnod; inod++) {
			nodpressprev[inod] = nodpress[inod];
			lambdaprev[inod] = lambda[inod];
		}
		/////////////// STEP 1 - iterative update of pressures /////////////
		for (iter1 = 1; iter1 <= ninner; iter1++) {
			maxerr1 = 0.;
			for (inod = 1; inod <= nnod; inod++) if (knowntyp[inod] == 1 || knowntyp[inod] == 2) {
				resid = nodeinflow[inod];
				for (i = 1; i <= nodtyp[inod]; i++) {
					iseg = nodseg[i][inod];
					inod1 = nodnod[i][inod];
					resid += nodpress[inod1] * cond[iseg];
				}
				if (condsum[inod] > 0.) {
					resid = resid / condsum[inod] - nodpress[inod];
					maxerr1 = DMAX(maxerr1, fabs(resid));
					nodpress[inod] += omega1 * resid;
				}
			}
			if (maxerr1 < eps) goto converged1;
		}
	converged1:;
		/////////////// STEP 2 - iterative update of lambdas /////////////
		for (iter2 = 1; iter2 <= ninner; iter2++) {		//twice as many iterations here
			maxerr2 = 0.;
			for (inod = 1; inod <= nnod; inod++) if (knowntyp[inod] == 1 || knowntyp[inod] == 2) {
				resid = -kpress * length_weight[inod] * (nodpress[inod] - presstarget1);
				for (i = 1; i <= nodtyp[inod]; i++) {
					iseg = nodseg[i][inod];
					inod1 = nodnod[i][inod];
					if (knowntyp[inod1] == 1 || knowntyp[inod1] == 2) resid -= cond[iseg] * lambda[inod1];
					if (ista[iseg] == inod) resid += ktau * hfactor1[iseg];
					else resid -= ktau * hfactor1[iseg];
					resid += ktau * hfactor2[iseg] * nodpress[inod1];
				}
				resid += condsum[inod] * lambda[inod];
				resid -= ktau * hfactor2sum[inod] * nodpress[inod];
				if (condsum[inod] > 0.) {
					resid = -resid / condsum[inod];
					maxerr2 = DMAX(maxerr2, fabs(resid));
					lambda[inod] += omega1 * resid;
				}
			}
			if (maxerr2 < 100. * eps) goto converged2;		//lambda values are higher, so use higher error bound
		}
	converged2:;
		/////////////// STEP 3 - iterative update of unknown boundary pressures /////////////
		for (iter3 = 1; iter3 <= ninner; iter3++) {
			maxerr3 = 0.;
			for (inod = 1; inod <= nnod; inod++) if (knowntyp[inod] == 3) {
				resid = -kpress * length_weight[inod] * (nodpress[inod] - presstarget1);
				if (nodtyp[inod] != 1) printf("*** Error: node %i should be nodtyp = 1\n",inod);
				iseg = nodseg[1][inod];
				inod1 = nodnod[1][inod];
				if (knowntyp[inod1] == 1 || knowntyp[inod1] == 2) resid -= cond[iseg] * lambda[inod1];
				if (ista[iseg] == inod) resid += ktau * hfactor1[iseg];
				else resid -= ktau * hfactor1[iseg];
				resid += ktau * hfactor2[iseg] * (nodpress[inod1] - nodpress[inod]);
				lhscoeff = kpress * length_weight[inod] + ktau * hfactor2[iseg];
				if (lhscoeff > 0.) {
					resid = resid / lhscoeff;
					maxerr3 = DMAX(maxerr3, fabs(resid));
					nodpress[inod] += omega1 * resid;
				}
			}
			if (maxerr3 < eps) goto converged3;
		}
	converged3:;
		/////////////////////////////////////////////////////////////////////////////////
		//if(iter % 100 == 0) printf("(iter, max err)  1: (%i, %10e), 2: (%i, %10e), 3: (%i, %10e)\n",iter1, maxerr1, iter2, maxerr2, iter3, maxerr3);
		dnodpressmax = 0.;
		dlambdamax = 0.;
		for (inod = 1; inod <= nnod; inod++) if (knowntyp[inod] == 3) {
			dnodpress = nodpress[inod] - nodpressprev[inod];
			if (fabs(dnodpress) >= dnodpressmax) {
				dnodpressmax = fabs(dnodpress);
				errnode_p = inod;
			}
			nodpress[inod] = nodpressprev[inod] + omega2 * dnodpress;	//underrelax here
		}
		for (inod = 1; inod <= nnod; inod++) if (knowntyp[inod] == 1 || knowntyp[inod] == 2) {
			dlambda = lambda[inod] - lambdaprev[inod];
			if (fabs(dlambda) >= dlambdamax) {
				dlambdamax = fabs(dlambda);
				errnode_l = inod;
			}
			lambda[inod] = lambdaprev[inod] + omega2 * dlambda;	//underrelax here
		}
		//if (iter % 100 == 0) printf("Outer loop %i: dnodpressmax = %g at node %i, dlambdamax = %g at node %i\n", iter, dnodpressmax, errnode_p, dlambdamax, errnode_l);
		if (dnodpressmax < eps && dlambdamax < eps) goto converged;
	}
	printf("*** Warning: linear iteration not converged, dnodpressmax = %g, dlambdamax = %g\n", dnodpressmax, dlambdamax);
converged:;
	printf("relax: %i %e %e\n", iter, dnodpressmax, dlambdamax);
}
