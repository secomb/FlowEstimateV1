/************************************************************************
relax_method - for FlowEstimateV1 - TWS May 2019
Use successive overrelaxation approach to solve linear system
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void relax_method(double *nodpress, double *lambda)
{
	extern int nnod, nseg, nitmax2;
	extern int *nodtyp, *knowntyp, *ista, *iend;
	extern int **nodseg, **nodnod;	
	extern double kpress, presstarget1, omega1, eps;
	extern double *cond, *condsum, *length_weight, *lseg, *nodeinflow;
	extern double *hfactor1, *hfactor2, *hfactor2sum;
	
	int i, inod, inod1, iseg, iter;
	double resid, lhscoeff, dnodpressmax, dlambdamax;

	for (iter = 1; iter <= nitmax2; iter++) {
		dnodpressmax = 0.;
		dlambdamax = 0.;
		/////////////// update of interior and known boundary flow pressures /////////////
		for (inod = 1; inod <= nnod; inod++) if (knowntyp[inod] == 1 || knowntyp[inod] == 2) {
			resid = nodeinflow[inod];
			for (i = 1; i <= nodtyp[inod]; i++) {
				iseg = nodseg[i][inod];
				inod1 = nodnod[i][inod];
				resid += nodpress[inod1] * cond[iseg];
			}
			if (condsum[inod] > 0.) {
				resid = resid / condsum[inod] - nodpress[inod];
				dnodpressmax = DMAX(dnodpressmax, fabs(resid));
				nodpress[inod] += omega1 * resid;
			}
		}
		/////////////// update of lambdas /////////////
		for (inod = 1; inod <= nnod; inod++) if (knowntyp[inod] == 1 || knowntyp[inod] == 2) {
			resid = -kpress * length_weight[inod] * (nodpress[inod] - presstarget1);
			for (i = 1; i <= nodtyp[inod]; i++) {
				iseg = nodseg[i][inod];
				inod1 = nodnod[i][inod];
				if (knowntyp[inod1] == 1 || knowntyp[inod1] == 2) resid -= cond[iseg] * lambda[inod1];
				if (ista[iseg] == inod) resid += hfactor1[iseg];
				else resid -= hfactor1[iseg];
				resid += hfactor2[iseg] * nodpress[inod1];
			}
			resid += condsum[inod] * lambda[inod];
			resid -= hfactor2sum[inod] * nodpress[inod];
			if (condsum[inod] > 0.) {
				resid = -resid / condsum[inod];
				dlambdamax = DMAX(dlambdamax, fabs(resid));
				lambda[inod] += omega1 * resid;
			}
		}
		/////////////// update of unknown boundary pressures /////////////
		for (inod = 1; inod <= nnod; inod++) if (knowntyp[inod] == 3) {
			resid = -kpress * length_weight[inod] * (nodpress[inod] - presstarget1);
			if (nodtyp[inod] != 1) printf("*** Error: node %i should be nodtyp = 1\n",inod);
			iseg = nodseg[1][inod];
			inod1 = nodnod[1][inod];
			if (knowntyp[inod1] == 1 || knowntyp[inod1] == 2) resid -= cond[iseg] * lambda[inod1];
			if (ista[iseg] == inod) resid += hfactor1[iseg];
			else resid -= hfactor1[iseg];
			resid += hfactor2[iseg] * (nodpress[inod1] - nodpress[inod]);
			lhscoeff = kpress * length_weight[inod] + hfactor2[iseg];
			if (lhscoeff > 0.) {
				resid = resid / lhscoeff;
				dnodpressmax = DMAX(dnodpressmax, fabs(resid));
				nodpress[inod] += 0.1 * resid;		//underrelax here
			}
		}
		/////////////////////////////////////////////////////////////////////////////////
		if (dnodpressmax < eps && dlambdamax < eps) goto converged;
	}
	printf("*** Warning: linear iteration not converged\n");
converged:;
	printf("relax: %i %e %e\n", iter, dnodpressmax, dlambdamax);
}
