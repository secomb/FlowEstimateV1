/************************************************************************
flowdirections - for FlowEstimateV1
BCF June 2010
Revised TWS, August 2018
Flows in nl/min, pressures in mmHg, viscosity in cP
Lengths, diameters in microns, times in s
Use double precision pressures for networks with many segments
Program to iteratively estimate flow directions. Outline of procedure:
(1) Start with ktau = 0
(2) Determine flows
(3) If q[iseg] is negative, set the flow direction of iseg to -1
	This means that the shear stress is negative, so the target shear stress is set to be negative
(4) Successively double ktau, and repeat (2) and (3)
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

void flow(void);

void flowdirections()
{
	extern int nnod, nseg, nnodbc, ktausteps, maxinsideit;
	extern int *actual_direction, *flow_direction;
	extern float *q;
	extern double ktau, *nodeinflow, *lambda, *nodpress;
	extern FILE *ofp1;

	int insideit, iseg, i, numdirectionchange, numreversedirection, currentdirection;

	ofp1 = fopen("Run_summary.out", "w");
	fprintf(ofp1, "ktau   press_rms   shear_rms   total_dev\n");
	fclose(ofp1);

	for (i = 1; i <= ktausteps; i++) {			//keep doubling ktau
		if (i == 1) ktau = 0.;
		if (i == 2) ktau = 0.001;
		else ktau = ktau * 2;
		printf("************************ ktau = %6f *********************************\n", ktau);
		insideit = 1;
		do {			// Iterate with increasing ktau values
			if (i <= 2 && insideit == 1) {		// randomize initial flow directions
				//srand(clock());		// generate a seed
				srand(11111);// srand(27182);// srand(31415);	//use a fixed seed
				for (iseg = 1; iseg <= nseg; iseg++) {
					if (rand() % 2 == 1) flow_direction[iseg] = 1;
					else flow_direction[iseg] = -1;
				}
			}
			///////////////////////////////////
			flow();			//calculate flows//
			///////////////////////////////////
			numdirectionchange = 0;
			for (iseg = 1; iseg <= nseg; iseg++) {
				currentdirection = flow_direction[iseg];		//current flow direction
				if (q[iseg] < 0.) flow_direction[iseg] = -1;	//new flow direction
				else flow_direction[iseg] = 1;
				if (currentdirection * flow_direction[iseg] == -1) {//flow changed direction
					numdirectionchange++;
					printf("%i ", iseg);
				}
			}
			insideit++;
			if (numdirectionchange > 0) printf("\nFlow direction changes: %4i\n", numdirectionchange);
		} while (insideit <= maxinsideit && numdirectionchange > 0);
	}
	printf("**************************************************************************\n");
}