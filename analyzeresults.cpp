/************************************************************************
analyzeresults for FlowEst2018a
TWS, August 2018.
*************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "nrutil.h"

void histogram(float *var, int n, const char filename[]);
float viscor(float d, float h);

void analyzeresults()
{
	extern int nnod, nseg, varyviscosity;
	extern int *nodname, *nodtyp, *knowntyp, *flow_direction, *actual_direction, *segname;	
	extern float totallength, constvisc;
	extern float *histogramdisplay, *q, *qq, *diam, *tau, *segpress, *hd;
	extern double *nodpress, *length_weight, *lseg, *cond;
	
	int iseg, inod, numbernegativeflows = 0, numberreversedflows = 0;
	int minpressnod, maxpressnod, maxflowseg, minflowseg, maxshearseg, minshearseg;
	float meanflow, meannodpress, meansegpress, meanshear, nodpressdeviation, segpressdeviation, sheardeviation;
	float maxflow, minflow, maxpress, minpress, maxshear, minshear;
	FILE *ofp;

	for (inod = 1; inod <= nnod; inod++) histogramdisplay[inod] = nodpress[inod];
	histogram(histogramdisplay, nnod, "histogram-pressures.out");
	for (iseg = 1; iseg <= nseg; iseg++) histogramdisplay[iseg] = log10(qq[iseg] + 1.e-6);
	histogram(histogramdisplay, nseg, "histogram-logflows.out");
	for (iseg = 1; iseg <= nseg; iseg++) histogramdisplay[iseg] = log10(fabs(tau[iseg]) + 1.e-6);
	histogram(histogramdisplay, nseg, "histogram-logstress.out");

	for (iseg = 1; iseg <= nseg; iseg++) {
		if (q[iseg] < 0) flow_direction[iseg] = -1;
		else flow_direction[iseg] = 1;
		if (flow_direction[iseg] * actual_direction[iseg] == -1) numberreversedflows++;
	}
	printf("Reversed flows (final): %4i\n", numberreversedflows);

	meannodpress = 0.;
	minpress = 1.e6;
	maxpress = -1.e6;
	for (inod = 1; inod <= nnod; inod++) {		// calculate average and max/min pressures
		meannodpress += nodpress[inod] * length_weight[inod];
		if (nodpress[inod] < minpress) {
			minpress = nodpress[inod];
			minpressnod = inod;
		}
		if (nodpress[inod] > maxpress) {
			maxpress = nodpress[inod];
			maxpressnod = inod;
		}
	}
	meannodpress = meannodpress / totallength;
	nodpressdeviation = 0.;
	for (inod = 1; inod <= nnod; inod++) nodpressdeviation += length_weight[inod] * SQR(nodpress[inod] - meannodpress);
	nodpressdeviation = sqrt(nodpressdeviation / totallength);

	meanflow = 0.;
	meansegpress = 0.;
	meanshear = 0.;
	maxflow = 0.;
	minflow = 1.e6;
	maxshear = 0.;
	minshear = 1.e6;
	for (iseg = 1; iseg <= nseg; iseg++) {					//mean calculations
		meanflow += fabs(q[iseg]) * lseg[iseg];
		meansegpress += segpress[iseg] * lseg[iseg];
		meanshear += fabs(tau[iseg]) * lseg[iseg];
		if (fabs(maxflow) < qq[iseg]) {
			maxflow = qq[iseg];
			maxflowseg = iseg;
		}
		if (fabs(minflow) > qq[iseg]) {
			minflow = qq[iseg];
			minflowseg = iseg;
		}
		if (fabs(maxshear) < fabs(tau[iseg])) {
			maxshear = fabs(tau[iseg]);
			maxshearseg = iseg;
		}
		if (fabs(minshear) > fabs(tau[iseg])) {
			minshear = fabs(tau[iseg]);
			minshearseg = iseg;
		}
	}
	meansegpress = meansegpress / totallength;
	meanflow = meanflow / totallength;
	meanshear = meanshear / totallength;
	segpressdeviation = 0.;
	sheardeviation = 0.;
	for (iseg = 1; iseg <= nseg; iseg++) {					// standard deviation calculations
		segpressdeviation += lseg[iseg] * SQR(segpress[iseg] - meansegpress);
		sheardeviation += lseg[iseg] * SQR(fabs(tau[iseg]) - meanshear);
	}
	segpressdeviation = sqrt(segpressdeviation / totallength);
	sheardeviation = sqrt(sheardeviation / totallength);

	ofp = fopen("Results.out", "w");
	fprintf(ofp, "All quantities weighted by segment length\n");
	fprintf(ofp, "Segment flow mean = %6f\n", meanflow);
	fprintf(ofp, "Nodal pressure mean +- s.d.: %6f +- %6f\n", meannodpress, nodpressdeviation);
	fprintf(ofp, "Segment pressure mean +- s.d.: %6f +- %6f\n", meansegpress, segpressdeviation);
	fprintf(ofp, "Shear stress mean +- s.d.: %6f +- %6f\n", meanshear, sheardeviation);
	fprintf(ofp, "Maximum pressure: %6f at node %4i\n", maxpress, maxpressnod);
	fprintf(ofp, "Minimum pressure: %6f at node %4i\n", minpress, minpressnod);
	fprintf(ofp, "Maximum flow: %6f at segment %4i\n", maxflow, maxflowseg);
	fprintf(ofp, "Minimum flow: %6f at segment %4i\n", minflow, minflowseg);
	fprintf(ofp, "Maximum shear stress: %6f at segment %4i\n", maxshear, maxshearseg);
	fprintf(ofp, "Minimum shear stress: %6f at segment %4i\n", minshear, minshearseg);
	fprintf(ofp, "segment name flow    length    conductance pressure shear stress dishem visc\n");
	for (iseg = 1; iseg <= nseg; iseg++) {
		fprintf(ofp, "%4i %4i %6f %6f %10.2e %6f %6f %6f ", iseg, segname[iseg], qq[iseg], lseg[iseg], cond[iseg], segpress[iseg], fabs(tau[iseg]), hd[iseg]);
		if (varyviscosity == 1) fprintf(ofp, "%6f\n", viscor(diam[iseg], hd[iseg]));
		else fprintf(ofp, "%6f\n", constvisc);
	}
	fclose(ofp);
}