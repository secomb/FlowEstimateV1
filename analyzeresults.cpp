/************************************************************************
analyzeresults for FlowEstV1a
TWS, August 2018, May 2019.
*************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "nrutil.h"

void histogram(float *var, float *weight, int n, const char filename[]);
float viscor(float d, float h);

void analyzeresults()
{
	extern int nnod, nseg, varyviscosity, phaseseparation, *ista, *iend, *segclasstyp, *nodclasstyp;
	extern int *nodname, *nodtyp, *nodout, *knowntyp, *flow_direction, *actual_direction, *segname, *segtyp;
	extern int **nodnod, **nodseg;
	extern float totallength, constvisc, pi1;
	extern float *histogramdisplay, *histogramweight, *q, *qq, *diam, *tau, *segpress, *hd;
	extern double *nodpress, *length_weight, *lseg, *cond;
	
	int i, iseg, inod, numbernegativeflows = 0, numberreversedflows = 0;
	int minpressnod, maxpressnod, maxflowseg, minflowseg, maxshearseg, minshearseg;
	int inod1, inod2, inod_test, iseg_test, previnod_test, inod1_orig, ncount;
	int *nlength;
	float maxflow, minflow, maxpress, minpress, maxshear, minshear, velocity, meannodpress, nodpressdeviation;
	float lengthsum;
	float *totallength1, *meanflow, *flowdeviation, *meansegpress, *segpressdeviation, *meanshear;
	float *sheardeviation, *meanHD, *HDdeviation, *meanshearrate, *shearratedeviation;
	float *meanvelocity, *velocitydeviation, *meanRBCflux, *RBCfluxdeviation, *meandiam, *diamdeviation;
	float *meanlength, *lengthdeviation;
	FILE *ofp;

	for (iseg = 1; iseg <= nseg; iseg++) {
		segpress[iseg] = (nodpress[ista[iseg]] + nodpress[iend[iseg]]) / 2;
		qq[iseg] = fabs(q[iseg]);
		if (q[iseg] < 0) flow_direction[iseg] = -1;
		else flow_direction[iseg] = 1;
		if (flow_direction[iseg] * actual_direction[iseg] == -1) numberreversedflows++;
	}
	//histograms
	for (inod = 1; inod <= nnod; inod++) {
		histogramdisplay[inod] = nodpress[inod];
		histogramweight[inod] = length_weight[inod];
	}
	histogram(histogramdisplay, histogramweight, nnod, "Current/histogram-pressures.out");

	i = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if(qq[iseg] > 0.) {
		i++;
		histogramdisplay[i] = log10(qq[iseg]);
		histogramweight[i] = lseg[iseg];
	}
	histogram(histogramdisplay, histogramweight, i, "Current/histogram-logflows.out");

	i = 0;	//velocities in mm/s
	for (iseg = 1; iseg <= nseg; iseg++) if (qq[iseg] > 0.) {
		i++;
		histogramdisplay[i] = log10(4000. / 60. * qq[iseg] / pi1 / SQR(diam[iseg]));
		histogramweight[i] = lseg[iseg];
	}
	histogram(histogramdisplay, histogramweight, i, "Current/histogram-logvelocities.out");

	i = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (qq[iseg] > 0.) {
		i++;
		histogramdisplay[i] = log10(fabs(tau[iseg]));
		histogramweight[i] = lseg[iseg];
	}
	histogram(histogramdisplay, histogramweight, i, "Current/histogram-logstress.out");

	if (phaseseparation) {
		i = 0;
		for (iseg = 1; iseg <= nseg; iseg++) if (qq[iseg] > 0.) {
			i++;
			histogramdisplay[i] = hd[iseg];
			histogramweight[i] = lseg[iseg];
		}
		histogram(histogramdisplay, histogramweight, i, "Current/histogram-hematocrit.out");
	}

	meannodpress = 0.;
	nodpressdeviation = 0.;
	minpress = 1.e6;
	maxpress = -1.e6;
	for (inod = 1; inod <= nnod; inod++) {		// calculate average and max/min pressures
		meannodpress += length_weight[inod] * nodpress[inod];
		nodpressdeviation += length_weight[inod] * SQR(nodpress[inod]);
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
	nodpressdeviation = sqrt(nodpressdeviation / totallength - SQR(meannodpress));

/**********************************************************************************
Create statistics based on segment type
Node classification
1: diverging
2: converging
0: other
Segment classification
0: O = other     (any other at start or end of the segment)
1: 2->1, 2->0, 0->1, 0->0 M = mesh (c,d)
2: 2->2 V = venule   (c=converging,c)
3: 1->2, 0->2, 1->0 C = capillary (d,c)
4: 1->1 A = arteriole (d=diverging node,d)
*********************************************************************************/
	//Setup nodclasstype 
	for (inod = 1; inod <= nnod; inod++) {
		nodclasstyp[inod] = 0;			//type 0 - not classifiable
		if (nodtyp[inod] > 2) {
			if (nodout[inod] == nodtyp[inod] - 1) nodclasstyp[inod] = 1;	//diverging
			if (nodout[inod] == 1) nodclasstyp[inod] = 2;	//converging			
		}
		else if (nodtyp[inod] == 1) {
			if (nodout[inod] == 1) nodclasstyp[inod] = 1;	//classify inflow boundary node as diverging
			if (nodout[inod] == 0) nodclasstyp[inod] = 2;	//classify outflow boundary node as converging
		}
	}
	//Classify segments as arteriole, capillary, venule
	for (iseg = 1; iseg <= nseg; iseg++) if (qq[iseg] > 0.) {
		segclasstyp[iseg] = 0;
		if (q[iseg] > 0.) {
			inod1 = ista[iseg];
			inod2 = iend[iseg];
		}
		else {
			inod1 = iend[iseg];
			inod2 = ista[iseg];
		}
		inod1_orig = inod1;	//flow is from this node
		if (nodtyp[inod1] == 2) {	//search upstream for a node with 1, 3 or more segments
			ncount = 0;
			previnod_test = inod2;
			while (nodtyp[inod1] == 2 && ncount < nseg) {
				for (i = 1; i <= 2; i++) if (nodnod[i][inod1] != previnod_test) inod_test = nodnod[i][inod1];
				previnod_test = inod1;
				inod1 = inod_test;
				ncount++;
			}
		}
		if (nodtyp[inod2] == 2) {//search downstream for a node with 1, 3 or more segments
			ncount = 0;
			previnod_test = inod1_orig;
			while (nodtyp[inod2] == 2 && ncount < nseg) {
				for (i = 1; i <= 2; i++) if (nodnod[i][inod2] != previnod_test) inod_test = nodnod[i][inod2];
				previnod_test = inod2;
				inod2 = inod_test;
				ncount++;
			}
		}
		segclasstyp[iseg] = 0;
		if (nodclasstyp[inod1] == 0) {
			if (nodclasstyp[inod2] == 0) segclasstyp[iseg] = 1;
			else if (nodclasstyp[inod2] == 1) segclasstyp[iseg] = 1;
			else if (nodclasstyp[inod2] == 2) segclasstyp[iseg] = 3;
		}
		else if (nodclasstyp[inod1] == 1) {
			if (nodclasstyp[inod2] == 0) segclasstyp[iseg] = 3;
			else if (nodclasstyp[inod2] == 1) segclasstyp[iseg] = 4;
			else if (nodclasstyp[inod2] == 2) segclasstyp[iseg] = 3;
		}
		else if (nodclasstyp[inod1] == 2) {
			if (nodclasstyp[inod2] == 0) segclasstyp[iseg] = 1;
			else if (nodclasstyp[inod2] == 1) segclasstyp[iseg] = 1;
			else if (nodclasstyp[inod2] == 2) segclasstyp[iseg] = 2;
		}
	}

	totallength1 = vector(0, 5);
	nlength = ivector(0, 5);
	meanlength = vector(0, 5);
	lengthdeviation = vector(0, 5);
	meandiam = vector(0, 5);
	diamdeviation = vector(0, 5);
	meanflow = vector(0, 5);
	flowdeviation = vector(0, 5);
	meansegpress = vector(0, 5);
	segpressdeviation = vector(0, 5);
	meanshear = vector(0, 5);
	sheardeviation = vector(0, 5);
	meanHD = vector(0, 5);
	HDdeviation = vector(0, 5);
	meanshearrate = vector(0, 5);
	shearratedeviation = vector(0, 5);
	meanvelocity = vector(0, 5);
	velocitydeviation = vector(0, 5);
	meanRBCflux = vector(0, 5);
	RBCfluxdeviation = vector(0, 5);

	for (i = 0; i <= 5; i++) {
		totallength1[i] = 0.;
		nlength[i] = 0;
		meandiam[i] = 0.;
		diamdeviation[i] = 0.;
		meanlength[i] = 0.;
		lengthdeviation[i] = 0.;
		meanflow[i] = 0.;
		flowdeviation[i] = 0.;
		meansegpress[i] = 0.;
		segpressdeviation[i] = 0.;
		meanshear[i] = 0.;
		sheardeviation[i] = 0.;
		meanHD[i] = 0.;
		HDdeviation[i] = 0.;
		meanshearrate[i] = 0.;
		shearratedeviation[i] = 0.;
		meanvelocity[i] = 0.;
		velocitydeviation[i] = 0.;
		meanRBCflux[i] = 0.;
		RBCfluxdeviation[i] = 0.;
	}
	maxflow = 0.;
	minflow = 1.e6;
	maxshear = 0.;
	minshear = 1.e6;

	for (iseg = 1; iseg <= nseg; iseg++) if (qq[iseg] > 0.) {		//mean calculations
		i = segclasstyp[iseg];
		meandiam[i] += lseg[iseg] * diam[iseg];
		meandiam[5] += lseg[iseg] * diam[iseg];
		diamdeviation[i] += lseg[iseg] * SQR(diam[iseg]);
		diamdeviation[5] += lseg[iseg] * SQR(diam[iseg]);
		meanflow[i] += lseg[iseg] * qq[iseg];
		meanflow[5] += lseg[iseg] * qq[iseg];
		flowdeviation[i] += lseg[iseg] * SQR(qq[iseg]);
		flowdeviation[5] += lseg[iseg] * SQR(qq[iseg]);
		meansegpress[i] += lseg[iseg] * segpress[iseg];
		meansegpress[5] += lseg[iseg] * segpress[iseg];
		segpressdeviation[i] += lseg[iseg] * SQR(segpress[iseg]);
		segpressdeviation[5] += lseg[iseg] * SQR(segpress[iseg]);
		meanshear[i] += lseg[iseg] * fabs(tau[iseg]);
		meanshear[5] += lseg[iseg] * fabs(tau[iseg]);
		sheardeviation[i] += lseg[iseg] * SQR(fabs(tau[iseg]));
		sheardeviation[5] += lseg[iseg] * SQR(fabs(tau[iseg]));
		meanHD[i] += lseg[iseg] * hd[iseg];
		meanHD[5] += lseg[iseg] * hd[iseg];
		HDdeviation[i] += lseg[iseg] * SQR(hd[iseg]);
		HDdeviation[5] += lseg[iseg] * SQR(hd[iseg]);
		velocity = 1.e6 / 60. * qq[iseg] / SQR(diam[iseg]);	//in um/s
		meanvelocity[i] += lseg[iseg] * 1.e-3 * velocity;	//convert to mm/s
		meanvelocity[5] += lseg[iseg] * 1.e-3 * velocity;	//convert to mm/s
		velocitydeviation[i] += lseg[iseg] * SQR(1.e-3 * velocity);
		velocitydeviation[5] += lseg[iseg] * SQR(1.e-3 * velocity);
		meanshearrate[i] += lseg[iseg] * velocity / diam[iseg];
		meanshearrate[5] += lseg[iseg] * velocity / diam[iseg];
		shearratedeviation[i] += lseg[iseg] * SQR(velocity / diam[iseg]);
		shearratedeviation[5] += lseg[iseg] * SQR(velocity / diam[iseg]);
		meanRBCflux[i] += lseg[iseg] * qq[iseg] * hd[iseg];
		meanRBCflux[5] += lseg[iseg] * qq[iseg] * hd[iseg];
		RBCfluxdeviation[i] += lseg[iseg] * SQR(qq[iseg] * hd[iseg]);
		RBCfluxdeviation[5] += lseg[iseg] * SQR(qq[iseg] * hd[iseg]);
		totallength1[i] += lseg[iseg];
		totallength1[5] += lseg[iseg];
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
	for (i = 0; i <= 5; i++) if(totallength1[i] > 0.) {
		meandiam[i] = meandiam[i] / totallength1[i];
		diamdeviation[i] = sqrt(diamdeviation[i] / totallength1[i] - SQR(meandiam[i]));
		meanflow[i] = meanflow[i] / totallength1[i];
		flowdeviation[i] = sqrt(flowdeviation[i] / totallength1[i] - SQR(meanflow[i]));
		meansegpress[i] = meansegpress[i] / totallength1[i];
		segpressdeviation[i] = sqrt(segpressdeviation[i] / totallength1[i] - SQR(meansegpress[i]));
		meanshear[i] = meanshear[i] / totallength1[i];
		sheardeviation[i] = sqrt(sheardeviation[i] / totallength1[i] - SQR(meanshear[i]));
		meanHD[i] = meanHD[i] / totallength1[i];
		HDdeviation[i] = sqrt(HDdeviation[i] / totallength1[i] - SQR(meanHD[i]));
		meanvelocity[i] = meanvelocity[i] / totallength1[i];
		velocitydeviation[i] = sqrt(velocitydeviation[i] / totallength1[i] - SQR(meanvelocity[i]));
		meanshearrate[i] = meanshearrate[i] / totallength1[i];
		shearratedeviation[i] = sqrt(shearratedeviation[i] / totallength1[i] - SQR(meanshearrate[i]));
		meanRBCflux[i] = meanRBCflux[i] / totallength1[i];
		RBCfluxdeviation[i] = sqrt(RBCfluxdeviation[i] / totallength1[i] - SQR(meanRBCflux[i]));
	}
	//calculation of segment lengths, combining segments at type 2 nodes
	for (iseg = 1; iseg <= nseg; iseg++) if (qq[iseg] > 0.) segclasstyp[iseg] += 10;	//mark as not processed yet
	for (iseg = 1; iseg <= nseg; iseg++) if (qq[iseg] > 0. && segclasstyp[iseg] >= 10) {
		inod1 = ista[iseg];
		inod2 = iend[iseg];
		lengthsum = lseg[iseg];
		segclasstyp[iseg] -= 10;			//mark as processed
		inod1_orig = inod1;
		if (nodtyp[inod1] == 2) {	//search for a node with 1, 3 or more segments
			ncount = 0;
			previnod_test = inod2;
			while (nodtyp[inod1] == 2 && ncount < nseg) {
				for (i = 1; i <= 2; i++) if (nodnod[i][inod1] != previnod_test) {
					inod_test = nodnod[i][inod1];
					iseg_test = nodseg[i][inod1];
					lengthsum += lseg[iseg_test];
					segclasstyp[iseg_test] -= 10;
				}
				previnod_test = inod1;
				inod1 = inod_test;
				ncount++;
			}
		}
		if (nodtyp[inod2] == 2) {//search for a node with 1, 3 or more segments
			ncount = 0;
			previnod_test = inod1_orig;
			while (nodtyp[inod2] == 2 && ncount < nseg) {
				for (i = 1; i <= 2; i++) if (nodnod[i][inod2] != previnod_test) {
					inod_test = nodnod[i][inod2];
					iseg_test = nodseg[i][inod2];
					lengthsum += lseg[iseg_test];
					segclasstyp[iseg_test] -= 10;
				}
				previnod_test = inod2;
				inod2 = inod_test;
				ncount++;
			}
		}
		i = segclasstyp[iseg];
		nlength[i]++;
		nlength[5]++;
		meanlength[i] += lengthsum;
		meanlength[5] += lengthsum;
		lengthdeviation[i] += SQR(lengthsum);
		lengthdeviation[5] += SQR(lengthsum);
	}
	if(fabs(meanlength[5] - totallength1[5]) > 1.) printf("*** Error: length mismatch %f %f\n", meanlength[5], totallength1[5]);
	for (i = 0; i <= 5; i++) if (nlength[i] > 0.) {
		meanlength[i] = meanlength[i] / nlength[i];
		lengthdeviation[i] = sqrt(lengthdeviation[i] / nlength[i] - SQR(meanlength[i]));
	}

	ofp = fopen("Current/Results.out", "w");
	fprintf(ofp, "All quantities (except length) weighted by segment length\n");
	fprintf(ofp, "Nodal pressure (mmHg) mean +- s.d.: %6f +- %6f CV = %6f\n", meannodpress, nodpressdeviation, nodpressdeviation / meannodpress);
	for (i = 0; i <= 5; i++) if (totallength1[i] > 0.) {
		switch (i) {
		case 0:
			fprintf(ofp, "\nOther_vessels\n");
			break;
		case 1:
			fprintf(ofp, "\nMesh_vessels\n");
			break;
		case 2:
			fprintf(ofp, "\nVenules\n");
			break;
		case 3:
			fprintf(ofp, "\nCapillaries\n");
			break;
		case 4:
			fprintf(ofp, "\nArterioles\n");
			break;
		case 5:
			fprintf(ofp, "\nAll_vessels\n");
			break;
		}
		fprintf(ofp, "Total_number (combining_type_2_nodes): %i\n", nlength[i]);
		fprintf(ofp, "Total_length (um): %2f\n", totallength1[i]);
		fprintf(ofp, "Segment_diam (um) mean+-s.d.: %6f +- %6f CV: %6f\n", meandiam[i], diamdeviation[i], diamdeviation[i] / meandiam[i]);
		fprintf(ofp, "Segment_length (um) mean+-s.d.: %6f +- %6f CV: %6f\n", meanlength[i], lengthdeviation[i], lengthdeviation[i] / meanlength[i]);
		fprintf(ofp, "Segment_flow (nl/min) mean+-s.d.: %6f +- %6f CV: %6f\n", meanflow[i], flowdeviation[i], flowdeviation[i] / meanflow[i]);
		fprintf(ofp, "Segment_pressure (mmHg) mean+-s.d.: %6f +- %6f CV: %6f\n", meansegpress[i], segpressdeviation[i], segpressdeviation[i] / meansegpress[i]);
		fprintf(ofp, "Shear_stress (dyn/cm2) mean+-s.d.: %6f +- %6f CV: %6f\n", meanshear[i], sheardeviation[i], sheardeviation[i] / meanshear[i]);
		fprintf(ofp, "Hematocrit (HD) mean+-s.d.: %6f +- %6f CV: %6f\n", meanHD[i], HDdeviation[i], HDdeviation[i] / meanHD[i]);
		fprintf(ofp, "Velocity (mm/s) mean+-s.d.: %6f +- %6f CV: %6f\n", meanvelocity[i], velocitydeviation[i], velocitydeviation[i] / meanvelocity[i]);
		fprintf(ofp, "Shear_rate (V/D,s-1) mean+-s.d.: %6f +- %6f CV: %6f\n", meanshearrate[i], shearratedeviation[i], shearratedeviation[i] / meanshearrate[i]);
		fprintf(ofp, "RBC_flux (nl/min) mean+-s.d.: %6f +- %6f CV: %6f\n", meanRBCflux[i], RBCfluxdeviation[i], RBCfluxdeviation[i] / meanRBCflux[i]);
	}

	fprintf(ofp, "Maximum_pressure (mmHg): %6f at_node %4i\n", maxpress, nodname[maxpressnod]);
	fprintf(ofp, "Minimum_pressure (mmHg): %6f at_node %4i\n", minpress, nodname[minpressnod]);
	fprintf(ofp, "Maximum_flow (nl/min): %6f at_segment %4i\n", maxflow, segname[maxflowseg]);
	fprintf(ofp, "Minimum_flow (nl/min): %6f at_segment %4i\n", minflow, segname[minflowseg]);
	fprintf(ofp, "Maximum_shear_stress (dyn/cm2): %6f at_segment %4i\n", maxshear, segname[maxshearseg]);
	fprintf(ofp, "Minimum_shear_stress (dyn/cm2): %6f at_segment %4i\n", minshear, segname[minshearseg]);

	fprintf(ofp, "\nsegment name flow    length    conductance pressure shear stress dishem visc\n");
	for (iseg = 1; iseg <= nseg; iseg++) {
		fprintf(ofp, "%4i %4i %6f %6f %10.2e %6f %6f %6f ", iseg, segname[iseg], q[iseg], lseg[iseg], cond[iseg], segpress[iseg], tau[iseg], hd[iseg]);
		if (varyviscosity == 1) fprintf(ofp, "%6f\n", viscor(diam[iseg], hd[iseg]));
		else fprintf(ofp, "%6f\n", constvisc);
	}
	fclose(ofp);
}
