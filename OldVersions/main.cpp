/************************************************************************
Main program for FlowEstimateV1
Program to estimate flows in microvascular networks
Brendan Fry, August 2009.
Modified algorithms by Bohan Li, December 2017.
Modified by Tim Secomb, August 2018.	
*************************************************************************
	Color coding for CMGUI- Three different routines for classifying
	type of boundary node:
	1. Pressure.  Compute average capillary pressure, and for larger vessels
	classify as arteriole (if higher pressure) or venule (if lower pressure).
	2. Bifurcation. Trace back pathway to nearest bifurcation.  Converging bifurcations
	are venules, and diverging bifurcations are arterioles.
	3. Directionality.  Determine flow direction of boundary segment.  Arterioles
	characterized as downward in z-direction, and venules characterized as upward.
	Primarily relevant to brain networks.
	Classifications:
	4: Inflow Venule
	5: Infow Arteriole
	6: Inflow Capillary
	7: Outflow Venule
	8: Outflow Arteriole
	9: Outflow Capillary
*************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"
#include <Windows.h>	//needed for CopyFile
//#include <boost/filesystem.hpp>	//needed for copy_file. October 2018

//namespace fs = boost::filesystem;

void input(void);
void analyzenet(void);
void setuparrays1(int nseg, int nnod);
void flow(void);
void cmgui(float *segvar);
void picturenetwork(float *nodvar, float *segvar, const char fname[]);
void writeflow();
void analyzeresults();
void putrank(void);

int mxx, myy, mzz, nseg, nnod, nnodfl, nnodbc, nnodbck, nodsegm, solvetyp;
int numberknownpress, numberunknown, matrixdim;
int inodbc, currentnod, ktausteps, maxinsideit;
int nitmax, nitmax1, nsegfl, nitmax2, ninner, nouter;
int varyviscosity, phaseseparation, varytargetshear;
int inod, outflow, pressclass, biclass, zclass, ncaps = 0, maxsegs = 5; // Added Jan 2019 for boundary classification 
int binod, biseg, loops; // Added Jan 2019 for boundary classification 
int *flow_direction, *actual_direction, *idx;
int *nk, *nodrank;
int *nodtyp, *nodout, *bcnodname, *bcnod, *bctyp;
int *nodname, *segname, *segtyp, *ista, *iend, *knowntyp;
int **segnodname, **nodseg, **nodnod;
int *nodelambda;
float pi1 = atan(1.0)*4.0, lb, maxl, alx, aly, alz;
float tol, omega, qtol, hdtol, optlam, constvisc, consthd, totallength, diamcrit;
float vplas, optw, diamthresh, mcvcorr, mcv;
float *diam, *q, *qq, *bcprfl, *nodvar, *segvar, *xsl0, *xsl1, *xsl2;
float *tau, *hd, *segpress, *bifpar, *bchd, *viscpar, *cpar, *fahrpr;
float *histogramdisplay;
float **cnode;
double sheartarget1, ktau, kpress, ktaustart, presstarget1, eps, omega1, omega2;
double qd; // Added Jan 2019 for boundary classification 
double *sheartarget, *nodpress, *cond, *shearfac, *lambda, *dd, *length_weight, *lseg, *nodeinflow;
double *precond, *bvector, *xvector, *nodpressprev, *lambdaprev, *condsum, *hfactor2sum;
double **hmat, **kmat, **fullmatrix;
double *hfactor1, *hfactor2;

FILE *ofp1;

int main(int argc, char *argv[])
{
	int insideit, iseg, i, numdirectionchange, numreversedirection, currentdirection, ncaps;
	float duration, pcaps;
	clock_t tstart, tfinish;
	bool NoOverwrite = false;

	//Create a Current subdirectory if it does not already exist. August 2017.
	DWORD ftyp = GetFileAttributesA("Current\\");
	if (ftyp != FILE_ATTRIBUTE_DIRECTORY) system("mkdir Current");

	//copy input data files to "Current" directory
	CopyFile("ContourParams.dat", "Current\\SoluteParams.dat", NoOverwrite);
	CopyFile("FlowEstParams.dat", "Current\\FlowEstParams.dat", NoOverwrite);
	CopyFile("Network.dat", "Current\\Network.dat", NoOverwrite);
	CopyFile("RheolParams.dat", "Current\\RheolParams.dat", NoOverwrite);

	solvetyp = 1;		//solvetyp: 1 = lu decomposition, 2 = sparse conjugate gradient, 3 = relaxation

	input();

	setuparrays1(nseg, nnod);

	analyzenet();

	for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = segname[iseg];
	cmgui(segvar);		//3D image of network
	picturenetwork(nodvar, segvar, "Current/NetNodesSegs.ps");	//2D projection of network
	
	tstart = clock();
	ofp1 = fopen("Current/Run_summary.out", "w");
	fprintf(ofp1, "ktau   press_rms   shear_rms   total_dev\n");
	fclose(ofp1);

	for (i = 1; i <= ktausteps; i++) {			//keep doubling ktau
		if (i == 1) ktau = ktaustart;
		else ktau = ktau * 2;

		printf("************************ ktau = %6f *********************************\n", ktau);
		insideit = 1;
		do {			// Iterate with increasing ktau values
			if (i <= 2 && insideit == 1) {		// randomize initial flow directions
				//srand(clock());		// generate a seed
				srand(11111);	//use a fixed seed
				for (iseg = 1; iseg <= nseg; iseg++) {
					if (rand() % 2 == 1) flow_direction[iseg] = 1;
					else flow_direction[iseg] = -1;
				}
			}
			flow();			//calculate flows//
			numdirectionchange = 0;
			for (iseg = 1; iseg <= nseg; iseg++) {
				currentdirection = flow_direction[iseg];		//current flow direction
				if (q[iseg] < 0.) flow_direction[iseg] = -1;	//new flow direction
				else flow_direction[iseg] = 1;
				if (currentdirection * flow_direction[iseg] == -1) numdirectionchange++;	//flow changed direction
			}
			insideit++;
			if (numdirectionchange > 0) printf("Flow direction changes: %4i\n", numdirectionchange);
		} while (insideit <= maxinsideit && numdirectionchange > 0);
	}
	printf("**************************************************************************\n");
	tfinish = clock();
	duration = (float)(tfinish - tstart) / CLOCKS_PER_SEC;
	printf("%2.1f seconds for flow estimation\n", duration);

	analyzeresults();	//statistics and histograms of resulting flows

	//classification of boundary segments
	//calculate mean pressure of capillaries
	pcaps = 0.;
	ncaps = 0;
	for (iseg = 1; iseg <= nseg; iseg++) {
		if (diam[iseg] <= diamcrit) {
			pcaps += segpress[iseg];
			ncaps++;
		}
	}
	pcaps = pcaps / ncaps;

	for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = 0; //initialize all segments to 0

	for (inodbc = 1; inodbc <= nnodbc; inodbc++) {
		inod = bcnod[inodbc];
		iseg = nodseg[1][inod];

		// Sort in to inflow and outflow, assuming everything is a capillary
		if (ista[iseg] == inod && q[iseg] > 0) outflow = 0;
		else if (iend[iseg] == inod && q[iseg] < 0) outflow = 0; //inflow
		else outflow = 1;//outflow

		pressclass = 3 * (outflow + 2);
		biclass = 3 * (outflow + 2);
		zclass = 3 * (outflow + 2);

		// For larger vessels, sort into inflowing and outflowing arterioles and venules
		if (diam[iseg] > diamcrit) {
			// Pressure Classification
			if (segpress[iseg] < pcaps) pressclass -= 2; //classify as venule
			else pressclass -= 1; //classify as arteriole

			// Bifurcation Classification
			binod = inod;
			biseg = iseg;
			loops = 0;

			while (nodtyp[binod] < 3 && loops < 5) { // Trace back to closest bifurcation
				if (ista[biseg] == binod) binod = iend[biseg];
				else if (iend[biseg] == binod) binod = ista[biseg];
				if (biseg == nodseg[1][binod]) biseg = nodseg[2][binod];
				else if (biseg == nodseg[2][binod]) biseg = nodseg[1][binod];
				loops++;
			}

			if (2 * nodout[binod] > nodtyp[binod]) biclass -= 1; // Diverging Bifurcation, classify as arteriole
			else if (loops == 5) biclass = 0;
			else biclass -= 2;

			// Direction Classification           
			if (cnode[3][iend[iseg]] - cnode[3][ista[iseg]] < 0){ //oriented in positive z direction
				if (q[iseg] > 0) zclass -= 2; //classify as venule
				else zclass -= 1; //classify as arteriole
			}
			else {
				if (q[iseg] > 0) zclass -= 1; //classify as arteriole
				else zclass -= 2; //classify as venule
			}
		}

		// Color CMGUI boundary nodes using pressure classification
		segvar[iseg] = (float)(pressclass);

		// Assign boundary node type using pressure classification
		knowntyp[inod] = pressclass;
	}

	writeflow();		//write new network.dat file
	cmgui(segvar);		//3D image of network
	picturenetwork(nodvar, segvar, "Current/NetNodesSegs.ps");	//2D projection of network

	return 0;
}