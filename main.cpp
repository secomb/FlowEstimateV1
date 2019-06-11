/************************************************************************
Main program for FlowEstimateV1
Program to estimate flows in microvascular networks
Brendan Fry, August 2009.
Modified algorithms by Bohan Li, December 2017.
Modified by Tim Secomb, August 2018.
Segment classification by Jeff Lee, Jan. 2019
*************************************************************************
Classification of boundary segments is based on: (1) inflow or outflow;
(2) diameter compared to diamcrit; (3) pressure compared to mean capillary pressure
Classification codes (revised May 2019 to give better color codes in cmgui):
	9: Inflow arteriole
	8: Outflow arteriole
	7: Inflow capillary
	6: Outflow capillary
	5: Inflow venule
	4: Outflow venule
*************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

#if defined(__linux__)
	// Requires c++17 support, should be included in all current linux releases
	#include <experimental/filesystem> 
	namespace fs = std::experimental::filesystem::v1;
#elif defined(__APPLE__)
	// Requires removal of the -lstdc++fs flag from makefile
	#include <filesystem>
	namespace fs = std::filesystem;
#elif defined(_WIN32)    //Windows version
	#include <Windows.h>
#endif

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
int insideit, numberknownpress, numberunknown, matrixdim;
int inodbc, currentnod, ktausteps, maxinsideit;
int nitmax, nitmax1, nsegfl, nitmax2, seed;
int varyviscosity, phaseseparation, varytargetshear;  
int *flow_direction, *actual_direction, *idx, *nk, *nodrank;
int *nodtyp, *nodout, *bcnodname, *bcnod, *bctyp, *known_flow_direction;
int *nodname, *segname, *segtyp, *ista, *iend, *knowntyp, *nodelambda;
int **segnodname, **nodseg, **nodnod;
float pi1 = atan(1.0)*4.0, lb, maxl, alx, aly, alz, known_flow_weight, qf, kappa;
float tol, omega, qtol, hdtol, optlam, constvisc, consthd, totallength, diamcrit;
float vplas, optw, diamthresh, mcvcorr, mcv;
float *diam, *q, *qq, *bcprfl, *nodvar, *segvar, *xsl0, *xsl1, *xsl2;
float *tau, *hd, *segpress, *bifpar, *bchd, *viscpar, *cpar, *fahrpr, *histogramdisplay;
float **cnode;
double sheartarget1, ktau, kpress, ktaustart, presstarget1, eps, omega1;
double *sheartarget, *nodpress, *cond, *shearfac, *lambda, *dd, *length_weight, *lseg, *nodeinflow;
double *precond, *bvector, *xvector,  *condsum;
double *hfactor1, *hfactor2, *hfactor2sum;
double **hmat, **kmat, **fullmatrix;

FILE *ofp1;

int main(int argc, char *argv[])
{
	int inod, iseg, i, numdirectionchange, ncaps, outflow, pressclass;
	float duration, pcaps;
	clock_t tstart, tfinish;
	FILE *ofp;

	//Create a Current subdirectory if needed. Copy data files to it.
	#if defined(__unix__)
		if (!fs::exists("Current")) fs::create_directory("Current");
		fs::copy_file("ContourParams.dat", fs::path("Current/ContourParams.dat"), fs::copy_options::overwrite_existing);
		fs::copy_file("FlowEstParams.dat", fs::path("Current/FlowEstParams.dat"), fs::copy_options::overwrite_existing);
		fs::copy_file("Network.dat", fs::path("Current/Network.dat"), fs::copy_options::overwrite_existing);
		fs::copy_file("RheolParams.dat", fs::path("Current/RheolParams.dat"), fs::copy_options::overwrite_existing);
	#elif defined(_WIN32)
		BOOL NoOverwrite = FALSE;
		DWORD ftyp = GetFileAttributesA("Current\\");
		if (ftyp != FILE_ATTRIBUTE_DIRECTORY) system("mkdir Current");
		CopyFile("ContourParams.dat", "Current\\ContourParams.dat", NoOverwrite);
		CopyFile("FlowEstParams.dat", "Current\\FlowEstParams.dat", NoOverwrite);
		CopyFile("Network.dat", "Current\\Network.dat", NoOverwrite);
		CopyFile("RheolParams.dat", "Current\\RheolParams.dat", NoOverwrite);
	#endif
		
	solvetyp = 2;		//solvetyp: 1 = lu decomposition, 2 = sparse conjugate gradient

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

	kappa = 1.;		//factor to cancel bias in shear stress optimization

	for (i = 1; i <= ktausteps; i++) {			//keep doubling ktau
		if (i == 1) ktau = ktaustart;
		else ktau = ktau * 2;

		printf("************************ ktau = %6f *********************************\n", ktau);
		insideit = 1;
		do {
			if (i == 1 && insideit == 1) {		// randomize initial flow directions
				if(seed == 0) srand(clock());		// generate a random seed
				else srand(seed);	//use a fixed seed
				for (iseg = 1; iseg <= nseg; iseg++) {
					if (known_flow_direction[iseg] == 0) {				//not a known flow direction
						if (rand() % 2 == 1) flow_direction[iseg] = 1;	//target flow direction
						else flow_direction[iseg] = -1;
					}
					else flow_direction[iseg] = known_flow_direction[iseg];	//known flow direction
				}
			}
			///////////////////////////////////
			flow();	     //calculate flows   //
			///////////////////////////////////
			numdirectionchange = 0;
			for (iseg = 1; iseg <= nseg; iseg++){
				qf = q[iseg] * flow_direction[iseg];
				if (qf < 0.) {	//current flow is opposite to target
					if (known_flow_direction[iseg] == 0) {					//don't do this for known flow directions
						flow_direction[iseg] = -flow_direction[iseg];		//update target flow direction
						numdirectionchange++;								//flow changed direction
					}
					else printf("*** Warning: Flow in segment %i is not in specified direction\n", segname[iseg]);
				}
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

	//Classification of boundary segments
	pcaps = 0.;	//calculate mean pressure of flowing capillaries
	ncaps = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
		if (diam[iseg] <= diamcrit) {
			pcaps += segpress[iseg];
			ncaps++;
		}
	}
	pcaps = pcaps / ncaps; 
	for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = 0.;	//initialize interior segments
	for (inodbc = 1; inodbc <= nnodbc; inodbc++) {
		inod = bcnod[inodbc];
		iseg = nodseg[1][inod];
		if (ista[iseg] == inod && q[iseg] > 0) outflow = 0;
		else if (iend[iseg] == inod && q[iseg] < 0) outflow = 0;
		else outflow = 1;
		if (diam[iseg] < diamcrit) {
			if (outflow) pressclass = 6;		//outflow capillary
			else pressclass = 7;				//inflow capillary
		}
		else {
			if (segpress[iseg] > pcaps) {	
				if (outflow) pressclass = 8;	//outflow arteriole
				else pressclass = 9;			//inflow arteriole
			}
			else {
				if (outflow) pressclass = 4;	//outflow venule
				else pressclass = 5;			//inflow venule
			}
		}		
		segvar[iseg] = (float)(pressclass - 3);		// Color code boundary segments by classification
		knowntyp[inod] = pressclass;			// Assign boundary node type by classification
	}

	writeflow();		//write new network.dat file

#if defined(__linux__) 
//linux code goes here
#elif defined(_WIN32)			//Windows version
	CopyFile("NetworkNew.dat", "Current\\NetworkNew.dat", NoOverwrite);
#endif
	
	cmgui(segvar);		//write input files for cmgui visualization
	
	picturenetwork(nodvar, segvar, "Current/BoundaryNodeTypes.ps");	//2D projection of network

	return 0;
}