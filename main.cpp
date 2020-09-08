/************************************************************************
Main program for FlowEstimateV1
Program to estimate flows in microvascular networks
Brendan Fry, August 2009.
Modified algorithms by Bohan Li, December 2017.
Modified by Tim Secomb, August 2018.
Segment classification by Grace Lee, Jan. 2019
*************************************************************************
Classification of boundary segments is based on: (1) inflow or outflow;
(2) diameter compared to diamcrit; (3) pressure compared to mean capillary pressure
Classification codes (revised May 2019 to give better color codes in cmgui):
	9: Inflow arteriole     Dark Red
	8: Outflow arteriole    Light Red
	7: Inflow capillary     Yellow
	6: Outflow capillary    Green
	5: Inflow venule        Light Blue
    4: Outflow venule       Dark Blue
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
	#include <filesystem> 
	namespace fs = std::filesystem;
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
void flowtest();
void cmgui(float *segvar);
void picturenetwork(float *nodvar, float *segvar, const char fname[]);
void writeflow();
void analyzeresults();
void putrank(void);

int mxx, myy, mzz, nseg, nnod, nnodfl, nnodbc, nnodbck, nodsegm, solvetyp;
int insideit, numberknownpress, numberunknown, matrixdim;
int inodbc, currentnod, ktausteps, maxinsideit;
int nitmax, nitmax1, nsegfl, nitmax2, nitmax3, seed;
int varyviscosity, phaseseparation, varytargetshear;  
int *flow_direction, *actual_direction, *idx, *nk, *nodrank, *segclasstyp, *nodclasstyp;
int *nodtyp, *nodout, *bcnodname, *bcnod, *bctyp, *known_flow_direction;
int *nodname, *segname, *segtyp, *ista, *iend, *knowntyp, *nodelambda;
int **segnodname, **nodseg, **nodnod;
float pi1 = atan(1.0)*4.0, lb, maxl, alx, aly, alz, known_flow_weight, qf, kappa;
float tol, omega, qtol, hdtol, optlam, constvisc, consthd, totallength, diamcrit, diamcorr;
float vplas, optw, diamthresh, mcvcorr, mcv;
float *diam, *q, *qq, *bcprfl, *nodvar, *segvar, *xsl0, *xsl1, *xsl2;
float *tau, *hd, *segpress, *bifpar, *bchd, *viscpar, *cpar, *fahrpr, *histogramdisplay, *histogramweight;
float **cnode;
double sheartarget1, ktau, kpress, ktaustart, presstarget1, eps, omega1, omega2;
double *sheartarget, *nodpress, *cond, *shearfac, *lambda, *dd, *length_weight, *lseg, *nodeinflow;
double *precond, *bvector, *xvector,  *condsum;
double *hfactor1, *hfactor2, *hfactor2sum;
double **hmat, **kmat, **fullmatrix;

float *tmpordart, *tmpordvein, *tmpq0;
FILE *ofp1;

int main(int argc, char *argv[])
{
	int inod, iseg, i, numdirectionchange, outflow, pressclass, ktausteps1;
	float duration, pcaps, totallength1;
	clock_t tstart, tfinish;

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
		
	input();

	setuparrays1(nseg, nnod);

    for (iseg = 1; iseg <=nseg; iseg++) diam[iseg] += diamcorr;

	if (seed == 0) srand(clock());		// generate a random seed, needed by flowtest
	else srand(seed);	//use a fixed seed

	flowtest();	//removes zero-flow segments: renumbers nodes and segments

	analyzenet();

	for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = segname[iseg];
	cmgui(segvar);		//3D image of network
	picturenetwork(nodvar, segvar, "Current/NetNodesSegs.ps");	//2D projection of network
	
	tstart = clock();
	ofp1 = fopen("Current/Run_summary.out", "w");
	fprintf(ofp1, "ktau   press_rms   shear_rms   total_dev   kappa\n");
	fclose(ofp1);

	kappa = 1.;		//initialize factor to cancel bias in shear stress optimization

	if (ktausteps > 1) ktausteps1 = ktausteps + 3;
	else ktausteps1 = 1;		//for quick tests, don't do extra iterations

	for (i = 1; i <= ktausteps1; i++) {			//keep doubling ktau, then three extra iterations
		if (i == 1) ktau = ktaustart;
		else if (i <= ktausteps) ktau = ktau * 2;

		printf("************************ ktau = %6f *********************************\n", ktau);
		insideit = 1;
		do {
			if (i == 1 && insideit == 1) {		// randomize initial flow directions
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

    int classdata = 1; //switch for known segment classification data
	pcaps = 0.;	//calculate length-weighted mean pressure of flowing capillaries
	totallength1 = 0.;
	for (iseg = 1; iseg <= nseg; iseg++) if (qq[iseg] > 0.) {
		if (diam[iseg] <= diamcrit) {
			pcaps += segpress[iseg] * lseg[iseg];
			totallength1 += lseg[iseg];
		}
	}
	pcaps = pcaps / totallength1;

    for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = 0.;

	for (i = 1; i <= nnodbck; i++) {
        for (inodbc=1; inodbc<=nnodbc; inodbc++) {
            if (bcnodname[i] == nodname[bcnod[inodbc]]) {
                inod = bcnod[inodbc];
                iseg = nodseg[1][inod];
                printf("%i %i %f\n", inod, iseg, bcprfl[i]);
                segvar[iseg] = bcprfl[i];
            }
        }
    }
    
	for (inodbc = 1; inodbc <= nnodbc; inodbc++) {
		inod = bcnod[inodbc];
		iseg = nodseg[1][inod];

		if (ista[iseg] == inod && q[iseg] > 0) outflow = 0;
		else if (iend[iseg] == inod && q[iseg] < 0) outflow = 0;
		else outflow = 1;

        if (classdata == 1) {
            if (tmpordart[iseg] == 0) {	
		    	if (outflow) pressclass = 8;	//outflow arteriole
		    	else pressclass = 9;			//inflow arteriole
		    }
		    else if (tmpordvein[iseg] == 0) {
		    	if (outflow) pressclass = 4;	//outflow venule
		    	else pressclass = 5;			//inflow venule
			}
        }
		else if (diam[iseg] < diamcrit) {
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

#if defined(__unix__)			//Apple, linux
	fs::copy_file("NetworkNew.dat", fs::path("Current/NetworkNew.dat"), fs::copy_options::overwrite_existing);
#elif defined(_WIN32)			//Windows version
	CopyFile("NetworkNew.dat", "Current/NetworkNew.dat", NoOverwrite);
#endif
	
	cmgui(segvar);		//write input files for cmgui visualization
	
	picturenetwork(nodvar, segvar, "Current/BoundaryNodeTypes.ps");	//2D projection of network
    
    float flow_direction_temp=1., velocity_temp;

    ofp1 = fopen("Current/data.out", "w");
    fprintf(ofp1, "segName diam flow press shearstress velocity\n");
    for (iseg = 1; iseg <= nseg; iseg++){
        if (cnode[3][ista[iseg]] < cnode[3][iend[iseg]]) flow_direction_temp = -1.;
        else flow_direction_temp = 1.;
    
        velocity_temp = 4000. / 60. * q[iseg] / pi1 / SQR(diam[iseg]);
        fprintf(ofp1, "%i %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", segname[iseg], diam[iseg], q[iseg], segpress[iseg], tau[iseg], flow_direction_temp * velocity_temp, velocity_temp);
    }
    fclose(ofp1); 

	return 0;
}
