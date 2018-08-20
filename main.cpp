/************************************************************************
Main program for FlowEstimateV1
Program to estimate flows in microvascular networks
Brendan Fry, August 2009.
Modified algorithms by Bohan Li, December 2017.
Modified by Tim Secomb, August 2018.
*************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

void input(void);
void analyzenet(void);
void setuparrays1(int nseg, int nnod);
void cmgui(float *segvar);
void picturenetwork(float *nodvar, float *segvar, const char fname[]);
void flowdirections(void);
void writeflow();
void analyzeresults();

int mxx, myy, mzz, nseg, nnod, nnodfl, nnodbc, nnodbck, nodsegm, solvetyp;
int numberknownpress, numberunknown, matrixdim;
int inodbc, currentnod, ktausteps, maxinsideit;
int nitmax, nitmax1, nsegfl;
int varyviscosity, phaseseparation;
int *flow_direction, *actual_direction, *idx;
int *nk, *nodrank;
int *nodtyp, *nodout, *bcnodname, *bcnod, *bctyp;
int *nodname, *segname, *segtyp, *ista, *iend, *knowntyp;
int **segnodname, **nodseg, **nodnod;
int *nodelambda;
float pi1 = atan(1.0)*4.0, lb, maxl, alx, aly, alz;
float tol, omega, qtol, hdtol, optlam, constvisc, consthd, totallength;
float vplas, optw, diamthresh, mcvcorr, mcv;
float *diam, *q, *qq, *bcprfl, *nodvar, *segvar, *xsl0, *xsl1, *xsl2;
float *tau, *hd, *segpress, *bifpar, *bchd, *viscpar, *cpar, *fahrpr;
float *histogramdisplay;
float **cnode;
double presstarget1, sheartarget1, ktau, kpress;
double *sheartarget, *nodpress, *cond, *shearfac, *lambda, *dd, *length_weight, *lseg, *nodeinflow;
double *precond, *bvector, *xvector;
double **hmat, **kmat, **fullmatrix;
FILE *ofp1;

int main(int argc, char *argv[])
{
	int iseg;
	float duration;
	clock_t tstart, tfinish;

	solvetyp = 2;		//solvetyp: 1 = lu decomposition, 2 = sparse conjugate gradient

	input();

	setuparrays1(nseg, nnod);

	analyzenet();

	tstart = clock();
	flowdirections();	//run flow estimation algorithm
	tfinish = clock();
	duration = (float)(tfinish - tstart) / CLOCKS_PER_SEC;
	printf("%2.1f seconds for flow estimation\n", duration);

	writeflow();		//write new network.dat file

	analyzeresults();	//statistics and histograms of resulting flows

	//shear stress, positive or negative according to component of flow in the x1 direction, 
	for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = SIGN(log(FMAX(1., fabs(tau[iseg]))), (cnode[1][iend[iseg]] - cnode[1][ista[iseg]]) * q[iseg]);
	cmgui(segvar);		//3D image of network
	picturenetwork(nodvar, segvar, "NetNodesSegs.ps");	//2D projection of network

 	return 0;
}