/************************************************************************
setuparrays1 - for flowest09.  TWS August 09
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setuparrays1(int nseg, int nnod)
{
	extern int nodsegm, nnodbc, solvetyp;
	extern int *nodout, *nodtyp, *ista, *iend;
	extern int *flow_direction, *known_flow_direction;
	extern int *knowntyp, *nodelambda, *bcnod;
	extern int *nk, *nodrank, *segclasstyp, *nodclasstyp;
	extern int **nodnod, **nodseg;
	extern float *tau, *segpress, *qq, *nodvar, *segvar, *bchd;
	extern float *histogramdisplay, *histogramweight;
	extern double *length_weight, *lambda, *nodeinflow;
	extern double *cond, *nodpress, *condsum, *hfactor2sum;
	extern double *sheartarget, *shearfac, *hfactor1, *hfactor2;
	extern double **hmat, **kmat;

	ista = ivector(1, nseg);
	iend = ivector(1, nseg);
	segclasstyp = ivector(1, nseg);
	nodclasstyp = ivector(1, nnod);
	nodout = ivector(1, nnod);
	nodtyp = ivector(1, nnod);
	nodnod = imatrix(1, nodsegm, 1, nnod);
	nodseg = imatrix(1, nodsegm, 1, nnod);
	qq = vector(1, nseg);
	tau = vector(1, nseg);
	segpress = vector(1, nseg);
	bchd = vector(1, nnod);		//larger than needed, but we don't know how many boundary nodes at this point
	nodpress = dvector(1, nnod);
	lambda = dvector(1, nnod);
	cond = dvector(1, nseg);
	condsum = dvector(1, nnod);
	hfactor2sum = dvector(1, nnod);
	nodeinflow = dvector(1, nnod);
	sheartarget = dvector(1, nseg);
	shearfac = dvector(1, nseg);
	hfactor1 = dvector(1, nseg);
	hfactor2 = dvector(1, nseg);
	length_weight = dvector(1, nnod); // determines weighted penalty sum for nodes

	nodelambda = ivector(1, nnod);
	bcnod = ivector(1, nnod);	//larger than needed, but we don't know how many boundary nodes at this point
	flow_direction = ivector(1, nseg);		// determines the direction of flow in a segment (1 or -1)
	known_flow_direction = ivector(1, nseg);
	knowntyp = ivector(1, nnod);
	segvar = vector(1, nseg);
	nodvar = vector(1, nnod);
	nk = ivector(1, nnod);		// From putrank
	nodrank = ivector(1, nnod);
	hmat = dmatrix(0, nodsegm, 1, nnod); // nodsegm = max number of segments allowed per node
	kmat = dmatrix(0, nodsegm, 1, nnod);
	histogramdisplay = vector(1, LMAX(nseg, nnod)); 	// added for histogram display
	histogramweight = vector(1, LMAX(nseg, nnod)); 	// added for histogram display
}
