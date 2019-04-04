/************************************************************************
writeflow - for FlowEst2018a
rewrite network.dat with updated flows
TWS 2018
Pressures in mmHg, flows in um^3/s, viscosity in cP
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void writeflow()
{
	extern int nseg, nnod, nnodbc, *segname, *segtyp, *nodname;
	extern int *bcnod, *knowntyp, **segnodname, **nodseg, *ista;
	extern float *diam, *q, *hd, *bcprfl;
	extern double *nodpress;

	int i, iseg, inod, inodbc, max = 200, type = 2;
	float qinput, hdinput;
	FILE *ifp, *ofp;
	char bb[200];

	//network data file
	ifp = fopen("Network.dat", "r");
	ofp = fopen("NetworkNew.dat", "w");
	for (i = 1; i <= 8; i++) {
		fgets(bb, max, ifp);
		fprintf(ofp, "%s", bb);
	}
	for (iseg = 1; iseg <= nseg; iseg++) {
		fscanf(ifp, "%i %i %i %i %f %f %f",
			&segname[iseg], &segtyp[iseg], &segnodname[1][iseg], &segnodname[2][iseg], &diam[iseg], &qinput, &hdinput);
		fgets(bb, max, ifp);
		fprintf(ofp, "%i %i %i %i %f %f %f",
			segname[iseg], segtyp[iseg], segnodname[1][iseg], segnodname[2][iseg], diam[iseg], q[iseg], hd[iseg]);
		fprintf(ofp, "%s", bb);
	}
	for (inod = 1; inod <= nnod+2; inod++) {	//copy 2 header lines also
		fgets(bb, max, ifp);
		fprintf(ofp, "%s", bb);
	}
	fprintf(ofp, "%i Number of boundary nodes\n", nnodbc);
	fprintf(ofp, "Node Bctype Press / Flow\n");
	for (inodbc = 1; inodbc <= nnodbc; inodbc++) {
		inod = bcnod[inodbc];
		iseg = nodseg[1][inod];
		if (inod == ista[iseg]) {
			fprintf(ofp, "%i %i %f 0.400 50.000\n", nodname[inod], knowntyp[inod], q[iseg]);
		}
		else {
			fprintf(ofp, "%i %i %f 0.400 50.000\n", nodname[inod], knowntyp[inod], -q[iseg]);
		}
	}
	fclose(ifp);
	fclose(ofp);
}