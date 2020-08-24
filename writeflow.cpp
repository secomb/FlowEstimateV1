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
	extern int nodsegm;
	extern int nseg, nnod, nnodbc, *segname, *segtyp, *nodname, *ista;
	extern int *bcnod, *knowntyp, **segnodname, **nodseg;
	extern float *diam, *q, *hd, *bcprfl;
	extern float **cnode;
	extern double *nodpress;

	int i, iseg, inod, inodbc, max = 200;
	FILE *ifp, *ofp;
	char bb[200];

	ifp = fopen("Network.dat", "r");
	ofp = fopen("NetworkNew.dat", "w");
	for (i = 1; i <= 6; i++) {
		fgets(bb, max, ifp);
		fprintf(ofp, "%s", bb);
	}
	fprintf(ofp, "%i total number of segments\n", nseg); 		//updated number of segments in vessel network
	fprintf(ofp, "name type from  to  diam flow HD\n");
	for (iseg = 1; iseg <= nseg; iseg++) fprintf(ofp, "%i %i %i %i %7.3f %7.3f %7.4f\n", //updated type, flow and hematocrit for each segment
			segname[iseg], segtyp[iseg], segnodname[1][iseg], segnodname[2][iseg], diam[iseg], q[iseg], hd[iseg]);
	fprintf(ofp, "%i total number of nodes\n", nnod);		//updated number of nodes in vessel network
	fprintf(ofp, "name    x       y       z\n");
	for (i = 1; i <= nnod; i++) fprintf(ofp, "%i %g %g %g\n", nodname[i], cnode[1][i], cnode[2][i], cnode[3][i]);
	fprintf(ofp, "%i Number of boundary nodes\n", nnodbc);
	fprintf(ofp, "Node Bctype Press / Flow\n");
	for (inodbc = 1; inodbc <= nnodbc; inodbc++) {
		inod = bcnod[inodbc];
		iseg = nodseg[1][inod];
		//node that classification codes 4 to 9 are placed in the Bctype field
		if (inod == ista[iseg]) fprintf(ofp, "%i %i %f 0.400 50.000\n", nodname[inod], knowntyp[inod], q[iseg]);
		else fprintf(ofp, "%i %i %f 0.400 50.000\n", nodname[inod], knowntyp[inod], -q[iseg]);
	}
	fclose(ifp);
	fclose(ofp);
}
