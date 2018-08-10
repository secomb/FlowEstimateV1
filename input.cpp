/******************************************************
input - reads input files.  For FlowEst2018. TWS August 2018
Note that input files refer to segname and nodname, but
all arrays use iseg and inod, which are assigned when
reading the file, as indices.
*******************************************************/
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void input(void)
{
	extern int nmaxvessel, nmaxtissue, nmax, ktausteps, maxinsideit;
	extern int mxx, myy, mzz, nseg, nnod, nnodbck, nodsegm;
	extern int nitmax, nitmax1, varyviscosity, phaseseparation;
	extern int *bcnodname, *bctyp, *actual_direction;
	extern int *segname, *segtyp, *nodname;
	extern int **segnodname, **nodseg;

	extern float lb, maxl, pi1, alx, aly, alz;
	extern float tol, omega, qtol, hdtol, optw, optlam, constvisc, vplas, mcvcorr, mcv, consthd;
	extern float *bifpar, *cpar, *viscpar, *fahrpr;
	extern float *diam, *q, *hd, *bcprfl, *bchd, *xsl0, *xsl1, *xsl2;
	extern float **cnode;

	extern double presstarget1, sheartarget1, *lseg;

	int i, iseg, max = 100;
	FILE *ifp;
	char bb[100];

	ifp = fopen("Network.dat", "r");		// read network.dat
	if (ifp == NULL) perror("Error opening file");
	else {
		fgets(bb, max, ifp);
		printf("%s\n", bb);
		fscanf(ifp, "%f %f %f%*[^\n]", &alx, &aly, &alz);	//dimensions of box in microns; vertex must be at origin
		fscanf(ifp, "%i %i %i%*[^\n]", &mxx, &myy, &mzz);
		fscanf(ifp, "%f%*[^\n]", &lb);
		fscanf(ifp, "%f%*[^\n]", &maxl);
		fscanf(ifp, "%i%*[^\n]", &nodsegm);
		fscanf(ifp, "%i%*[^\n]", &nseg); 		//number of segments in vessel network
		fgets(bb, max, ifp);
		fgets(bb, max, ifp);
		segname = ivector(1, nseg);
		segtyp = ivector(1, nseg);
		segnodname = imatrix(1, 2, 1, nseg);
		diam = vector(1, nseg);
		q = vector(1, nseg);
		hd = vector(1, nseg);
		lseg = dvector(1, nseg);
		actual_direction = ivector(1, nseg);		// gives actual flow directions (for comparison with flow_direction array)
		for (iseg = 1; iseg <= nseg; iseg++) {		//segment properties: name type nodefrom nodeto diameter flow hematocrit
			fscanf(ifp, "%i %i %i %i %f %f %f%*[^\n]", &segname[iseg], &segtyp[iseg], &segnodname[1][iseg],
				&segnodname[2][iseg], &diam[iseg], &q[iseg], &hd[iseg]);
			if (q[iseg] >= 0.) actual_direction[iseg] = 1;
			else actual_direction[iseg] = -1;
		}
		fscanf(ifp, "%i%*[^\n]", &nnod);		//number of nodes in vessel network
		fgets(bb, max, ifp);
		fgets(bb, max, ifp);
		//coordinates of nodes
		nodname = ivector(1, nnod);
		cnode = matrix(1, 3, 1, nnod);
		for (i = 1; i <= nnod; i++) fscanf(ifp, "%i %f %f %f%*[^\n]", &nodname[i], &cnode[1][i], &cnode[2][i], &cnode[3][i]);
		//boundary nodes **** with known boundary conditions ****
		fscanf(ifp, "%i%*[^\n]", &nnodbck);
		fgets(bb, max, ifp);
		fgets(bb, max, ifp);
		bcnodname = ivector(1, nnodbck);
		bctyp = ivector(1, nnodbck);
		bcprfl = vector(1, nnodbck);
		for (i = 1; i <= nnodbck; i++) fscanf(ifp, "%i %i %f%*[^\n]", &bcnodname[i], &bctyp[i], &bcprfl[i]);
		fclose(ifp);

		//Read parameters for drawing segments and nodes
		xsl0 = vector(1, 3);
		xsl1 = vector(1, 3);
		xsl2 = vector(1, 3);
		ifp = fopen("ContourParams.dat", "r");
		fscanf(ifp, "%f %f %f", &xsl0[1], &xsl0[2], &xsl0[3]);
		fgets(bb, max, ifp);
		fscanf(ifp, "%f %f %f", &xsl1[1], &xsl1[2], &xsl1[3]);
		fgets(bb, max, ifp);
		fscanf(ifp, "%f %f %f", &xsl2[1], &xsl2[2], &xsl2[3]);
		fgets(bb, max, ifp);
		fclose(ifp);

		// Read parameters for bifurcation hematocrit and viscosity (dishem and viscor)
		bifpar = vector(1, 3);
		cpar = vector(1, 4);
		viscpar = vector(1, 6);
		fahrpr = vector(1, 4);
		ifp = fopen("RheolParams.dat", "r");
		fscanf(ifp, "%f %f %f", &bifpar[1], &bifpar[2], &bifpar[3]);
		fgets(bb, max, ifp);
		fscanf(ifp, "%f %f %f %f", &cpar[1], &cpar[2], &cpar[3], &cpar[4]);
		fgets(bb, max, ifp);
		fscanf(ifp, "%f %f %f %f %f %f", &viscpar[1], &viscpar[2], &viscpar[3], &viscpar[4], &viscpar[5], &viscpar[6]);
		fgets(bb, max, ifp);
		fscanf(ifp, "%f %f %f %f", &fahrpr[1], &fahrpr[2], &fahrpr[3], &fahrpr[4]);
		fgets(bb, max, ifp);
		fscanf(ifp, "%i %f %f", &nitmax, &tol, &omega);
		fgets(bb, max, ifp);
		fscanf(ifp, "%i %f %f", &nitmax1, &qtol, &hdtol);
		fgets(bb, max, ifp);
		fscanf(ifp, "%f %f", &optw, &optlam);
		fgets(bb, max, ifp);
		fscanf(ifp, "%f %f %f", &constvisc, &vplas, &mcv);
		fgets(bb, max, ifp);
		fscanf(ifp, "%f", &consthd);
		fgets(bb, max, ifp);
		fscanf(ifp, "%i", &varyviscosity);
		fgets(bb, max, ifp);
		fscanf(ifp, "%i", &phaseseparation);
		fgets(bb, max, ifp);
		fclose(ifp);
		mcvcorr = pow(92. / mcv, 0.33333);
		if (phaseseparation == 0) for (i = 1; i <= nseg; i++) hd[i] = consthd;

		ifp = fopen("FlowEstParams.dat", "r");
		fgets(bb, max, ifp);
		fscanf(ifp, "%lf %*[^\n]", &presstarget1);
		fscanf(ifp, "%lf %*[^\n]", &sheartarget1);
		fscanf(ifp, "%i %*[^\n]", &ktausteps);
		fscanf(ifp, "%i %*[^\n]", &maxinsideit);
		fclose(ifp);
	}
}