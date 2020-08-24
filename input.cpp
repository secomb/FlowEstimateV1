/******************************************************
input - reads input files.  For FlowEst2018. TWS August 2018
Note that input files refer to segname and nodname, but
all arrays use iseg and inod, which are assigned when
reading the file, as indices.
------------------------------------------------------
Sample FlowEstParams
1		solvetype (1 = LU, 2 = CG, 3 = Relaxation)
31.		target pressure (mmHg)
0		varytargetshear (0 or 1)
63.5	target shear stress for fixed target (dyn/cm2)
0.001	ktaustart, starting value of ktau
10		number of ktau steps, doubling each time
5		max number of runs for each tau value
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
	extern int nmaxvessel, nmaxtissue, nmax, ktausteps, maxinsideit, solvetyp;
	extern int mxx, myy, mzz, nseg, nnod, nnodbck, nodsegm, seed;
	extern int nitmax, nitmax1, varyviscosity, phaseseparation, varytargetshear;
	extern int nitmax2, nitmax3;
	extern int *bcnodname, *bctyp, *actual_direction;
	extern int *segname, *segtyp, *nodname;
	extern int **segnodname, **nodseg;

	extern float lb, maxl, pi1, alx, aly, alz, known_flow_weight;
	extern float tol, omega, qtol, hdtol, optw, optlam, constvisc, vplas, mcvcorr, mcv, consthd, diamcrit, diamcorr;
	extern float *bifpar, *cpar, *viscpar, *fahrpr;
	extern float *diam, *q, *hd, *bcprfl, *bchd, *xsl0, *xsl1, *xsl2;
	extern float **cnode;
    
    extern float *tmpordart, *tmpordvein;

	extern double *lseg, presstarget1, sheartarget1, ktaustart, eps, omega1, omega2;

	int i, iseg, max = 200;
	FILE *ifp;
	char bb[200];

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

        tmpordart = vector(1, nseg);
        tmpordvein = vector(1, nseg);
		for (iseg = 1; iseg <= nseg; iseg++) {		//segment properties: name type nodefrom nodeto diameter flow hematocrit
			fscanf(ifp, "%i %i %i %i %f %f %f %f %f%*[^\n]", &segname[iseg], &segtyp[iseg], &segnodname[1][iseg],
				&segnodname[2][iseg], &diam[iseg], &q[iseg], &hd[iseg], &tmpordart[iseg], &tmpordvein[iseg]);
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
		fscanf(ifp, "%f %f %f%*[^\n]", &xsl0[1], &xsl0[2], &xsl0[3]);
		fscanf(ifp, "%f %f %f%*[^\n]", &xsl1[1], &xsl1[2], &xsl1[3]);
		fscanf(ifp, "%f %f %f%*[^\n]", &xsl2[1], &xsl2[2], &xsl2[3]);
		fclose(ifp);

		// Read parameters for bifurcation hematocrit and viscosity (dishem and viscor)
		bifpar = vector(1, 3);
		cpar = vector(1, 4);
		viscpar = vector(1, 6);
		fahrpr = vector(1, 4);
		ifp = fopen("RheolParams.dat", "r");
		fscanf(ifp, "%f %f %f%*[^\n]", &bifpar[1], &bifpar[2], &bifpar[3]);
		fscanf(ifp, "%f %f %f %f%*[^\n]", &cpar[1], &cpar[2], &cpar[3], &cpar[4]);
		fscanf(ifp, "%f %f %f %f %f %f%*[^\n]", &viscpar[1], &viscpar[2], &viscpar[3], &viscpar[4], &viscpar[5], &viscpar[6]);
		fscanf(ifp, "%f %f %f %f%*[^\n]", &fahrpr[1], &fahrpr[2], &fahrpr[3], &fahrpr[4]);
		fscanf(ifp, "%i %f %f%*[^\n]", &nitmax, &tol, &omega);	//not used
		fscanf(ifp, "%i %f %f%*[^\n]", &nitmax1, &qtol, &hdtol);	//not used
		fscanf(ifp, "%f %f%*[^\n]", &optw, &optlam);
		fscanf(ifp, "%f %f %f%*[^\n]", &constvisc, &vplas, &mcv);
		fscanf(ifp, "%f%*[^\n]", &consthd);
		fscanf(ifp, "%i%*[^\n]", &varyviscosity);
		fscanf(ifp, "%i%*[^\n]", &phaseseparation);
		fclose(ifp);
		mcvcorr = pow(92. / mcv, 0.33333);
		if (phaseseparation == 0) for (i = 1; i <= nseg; i++) hd[i] = consthd;

		ifp = fopen("FlowEstParams.dat", "r");
		fgets(bb, max, ifp);
		fscanf(ifp, "%i %*[^\n]", &solvetyp);
		fscanf(ifp, "%lf %*[^\n]", &presstarget1);
		fscanf(ifp, "%i %*[^\n]", &varytargetshear);
		fscanf(ifp, "%lf %*[^\n]", &sheartarget1);
		fscanf(ifp, "%lf %*[^\n]", &ktaustart);
		fscanf(ifp, "%i %*[^\n]", &ktausteps);
		fscanf(ifp, "%i %*[^\n]", &maxinsideit);
		fscanf(ifp, "%lf %*[^\n]", &eps);
		fscanf(ifp, "%i %*[^\n]", &nitmax2);
		fscanf(ifp, "%i %*[^\n]", &nitmax3);
		fscanf(ifp, "%lf %lf %*[^\n]", &omega1, &omega2);
		fscanf(ifp, "%f %*[^\n]", &diamcrit);
		fscanf(ifp, "%f %*[^\n]", &known_flow_weight);
		fscanf(ifp, "%i %*[^\n]", &seed);
        fscanf(ifp, "%f %*[^\n]", &diamcorr);
		fclose(ifp);
	}
}
