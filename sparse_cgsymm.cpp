/************************************************************************
cgsymm - for FlowEst2018a  TWS August 2018
Conjugate gradient method for symmetric system of linear equations:  a(i,j)*x(j) = b(i)
Modified to use external matrix multiplier
x is vector of [p lambda]
Method uses following code from wikipedia:
-------------------------------------
function [x] = conjgrad(a,b,x)
r=b-a*x;
p=r;
rsold=r'*r;
for i=1:size(a)(1)
v=a*p;
alpha=(r'*r)/(p'*v);
x=x+alpha*p;
r=r-alpha*v;
rsnew=r'*r;
if sqrt(rsnew)<1e-10
break;
end
p=r+rsnew/rsold*p;
rsold=rsnew;
end
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void Amultiply(double *input, double *output);

double dot(double *a, double *b, int n)
{
	double x = 0.;
	int i;
	for (i = 1; i <= n; i++) x += a[i] * b[i];
	return x;
}

double sparse_cgsymm(double *b, double *x, int n, double eps, int itmax)
{
	int i, jj;
	double *r, *p, *v, rsold, rsnew, pv;

	r = dvector(1, n);
	v = dvector(1, n);
	p = dvector(1, n);

	Amultiply(x, r);	//r=Ax
	for (i = 1; i <= n; i++) {
		r[i] = b[i] - r[i];
		p[i] = r[i]; //  initialize p as r
	}
	rsold = dot(r, r, n);
	jj = 0;
	do {
		Amultiply(p, v);	//v=Ap
		pv = dot(p, v, n);
		for (i = 1; i <= n; i++) {
			x[i] += rsold / pv * p[i];
			r[i] -= rsold / pv * v[i];
		}
		rsnew = dot(r, r, n);
		for (i = 1; i <= n; i++) p[i] = r[i] + rsnew / rsold * p[i];
		jj++;
		rsold = rsnew;
		//if (jj % 1000 == 0) printf("cgsymm: %i %e\n", jj, rsnew);
	} while (rsnew > n*eps*eps && jj < itmax);
	printf("cgsymm: %i %e\n", jj, rsnew);
	free_dvector(r, 1, n);
	free_dvector(v, 1, n);
	free_dvector(p, 1, n);

	return rsnew;
}