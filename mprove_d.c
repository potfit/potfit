/**** rewritten for double precision and zero-offset vectors and matrices ****
***** by Peter Brommer, ITAP, 2002-10-10                                  ***/

#define NRANSI
#include "nrutil.h"
void lubksb_d(double **a, int n, int *indx, double b[]);

void mprove_d(double **a, double **alud, int n, int indx[], double b[], double x[])
{
	int j,i;
	double sdp;
	double *r;

	r=dvector(0,n-1); /* zo */
	for (i=0;i<n;i++) { /* zo */
		sdp = -b[i];
		for (j=0;j<n;j++) sdp += a[i][j]*x[j]; /* zo */
		r[i]=sdp;
	}
	lubksb_d(alud,n,indx,r);
	for (i=0;i<n;i++) x[i] -= r[i]; /* zo */
	free_dvector(r,0,n-1); /* zo */
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */
