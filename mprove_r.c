/**** rewritten for double precision and zero-offset vectors and matrices ****
***** by Peter Brommer, ITAP, 2002-10-10                                  ****
***** adapted to real variables (ITAP standard) by PB, ITAP, 2002-10-24   ***/


#define NRANSI
#include "potfit.h"
#include "nrutil_r.h"

void lubksb_r(real **a, int n, int *indx, real b[]);

void mprove_r(real **a, real **alud, int n, int indx[], real b[], real x[])
{
	int j,i;
	real sdp;
	real *r;

	r=dvector(0,n-1); /* zo */
	for (i=0;i<n;i++) { /* zo */
		sdp = -b[i];
		for (j=0;j<n;j++) sdp += a[i][j]*x[j]; /* zo */
		r[i]=sdp;
	}
	lubksb_r(alud,n,indx,r);
	for (i=0;i<n;i++) x[i] -= r[i]; /* zo */
	free_dvector(r,0,n-1); /* zo */
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */
