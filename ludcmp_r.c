/****************************************************************
*
*  ludcmp_r.c: Routines used for LU-decomposition of a matrix
*
*****************************************************************/
/****************************************************************
* $Revision: 1.3 $
* $Date: 2003/03/19 09:05:36 $
*****************************************************************/
/**** rewritten for double precision and zero-offset vectors and matrices ****
***** by Peter Brommer, ITAP, 2002-10-10                                  ***/
/**** adapted to real variables (ITAP standard) by PB, ITAP, 2002-10-24   ***/

#include <math.h>
#define NRANSI
#include "potfit.h"
#include "nrutil_r.h"
#define TINY 1.0e-20;

void ludcmp_r(real **a, int n, int *indx, real *d)
{
	int i,imax,j,k;
	real big,dum,sum,temp;
	real *vv;

	vv=dvector(0,n-1);  /* zo */
	*d=1.0;
	for (i=0;i<n;i++) { /* zo */
		big=0.0;
		for (j=0;j<n;j++) /* zo */
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) { /* zo */
		for (i=0;i<j;i++) { /* zo */
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j]; /* zo */
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {   /* zo */
			sum=a[i][j];
			for (k=0;k<j;k++)  /* zo */
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {  /* zo */
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum; /* zo */
		}
	}
	free_dvector(vv,0,n-1);
}
#undef TINY
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */
