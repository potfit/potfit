/**** rewritten for double precision and zero-offset vectors and matrices ****
***** adapted to Powell requrirements (return vector instead of value)...
***** adapted to real variables (ITAP standard) by PB, ITAP, 2002-10-24   ***/
/**** by Peter Brommer, ITAP, 2002-10-10                                  ***/


#define NRANSI
#include <math.h>
#include "potfit.h"
#include "nrutil_r.h"
#define TOL 2.0e-4

int ncom,mcom;
real *pcom,*xicom;
real (*nrfunc)(real [], real []);

void mnbrak_r(real *ax, real *bx, real *cx, real *fa, real *fb, 
	real *fc, real (*func)(real));
real brent_r(real ax, real bx, real cx, real fbx, real tol, 
	real *xmin, real *xmin2, real *fxmin, real *fxmin2);
real f1dim_r(real x);

real linmin_r(real p[], real xi[], real fxi1, int n, int m, real *x1, 
		real *x2, real *fret1, real *fret2, 
		real (*func)(real [], real []))
/* takes vector xi (direction of search), p (originating point), 
   n,m (dimensions),
   x, x2 (two best locations),
   fret1, fret2 (return vectors) and the function to minimize as arguments */
	{
	int j;
	real xx,fx,fb,bx,ax;
	real fa=fxi1;
	real xmin;
	real xmin2;
	ncom=n;
	mcom=m;
	pcom=dvector(0,n-1);
	xicom=dvector(0,n-1);
	nrfunc=func;
	for (j=0;j<n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;                 /*do not change without correcting fa,*/
				/*saves 1 fcalc...*/
	xx=1.0;
	mnbrak_r(&ax,&xx,&bx,&fa,&fx,&fb,f1dim_r);
		fx=brent_r(ax,xx,bx,fx,TOL,&xmin,&xmin2,fret1,fret2);
	for (j=0;j<n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	*x1=xmin;
	*x2=xmin2;
	free_dvector(xicom,0,n-1);
	free_dvector(pcom,0,n-1);
	return fx;
}
#undef TOL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */
