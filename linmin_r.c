/**** rewritten for double precision and zero-offset vectors and matrices ****
***** adapted to Powell requrirements (return vector instead of value)...
***** adapted to real variables (ITAP standard) by PB, ITAP, 2002-10-24   ***/
/**** by Peter Brommer, ITAP, 2002-10-10                                  ***/


#define NRANSI
#include <math.h>
#include "potfit.h"
#include "nrutil_r.h"
#define TOL 2.0e-4

real *xicom,*delcom;

void mnbrak_r(real *ax, real *bx, real *cx, real *fa, real *fb, 
	real *fc, real (*func)(real));
real brent_r(real ax, real bx, real cx, real fbx, real tol, 
	real *xmin, real *xmin2, real *fxmin, real *fxmin2);
real f1dim_r(real x);

real linmin_r(real xi[], real del[], real fxi1, int n, int m, real *x1, 
		real *x2, real *fret1, real *fret2 )
/* takes vector xi (direction of search), p (originating point), 
   n,m (dimensions),
   x, x2 (two best locations),
   fret1, fret2 (return vectors) as arguments */
	{
	int j;
	real xx,fx,fb,bx,ax;
	real fa=fxi1;
	real xmin;
	real xmin2;
        xicom=xi;
        delcom=del;
	ax=0.0;                 /*do not change without correcting fa,*/
				/*saves 1 fcalc...*/
	xx=1.0;
	mnbrak_r(&ax,&xx,&bx,&fa,&fx,&fb,f1dim_r);
	fx=brent_r(ax,xx,bx,fx,TOL,&xmin,&xmin2,fret1,fret2);
	for (j=0;j<ndimtot;j++) {
		del[j] *= xmin;
		xi[j] += del[j];
	}
	*x1=xmin;
	*x2=xmin2;
	return fx;
}
#undef TOL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */

