/****************************************************************
* 
*  brent_r.c: Minmization of a multivariable function according to 
*      Brent's algorithm.
*
*****************************************************************/
/****************************************************************
* $Revision: 1.5 $
* $Date: 2003/03/19 09:05:32 $
*****************************************************************/
/**** rewritten for double precision                                      ****
***** by Peter Brommer, ITAP, 2002-10-10                                  ***/
/**** Adapted to Powell requirements 2002-10-11 			  ***/
/**** adapted to real variables (ITAP standard) by PB, ITAP, 2002-10-24   ***/

#include <math.h>
#define NRANSI
#include "potfit.h"
#include "nrutil_r.h"
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

extern real *xicom,*delcom;

real brent_r(real ax, real bx, real cx, real fbx, real tol, real *xmin,
		real *xmin2, real *fxmin, real *fxmin2)
/* take bracket (a,b,c), f(b), tol, pointers to xmin, xmin2, vectors fxmin, fxmin2 */
{
	int iter,j;
	real a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	real e=0.0;   /*Distance moved on step before last */
	static real *vecu = NULL, *fxu = NULL;   /* Vector of location u */

	if (fxu  == NULL) fxu  = dvector(0,mdim-1);
	if (vecu == NULL) vecu = dvector(0,ndimtot-1);
	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=fbx;
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			*xmin2=w;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		for (j=0;j<ndimtot;j++) vecu[j]=xicom[j]+u*delcom[j];/*set vecu*/
		fu=(*calc_forces)(vecu,fxu);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			for (j=0;j<mdim;j++) {      /* shifting fxmin2, */
						    /* fxmin, fxu */
				fxmin2[j]=fxmin[j];
				fxmin[j]=fxu[j];
			}		
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				for (j=0;j<mdim;j++) fxmin2[j]=fxu[j]; 
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	nrerror("Too many iterations in brent");
	*xmin=x;
	*xmin2=w;
	return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */
