/**** rewritten for double precision and zero-offset vectors and matrices ****
***** adapted to Powell requrirements (return vector instead of value)...
***** by Peter Brommer, ITAP, 2002-10-10                                  ***/


#define NRANSI
#include <math.h>
#include "nrutil.h"
#define TOL 2.0e-4

int ncom,mcom;
double *pcom,*xicom;
double (*nrfunc)(double [], double [], int, int);

void mnbrak_d(double *ax, double *bx, double *cx, double *fa, double *fb, 
	double *fc, double (*func)(double));
double brent_d(double ax, double bx, double cx, double fbx, double tol, 
	double *xmin, double *xmin2, double *fxmin, double *fxmin2);
double f1dim_d(double x);

double linmin_d(double p[], double xi[], double fxi1, int n, int m, double *x1, 
		double *x2, double *fret1, double *fret2, 
		double (*func)(double [], double [], int ,int))
/* takes vector xi (direction of search), p (originating point), n,m (dimensions),
   x, x2 (two best locations),
   fret1, fret2 (return vectors) and the function to minimize with input and result
   vectors as arguments */
	{
	int j;
	double xx,fx,fb,bx,ax;
	double fa=fxi1;
	double xmin;
	double xmin2;
	ncom=n;
	mcom=m;
	pcom=dvector(0,n-1);
	xicom=dvector(0,n-1);
	nrfunc=func;
	for (j=0;j<n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;                 /*do not change without correcting fa*/
	xx=1.0;
	mnbrak_d(&ax,&xx,&bx,&fa,&fx,&fb,f1dim_d);
	/*only do brent if there's something to improve*/
/*	if ((fabs((fa-fx)/fx)>TOL/100.)||(fabs((fb-fx)/fx))>TOL/100.) */
		fx=brent_d(ax,xx,bx,fx,TOL,&xmin,&xmin2,fret1,fret2);
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
