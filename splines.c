/****************************************************************
 *
 *  spline.c: Contains all routines used for spline interpolation
 *      with equidistant or non-equidistant sampling points.
 *
 *****************************************************************/

/****************************************************************
* $Revision: 1.12 $
* $Date: 2004/08/16 13:02:51 $
*****************************************************************/


/******************************************************************************
 *
 * spline_ed: initializes second derivatives used for spline interpolation
 *            (equidistant x[i])
 *
 *****************************************************************************/



/*** rewritten for real variables (ITAP standard) and zero-offset vectors ****
 *** adapted to equidistant x                                             ****
 *** by Peter Brommer, ITAP 2002-11-27                                    ****/
#define NRANSI
#ifndef POTSCALE
#include "potfit.h"
#else
#include "potscale.h"
#endif
#include "nrutil_r.h"

void spline_ed(real xstep, real y[], int n, real yp1, real ypn, real y2[])
{
	int i,k;
      	real p,qn,un;
	static real *u=NULL;
        static int  nmax=0;

	if (n>nmax) {
	    u = (real *) realloc(u, (n-1)*sizeof(real));
	    nmax=n;
	}
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(xstep))*((y[1]-y[0])/(xstep)-yp1);
	}
	for (i=1;i<n-1;i++) {
	    /* sig=(x[i]-x[i-1])/(x[i+1]-x[i-1])=.5;*/
		p=0.5*y2[i-1]+2.0;
		y2[i]=(-0.5)/p;
		u[i]=(y[i+1]-y[i])/xstep - (y[i]-y[i-1])/(xstep);
		u[i]=(6.0*u[i]/(2*xstep)-0.5*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(xstep))*(ypn-(y[n-1]-y[n-2])/(xstep));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */

/******************************************************************************
 *
 * splint_ed: interpolates the function with splines 
 *            (equidistant x[i]) 
 *
 *****************************************************************************/

/**** rewritten for real variables (ITAP standard) and zero-offset vectors ****
 **** adapted to equidistant x[i]                                          ****
 **** Peter Brommer, ITAP, 2002-11-27                                      ***/

real splint_ed(pot_table_t *pt, real *xi, int col, real r)
{
  real a, b, istep, rr, p1, p2, d21, d22;
  int k;

  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0) {printf("%f %f %d\n",r,pt->begin[col],col);error("short distance");}

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (rr * istep);
  b     = (rr - k * pt->step[col]) * istep;
  k    += pt->first[col];
  a     = 1.0 - b;
  p1    = xi[k];
  d21   = pt->d2tab[k++];
  p2    = xi[k];
  d22   = pt->d2tab[k];
  

  return a * p1 + b * p2 +
      ((a*a*a - a) * d21 + (b*b*b - b) * d22) / (6.0 * istep * istep);
}
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */

/******************************************************************************
 *
 * splint_comb_ed: calculates spline interpolation of a function (return value)
 *            and its gradiend (grad), equidistant x[i]
 *
 *****************************************************************************/

real splint_comb_ed(pot_table_t *pt, real *xi, int col, real r, real *grad)
{
  real a, b, istep, rr, p1, p2, d21, d22;
 int k;

  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0) error("short distance! in splint_comb_ed");

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (rr * istep);
  b     = (rr - k * pt->step[col]) * istep;
  k    += pt->first[col];
  a     = 1.0 - b;
  p1    = xi[k];
  d21   = pt->d2tab[k++];
  p2    = xi[k];
  d22   = pt->d2tab[k];
  *grad = (p2-p1)*istep + 
      ((3*(b*b) - 1) * d22 - (3*(a*a) - 1) * d21) / (6.0 * istep);
  return a * p1 + b * p2 +
      ((a*a*a - a) * d21 + (b*b*b - b) * d22) / (6.0 * istep * istep);
}

/******************************************************************************
 *
 * splint_grad_ed: calculates the first derivative from spline interpolation
 *            (equidistant x[i])
 *
 *****************************************************************************/


real splint_grad_ed(pot_table_t *pt, real *xi, int col, real r)
{
  real a, b, istep, rr, p1, p2, d21, d22;
  int k;

  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0) error("short distance! in splint_grad_ed");

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (rr * istep);
  b     = (rr - k * pt->step[col]) * istep;
  k    += pt->first[col];
  /* Check if we are at the last index */
  if (k>=pt->last[col]) {
    k--;b+=1.0;
  }
  a     = 1.0 - b;
  p1    = xi[k];
  d21   = pt->d2tab[k++];
  p2    = xi[k];
  d22   = pt->d2tab[k];

  return (p2-p1)*istep + 
      ((3*(b*b) - 1) * d22 - (3*(a*a) - 1) * d21) / (6.0 * istep);

}

/******************************************************************************
 *
 * splint_dir: interpolates the function with splines 
 *            (equidistant AND NON-eq.dist x[i]) 
 *            with known index position
 *
 *****************************************************************************/

/**** rewritten for real variables (ITAP standard) and zero-offset vectors ****
 **** adapted to equidistant x[i]                                          ****
 **** Peter Brommer, ITAP, 2002-11-27                                      ***/
real splint_dir(pot_table_t *pt, real *xi, int col, int k, real b, real step)
{
  real a, istep, p1, p2, d21, d22;



  /* indices into potential table */
  a     = 1.0 - b;
  p1    = xi[k];
  d21   = pt->d2tab[k++];
  p2    = xi[k];
  d22   = pt->d2tab[k];
  

  return a * p1 + b * p2 +
    ((a*a*a - a) * d21 + (b*b*b - b) * d22) * (step * step) / 6.0 ;
}
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */

/****************************************************************************
 *
 * splint_comb_dir_ed: calculates spline interpolation of a function 
 *            (return value)
 *            and its gradiend (grad), equidistant and non-eqd x[i]
 *            with known index position
 *
 *****************************************************************************/

real splint_comb_dir(pot_table_t *pt, real *xi, int col, int k, real b, real step, real *grad)
{
  real a, istep, p1, p2, d21, d22;


  /* indices into potential table */

  a     = 1.0 - b;
  p1    = xi[k];
  d21   = pt->d2tab[k++];
  p2    = xi[k];
  d22   = pt->d2tab[k];
  *grad = (p2-p1)/step + 
      ((3*(b*b) - 1) * d22 - (3*(a*a) - 1) * d21) * step / 6.0;
  return a * p1 + b * p2 +
      ((a*a*a - a) * d21 + (b*b*b - b) * d22) *(step*step) / 6.0 ;
}

/******************************************************************************
 *
 * splint_grad_dir: calculates the first derivative 
 *            from spline interpolation (equidistant and NON-eqd. x[i])
 *            with known index position
 *
 *****************************************************************************/


real splint_grad_dir(pot_table_t *pt, real *xi, int col, int k, real b, real step)
{
  real a, istep, p1, p2, d21, d22;
  
  /* indices into potential table */
  a     = 1.0 - b;
  p1    = xi[k];
  d21   = pt->d2tab[k++];
  p2    = xi[k];
  d22   = pt->d2tab[k];

  return (p2-p1)/step + 
      ((3*(b*b) - 1) * d22 - (3*(a*a) - 1) * d21) * step / 6.0 ;

}

/******************************************************************************
 *
 * spline_ne  : initializes second derivatives used for spline interpolation
 *            (nonequidistant x[i])
 *
 *****************************************************************************/


/*** rewritten for real variables (ITAP standard) and zero-offset vectors ****
 *** by Peter Brommer, ITAP 2002-11-27                                    ****/


#define NRANSI

void spline_ne(real x[], real y[], int n, real yp1, real ypn, real y2[])
{
	int i,k;
	real p,qn,sig,un;
	static real *u=NULL;
	static int  nmax=0;

	if (n>nmax) {
	    u = (real *) realloc(u, (n-1)*sizeof(real));
	    nmax=n;
	}
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}
#undef NRANSI

/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */


/******************************************************************************
 *
 * splint_ne: interpolates the function with splines 
 *            (nonequidistant x[i])
 *
 *****************************************************************************/

/**** rewritten for real variables (ITAP standard)                         ****
 **** zero-offset-vectors                                                  ****
 **** Peter Brommer, ITAP, 2002-11-27                                      ***/

real splint_ne(pot_table_t *pt, real *xi, int col, real r)
{
/* 	void nrerror(char error_text[]); */
  int klo,khi,k;
  real h,b,a,rr,d22,d21,p1,p2,x1,x2;
  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0) {printf("%f %f %d\n",r,pt->begin[col],col);error("short distance");}
  
  klo=pt->first[col];
  khi=pt->last[col];

  /* Find index by bisection */
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (pt->xcoord[k] > r) khi=k;
    else klo=k;
  }
  x1    =pt->xcoord[klo];
  x2    =pt->xcoord[khi];
  h     =x2-x1;
  p1    =xi[klo];
  p2    =xi[khi];
  d21   =pt->d2tab[klo];
  d22   =pt->d2tab[khi];

/*   if (h == 0.0) error("Bad xa input to routine splint"); */
  b=(r-x1)/h;
  a=(1.-b);
  return a*p1+b*p2+
    ((a*a*a-a)*d21+(b*b*b-b)*d22)*(h*h)/6.0;

}
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */

/******************************************************************************
 *
 * splint_comb_ne: calculates spline interpolation of a function (return value)
 *            and its gradiend (grad), non-equidistant x[i]
 *
 *****************************************************************************/

real splint_comb_ne(pot_table_t *pt, real *xi, int col, real r, real *grad)
{
  int klo,khi,k;
  real h,b,a,rr,d22,d21,p1,p2,x1,x2;
  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0) error("short distance");
  
  klo=pt->first[col];
  khi=pt->last[col];

  /* Find index by bisection */
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (pt->xcoord[k] > r) khi=k;
    else klo=k;
  }
  x1    =pt->xcoord[klo];
  x2    =pt->xcoord[khi];
  h     =x2-x1;
  p1    =xi[klo];
  p2    =xi[khi];
  d21   =pt->d2tab[klo];
  d22   =pt->d2tab[khi];

  b=(r-x1)/h;
  a=(1.-b);

  *grad = (p2-p1)/h + 
      ((3*(b*b) - 1) * d22 - (3*(a*a) - 1) * d21)* h  / 6.0;

  return a*p1+b*p2+
    ((a*a*a-a)*d21+(b*b*b-b)*d22)*(h*h)/6.0;

}

/******************************************************************************
 *
 * splint_grad_ne: calculates the first derivative from spline interpolation
 *            (equidistant x[i])
 *
 *****************************************************************************/

real splint_grad_ne(pot_table_t *pt, real *xi, int col, real r)
{
  int klo,khi,k;
  real h,b,a,rr,d22,d21,p1,p2,x1,x2;
  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0) error("short distance");
  
  klo=pt->first[col];
  khi=pt->last[col];

  /* Find index by bisection */
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (pt->xcoord[k] > r) khi=k;
    else klo=k;
  }
  x1    =pt->xcoord[klo];
  x2    =pt->xcoord[khi];
  h     =x2-x1;
  p1    =xi[klo];
  p2    =xi[khi];
  d21   =pt->d2tab[klo];
  d22   =pt->d2tab[khi];

  b=(r-x1)/h;
  a=(1.-b);

  return (p2-p1)/h + 
      ((3*(b*b) - 1) * d22 - (3*(a*a) - 1) * d21)* h  / 6.0;

}

