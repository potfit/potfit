/******************************************************************************
* an example to minimize Wood's function using geodesiclm
*
* Author: Mingjian Wen (wenxx151@umn.edu), University of Minnesota
******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dtat_f(i,j) dtat[j][i]

/******************************************************************************
* function to compute force Wood's function 
******************************************************************************/
void func_(int *m, int* n, double* x, double* fvec)
{
  int i;
  
  fvec[0]=10.0*(x[1] - x[0]*x[0]);
  fvec[1]=1.0 - x[0];
  fvec[2]=sqrt(90.0)*(x[3] - x[2]*x[2]);
  fvec[3]=1.0 - x[2];
  fvec[4]=sqrt(10.0)*(x[1]+x[3] - 2.0);
  fvec[5]=(x[1] - x[3])/sqrt(10.0);

  /* minus reference force */
  for (i=0; i<*m; ++i)
    fvec[i] = fvec[i] - 0.0; 

}

/* jacobian */
void jacobian_(int* m, int* n, double* x,double* fjac){}
/* Avv */
void Avv_(int* m, int* n, double* x, double* v, double* acc){}
/* callback */
void callback_(int* m, int* n, double* x,double* v,double* a,double* fvec,
               double* fjac, double* acc, double* lam, double* dtd,
               double* fvec_new, int* accepted, int* info){}



/******************************************************************************
* main function 
******************************************************************************/

int main()
{

  /* variables */
  int i;

  /* input */
  int n, m;
  double *x;
  int analytic_jac, analytic_Avv, center_diff;
  double h1, h2;
  double** dtd;
  int damp_mode;
  int maxiter, maxfev, maxaev, maxjev;
  double  maxlam, minlam, artol, Cgoal, gtol, xtol, xrtol, ftol, frtol;
  int print_level, print_unit;
  int iaccel, ibroyden, ibold, imethod;
  double avmax, initialfactor, factoraccept, factorreject;

  /* output */
  double* fvec;
  double** fjac;
  int info;
  int niters, nfev, njev, naev;
  int converged;



/**************************************
* set parameters
**************************************/

  n = 4;
  m = 6;
 
  /* allocate memory */
  x = (double*)malloc(n*sizeof(double));
  /* we need to allocate memory, even though it is output */
  fvec = (double*)malloc(m*sizeof(double));
  dtd = (double**)malloc(n*sizeof(double*));
  dtd[0] = (double*)malloc(n*n*sizeof(double)); 
  for(i=1; i<n; ++i)
    dtd[i] = dtd[i-1] + n;
  fjac = (double**)malloc(m*sizeof(double*));
  fjac[0] = (double*)malloc(m*n*sizeof(double)); 
  for(i=1; i<m; ++i)
    fjac[i] = fjac[i-1] + n;


  /* get the optimizable parameters */
    x[0] = -3.0;
    x[1] = -1.0;
    x[2] = -3.0;
    x[3] = -1.0;
  
  analytic_jac=0;
  analytic_Avv=0;
  center_diff=0;  
  
  h1 = 1.0E-5;
  h2 = 0.2;
  factoraccept = 5;
  factorreject = 2;
  maxlam = 1.E7;
  
  imethod = 0;
  initialfactor = 1;
  damp_mode = 0;
  maxiter = 500;
  maxfev = 0;
  Cgoal = 1E-5;
  maxjev = 0;
  iaccel = 0;
  avmax = 0.75;
  maxaev = 0;
  print_level =5;
  print_unit = 6;
  ibold = 0;
  ibroyden = 0;
  artol = 1.E-5;
  gtol = 1.5E-8;
  xtol = 1.E-10;
  xrtol = 1.5E-8;
  frtol = 1.5E-8;
  ftol = 1.5E-8;
  minlam = -1;

/* input not initialized */
/*
  dtd =
*/

/**************************************
* call geodesiclm and write results
**************************************/
  geodesiclm_(func_, jacobian_, Avv_, 
              x, fvec, fjac[0], &n, &m, 
              callback_, &info, 
              &analytic_jac, &analytic_Avv, 
              &center_diff, &h1, &h2,
              dtd[0], &damp_mode, 
              &niters, &nfev, &njev, &naev, 
              &maxiter, &maxfev, &maxjev, &maxaev, &maxlam, &minlam, 
              &artol, &Cgoal, &gtol, &xtol, &xrtol, &ftol, &frtol, 
              &converged, 
              &print_level, &print_unit, 
              &imethod, &iaccel, &ibold, &ibroyden, 
              &initialfactor, &factoraccept, &factorreject, &avmax);

printf("The optimized parameters are: %5.3f %5.3f %5.3f %5.3f\n",x[0],x[1],x[2],x[3]);


  /* free memory */
  free(x);
  free(fvec);
  free(dtd[0]);
  free(dtd);
  free(fjac[0]);
  free(fjac);


  return 0;
}

