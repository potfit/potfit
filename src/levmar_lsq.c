/******************************************************************************
* functions to use geodesic Levenberg-Marquardt algorithm
* 
* Author: Mingjian Wen (wenxx151@umn.edu), University of Minnesota
*******************************************************************************/

#include "potfit.h"
#include "potential_output.h"
#include "force.h"
#include "optimize.h"

/* exchange access order for fortran code */
#define f_dtd(i,j) dtd[j][i]

/* function decleration */ 
typedef void(*fptr_func_)(int*, int*, double*, double*);
typedef void(*fptr_jacobian_)(int*, int*, double*, double*); 
typedef void(*fptr_Avv_)(int*, int*, double*, double*, double*);
typedef void(*fptr_callback_)(int*, int*, double*, double*, double*, double*,
                              double*, double*, double*, double*,
                              double*, int*, int*); 

void geodesiclm_(fptr_func_, fptr_jacobian_, fptr_Avv_, 
                double* x, double* fvec, double*, int* n, int* m, 
                fptr_callback_, int* info, 
                int* analytic_jac, int* analytic_Avv, 
                int* center_diff, double* h1, double* h2,
                double*, int* damp_mode, 
                int* niters, int* nfev, int* njev, int* naev, 
                int* maxiter, int* maxfev, int* maxjev, int* maxaev, double* maxlam, 
                double* minlam, double* artol, double* Cgoal, double* gtol, 
                double* xtol, double* xrtol, double* ftol, double* frtol, 
                int* converged, 
                int* print_level, int* print_unit, 
                int* imethod, int* iaccel, int* ibold, int* ibroyden, 
                double* initialfactor, double* factoraccept, 
                double* factorreject, double* avmax);

/******************************************************************************
* local fucntions
******************************************************************************/
/* jacobian */
void jacobian_(int* m, int* n, double* x, double* fjac){}
/* Avv */
void Avv_(int* m, int* n, double* x, double* v, double* acc){}

/* callback */
void callback_(int* m, int* n, double* x,double* v,double* a,double* fvec,
               double* fjac, double* acc, double* lam, double* dtd,
               double* fvec_new, int* accepted, int* info)
{
  /* write best fit so far */
  write_pot_table_potfit(g_files.tempfile);
}

/******************************************************************************
* Force parameters relation function 

* idx: index of optimizable parameters in potential table
******************************************************************************/
void force_param_relation_(int* m, int* n, double *p, double *force)
{
  int i;

  /* get access to global variable */
  double* xi = g_pot.opt_pot.table;
	int* idx = g_pot.opt_pot.idx;

  /* update the optimizable values */
  for (i=0; i<*n; ++i)
    xi[idx[i]] = p[i];

  /* init forces */
  for(i=0; i<*m; ++i)
    force[i] = 0;

  /* update force */
  calc_forces(xi, force, 0);

}

/******************************************************************************
* 
* potfit global variables:
* ndim: number of optimizable params
* mdim: general force dimension 
* idx: index of optimizable parameters in potential table  
******************************************************************************/

int run_levmar_lsq(double* xi)
{
  /* local variables */
  int i,j;
  /* input */
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

	/* get access to glboal variable */
	int* idx = g_pot.opt_pot.idx;
  int n = g_calc.ndim;
  int m = g_calc.mdim;

 /* allocate memory */
  x = (double*)malloc(n*sizeof(double));
  /* we need to allocate memory for fvec, even though it is output*/
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
  for (i=0; i<n; ++i)
    x[i] = xi[idx[i]];

  /* initialize dtd in order to use damp mode 1 */
  for (i=0; i<n; ++i)
    for(j=0; j<n; ++j)
      if(j==i) {
        f_dtd(i,j) = 1.0;
     } else
        f_dtd(i,j) = 0.0;
  
  /* TO DO, the parameters needs to be tuned futher; also needs to enable setting 
    them in the parameter file. */
  /* set geodesicLM parameters */
  analytic_jac=0;
  analytic_Avv=0;
  center_diff=0;  
  
  h1 = 1.0E-5;
  h2 = 0.1;
  factoraccept = 5;
  factorreject = 2;
  maxlam = 1.E15;
  minlam = -1;
  
  imethod = 2;
  initialfactor = 1;
  damp_mode = 0;
  maxiter = 10000;
  maxfev = 0;
  Cgoal = 0.5E-7;
  maxjev = 0;
  iaccel = 1;
  avmax = 0.75;
  maxaev = 0;
  print_level =2;
  print_unit = 6;
  ibold = 2;
  ibroyden = 0;
  artol = 1.E-8;
  gtol = 1.5E-8;
  xtol = 1.E-6;
  xrtol = 1.5E-6;
  frtol = 1.5E-6;
  ftol = 1.5E-6;

/* input not initialized */
/*
  dtd =
*/

/* call geodesiclm */
  geodesiclm_(force_param_relation_, jacobian_, Avv_, 
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


#if !defined(KIM)
#if defined(APOT)
  update_apot_table(xi);
#endif
#endif

  /* put the optimized params back to xi */
  for (i=0; i<n; ++i)
    xi[idx[i]] = x[i];

  /* free memory */
  free(x);
  free(fvec);
  free(dtd[0]);
  free(dtd);
  free(fjac[0]);
  free(fjac);

  return 0;
}
