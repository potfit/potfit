
/******************************************************************************
*
*  PoLSA - Powell Least Square Algorithm
* 
*  powell.c: Contains Powell optimization algorithm and some small subroutines
*
******************************************************************************/

/****************************************************************
* $Revision: 1.18 $
* $Date: 2004/02/25 16:39:09 $
*****************************************************************/

/******************************************************************************
*
* This utility will optimize a parameter vector xi[N] to minimize the sum of 
* squares U=sum_{i=1..M}(f_i(xi)-F_i)^2, where F_i is an array passed to the
* utility, and f_i(xi) is a function of the vector xi.
*
******************************************************************************/

#include <mkl_lapack.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "potfit.h"
#include "nrutil_r.h"
#include "powell_lsq.h"
#define EPS .0001
#define PRECISION 1.E-7
/* Well, almost nothing */
#define NOTHING 1.E-12
#define INNERLOOPS 801
#define TOOBIG 10000


void powell_lsq(real *xi)
{
  char uplo[1]="U"; 		/* char used in dsysvx */
  char fact[1]="N"; 		/* char used in dsysvx */
  int i,j,k,m=0,n=0;                    /* Simple counting variables */
  real *force_xi;                /* calculated force, alt */
  real **d;                      /* Direction vectors */
  real **gamma;                  /* Matrix of derivatives */
  real **lineqsys;               /* Lin.Eq.Sys. Matrix */
  real **les_inverse;		   /* LU decomp. of the lineqsys */
  real *delta;		   /* Vector pointing into correct dir'n */
  real *delta_norm;		   /* Normalized vector delta */
  real *fxi1,*fxi2;		   /* two latest force vectors */
  real *work;			/* work array to be used by dsysvx */
  int *iwork;
  int *perm_indx;		   /* Keeps track of LU pivoting */
  int worksize;			/* Size of work array (dsysvx) */
  real cond;			/* Condition number dsysvx */
  real *p,*q;                    /* Vectors needed in Powell's algorithm */
  real perm_sig;                 /* Signature of permutation in LU decomp */
  real F,F2,F3,df,xi1,xi2;       /* Fn values, changes, steps ... */
  real temp,temp2;               /* as the name indicates: temporary vars */
  real ferror,berror;		/* forward/backward error estimates */
  static char errmsg[256];	   /* Error message */
  FILE *ff;			   /* Exit flagfile */
  d=dmatrix(0,ndim-1,0,ndim-1);
  gamma=dmatrix(0,mdim-1,0,ndim-1);
  lineqsys=dmatrix(0,ndim-1,0,ndim-1);
  les_inverse=dmatrix(0,ndim-1,0,ndim-1);
  perm_indx=ivector(0,ndim-1);
  delta_norm=dvector(0,ndimtot-1); /*==0*/
  force_xi=dvector(0,mdim-1);
  p=dvector(0,ndim-1);
  q=dvector(0,ndim-1);
  delta=dvector(0,ndimtot-1); /* ==0*/
  fxi1=dvector(0,mdim-1);
  fxi2=dvector(0,mdim-1);
  worksize=64*ndim;
  work=(real *) malloc(worksize*sizeof(real));
  iwork=(int *) malloc(ndim*sizeof(int));

  /* calculate the first force */
  F=(*calc_forces)(xi,fxi1,0);

  printf("%d %f %f %f %f %f %f %d \n",
	 m,F,xi[0],xi[1],xi[2],xi[3],xi[4],fcalls);
  fflush(stdout);

  if (F<NOTHING) return;	/* If F is less than nothing, */
				/* what is there to do?*/

  (void) copy_vector(fxi1,force_xi,mdim);

  do { /*outer loop, includes recalculating gamma*/
    m=0;

    /* Init gamma */
    if (i=gamma_init(gamma, d, xi, fxi1)) {
      write_pot_table( &pair_pot, tempfile ); /*emergency writeout*/    
      sprintf(errmsg, "F does not depend on xi[%d], fit impossible!\n",
	      idx[i-1]);
      warning(errmsg);
      break;
    }

    (void) lineqsys_init(gamma,lineqsys,fxi1,p,ndim,mdim); /*init LES */
    F3=F;

    do { /*inner loop - only calculate changed rows/lines in gamma */
      /* (a) solve linear equation */

      /* All in one driver routine */
      j=1;			/* 1 rhs */

      /* Linear Equation Solution (lapack) */
      dsysvx(fact,uplo,&ndim,&j,&lineqsys[0][0],&ndim,&les_inverse[0][0],\
	     &ndim,perm_indx,p,&ndim,q,&ndim,&cond,&ferror,&berror,work,\
	     &worksize,iwork,&i);
      
#ifdef DEBUG
       printf("q0: %d %f %f %f %f %f %f %f %f\n",i, q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7]);
#endif

      /* (b) get delta by multiplying q with the direction vectors */
      matdotvec(d,q,delta,ndim,ndim);
      /*     and store delta */
      copy_vector(delta, delta_norm, ndimtot);

      F2=F;   /*shift F*/

      /* (c) minimize F(xi) along vector delta, return new F */
      F=linmin_r(xi,delta,F,ndim,mdim,&xi1,&xi2,fxi1,fxi2);

#ifdef DEBUG
      printf("%f %6g %f %f %d \n",F,cond,ferror,berror,i);
#endif

      /* (d) if error estimate is too high after minimization
	 in 5 directions: restart outer loop */
      if (ferror+berror>1. && m>5) break;

      /* (e) find optimal direction to replace */
      j=0;temp2=0.;
      for (i=0;i<ndim;i++) 
	if ((temp=fabs(p[i]*q[i]))>temp2) {j=i;temp2=temp;};

      /* (f) update gamma, but if fn returns 1, matrix will be sigular, 
	 break inner loop and restart with new matrix*/
      if (gamma_update(gamma,xi1,xi2,fxi1,fxi2,delta_norm,j,mdim,ndimtot,F)) {
	sprintf(errmsg,
		"Matrix gamma singular after step %d, restarting inner loop",m);
	warning(errmsg);
	break;
      } 

      /* (g) set new direction vector */
      for (i=0;i<ndim;i++) d[i][j]=delta_norm[idx[i]]; 
      
      /* (h) update linear equation system */
      lineqsys_update(gamma,lineqsys,fxi1,p,j,ndim,mdim);

      m++;           /*increment loop counter*/
      df=F2-F;

      /* loop at least 7*ndim times, but at most INNERLOOPS or until no
	 further improvement */ 
    } while ((m<7*ndim+1 || (m<=INNERLOOPS && df>PRECISION) ) &&
	     df<TOOBIG); 	/* inner loop */

    n++;			/* increment outer loop counter */

    /* Print the steps in current loop, F, a few values of xi, and
       total number of fn calls */
    printf("%d %f %f %f %f %f %f %d \n",
	   m,F,xi[0],xi[1],xi[2],xi[3],xi[4],fcalls);
    fflush(stdout);
    
    /* End fit if break flagfile exists */
    ff = fopen(flagfile,"r");
    if (NULL != ff) {
      printf("Fit terminated prematurely in presence of break flagfile %s!\n",
	     flagfile);
      break;
    }

    /* write temp file  */
    if (tempfile != "\0" ) write_pot_table( &pair_pot, tempfile );

    /*End fit if whole series didn't improve F */
  } while (F3-F>PRECISION/10.); /* outer loop */

  /* Free memory */
  free_dvector(delta,0,ndim-1);
  free_dvector(fxi1,0,mdim-1);
  free_dvector(fxi2,0,mdim-1);
  free_dmatrix(d,0,ndim-1,0,ndim-1);
  free_dmatrix(gamma,0,mdim-1,0,ndim-1);
  free_dmatrix(lineqsys,0,ndim-1,0,ndim-1);
  free_dmatrix(les_inverse,0,ndim-1,0,ndim-1);
  free_ivector(perm_indx,0,ndim-1);
  free_dvector(delta_norm,0,ndim-1);
  free_dvector(force_xi,0,mdim-1);
  free_dvector(p,0,ndim-1);
  free_dvector(q,0,ndim-1);
  free(work);
  free(iwork);
  return ;
}


/******************************************************************
*
* gamma_init: (Re-)Initialize gamma[j][i] (Gradient Matrix) after
*            (Re-)Start or whenever necessary by calculating numerical
*            gradients in coordinate directions. Includes re-setting the
*            direction vectors to coordinate directions.
*
******************************************************************/

int gamma_init(real **gamma, real **d, real *xi, real *force_xi)
{
  real *force;
  int i,j;			/* Auxiliary vars: Counters */
  real sum,temp;			/* Auxiliary var: Sum */
/*   Set direction vectors to coordinate directions d_ij=KroneckerDelta_ij */
  /*Initialize direction vectors*/
  for (i=0;i<ndim;i++) {
    for (j=0;j<ndim;j++) d[i][j]=(i==j)?1.:0.;
  }
/* Initialize gamma by calculating numerical derivatives    */
  force=dvector(0,mdim-1);
  for (i=0;i<ndim;i++) {           /*initialize gamma*/
    xi[idx[i]]+=EPS;                /*increase xi[idx[i]]...*/
    sum = 0.;
    (void) (*calc_forces)(xi,force,0);
    for (j=0;j<mdim;j++) {
      temp =(force[j]-force_xi[j])/EPS;
      gamma[j][i]=temp;
      sum += SQR(temp);
    }
    temp = sqrt(sum);
    xi[idx[i]]-=EPS;                /*...and decrease xi[idx[i]] again*/
/* scale gamma so that sum_j(gamma^2)=1                      */
    if (temp>NOTHING) { 
//      printf("%f\n",temp);
      for (j=0;j<mdim;j++) gamma[j][i] /= temp; /*normalize gamma*/
      d[i][i]/=temp;	/* rescale d */
    }
    else
      return i+1;		/* singular matrix, abort */
  }
  free_dvector(force,0,mdim-1);
/*     for (i=0;i<n;i++) { */
/* 	sum =0.; */
/* 	for (j=0;j<m;j++) sum += SQR(gamma[j][i]); */
/* 	temp=sqrt(sum); */
/* 	if (temp>NOTHING)  */
/* 	    for (j=0;j<m;j++) gamma[j][i] /= temp; /\*normalize gamma*\/ */
/* 	else */
/* 	    return i+1;		/\* singular matrix, abort *\/ */
/*     } */
  return 0;
}

/*******************************************************************
*
* gamma_update: Update column j of gamma ( to newly calculated 
*           numerical derivatives (calculated from fa, fb 
*           at a,b); normalize new vector.
*
*******************************************************************/

int  gamma_update(real **gamma, real a, real b, real *fa, real *fb, real *delta,
		int j, int m, int n, real fmin) {
    int i;
    real temp;
    real sum=0.;
    real mu=0.;
    for (i=0;i<m;i++) {
	    temp=((fa[i]-fb[i])/(a-b));
	    gamma[i][j]=temp;
	    mu+=temp*fa[i];
    }
    mu/=fmin;
    for (i=0;i<m;i++) {
      temp=gamma[i][j]-mu*fa[i];
      gamma[i][j]=temp;
      sum+=temp*temp;
    }
    temp=sqrt(sum);               /* normalization factor */
    if (temp>NOTHING) {
    	for (i=0;i<m;i++) gamma[i][j]/=temp; 
	for (i=0;i<n;i++) delta[i]/=temp;
    }
    else
	return 1;  /*Matrix will be singular: Restart!*/
    return 0;
}
		
/*******************************************************************
*
* lineqsys_init: Initialize LinEqSys matrix, vector p in 
*              lineqsys . q == p
*
*******************************************************************/

void lineqsys_init(real **gamma, real **lineqsys, real *deltaforce, 
	real *p, int n, int m){
    int i,j,k;			/* Auxiliary vars: Counters */
    real temp;
    /* calculating vector p (lineqsys . q == P in LinEqSys) */

    for (i=0;i<n;i++) {
	p[i]=0.;
	for (j=0;j<m;j++) {
	    p[i]-=gamma[j][i]*deltaforce[j];
	}
    }
    /* calculating the linear equation system matrix gamma^t.gamma */
   for (i=0;i<n;i++) {
      lineqsys[i][i]=0;
      for (j=0;j<m;j++)
	lineqsys[i][i]+=SQR(gamma[j][i]);
      for (k=i+1;k<n;k++) {
//      for (k=0;k<n;k++) {
	lineqsys[i][k]=0.;
//	lineqsys[k][i]=0.;
	for (j=0;j<m;j++) {
	  lineqsys[i][k]+=gamma[j][i]*gamma[j][k];
	}
	lineqsys[k][i]=lineqsys[i][k];
      }
    }
    return;
}

/*******************************************************************
*
* lineqsys_update: Update LinEqSys matrix row and column i, vector
*            p.
*
*******************************************************************/

void lineqsys_update(real **gamma, real **lineqsys, real *force_xi,
		     real *p, int i, int n, int m){
  int k;
#ifdef _OPENMP
#pragma omp parallel
#endif 
  {
/*     int j,k; */
/*     real temp; */
#ifdef _OPENMP
#pragma omp for
#endif
    for (k=0;k<n;k++) { 
      int j;
      real temp;
      p[k]=0.;  
      lineqsys[i][k]=0.;
      for (j=0;j<m;j++) {
	p[k]-=gamma[j][k]*force_xi[j];
	lineqsys[i][k]+=gamma[j][i]*gamma[j][k];
      }
      lineqsys[k][i]=lineqsys[i][k];
    }
  }
  
  return;
}
	


/*******************************************************************
*
*  copy_matrix: Copies data from Matrix a into matrix b
*              (matrix dimension n x m)
*
*******************************************************************/

/*>>>>>>>>>>>>>>>  INLINING? <<<<<<<<<<<<<<*/

void copy_matrix(real **a, real **b, int n, int m) {
    int i,j;
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
	    b[j][i]=a[j][i];
        }
    }
    return;
}

/******************************************************************
*
* copy_vector: Copies data from vector a into vector b (both dim n)
*
******************************************************************/

/*>>>>>>>>>>>>>>>  INLINING? <<<<<<<<<<<<<<*/

void copy_vector(real *a, real *b, int n) {
    int i;
    for (i=0; i<n; i++) b[i]=a[i];
    return;
}

/******************************************************************
*
* matdotvec: Calculates the product of matrix a (n x m) with column 
* vector x (dim m), Result y (dim n). (A . x = m)
*
******************************************************************/

void matdotvec(real **a, real *x, real *y, int n, int m){
    int i, j;
    for (i=0; i<n; i++) {
	y[idx[i]]=0.;
	for (j=0; j<m; j++) y[idx[i]]+=a[i][j]*x[j];
    }
    return;
}

/*******************************************************************
*
* normalize_vector: Normalizes vector to |vec|^2=1, returns norm of old vector
*
*******************************************************************/

real normalize_vector(real *v, int n){
    int j;
    real temp,sum=0.0;
    for (j=0;j<n;j++) sum+=SQR(v[j]);
    temp=sqrt(sum);
    for (j=0;j<n;j++) v[j]/=temp;
    return temp;
}

