
/******************************************************************************
*
*  PoLSA - Powell Least Square Algorithm
* 
*  powell.c: Contains Powell optimization algorithm and some small subroutines
*
******************************************************************************/

/****************************************************************
* $Revision: 1.16 $
* $Date: 2003/07/29 08:42:53 $
*****************************************************************/

/******************************************************************************
*
* This utility will optimize a parameter vector xi[N] to minimize the sum of 
* squares U=sum_{i=1..M}(f_i(xi)-F_i)^2, where F_i is an array passed to the
* utility, and f_i(xi) is a function of the vector xi.
*
******************************************************************************/

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
   int i,j,k,m=0,n=0;                    /* Simple counting variables */
    real *force_xi;                /* calculated force, alt */
    real **d;                      /* Direction vectors */
    real **gamma;                  /* Matrix of derivatives */
    real **lineqsys;               /* Lin.Eq.Sys. Matrix */
    real **les_inverse;		   /* LU decomp. of the lineqsys */
    real *delta;		   /* Vector pointing into correct dir'n */
    real *delta_norm;		   /* Normalized vector delta */
    real *fxi1,*fxi2;		   /* two latest force vectors */
    int *perm_indx;		   /* Keeps track of LU pivoting */
    real *p,*q;                    /* Vectors needed in Powell's algorithm */
    real perm_sig;                 /* Signature of permutation in LU decomp */
    real F,F2,F3,df,xi1,xi2;       /* Fn values, changes, steps ... */
    real temp,temp2;               /* as the name indicates: temporary vars */
    static char errmsg[256];	   /* Error message */
    FILE *ff;			   /* Exit flagfile */

    d=dmatrix(0,ndim-1,0,ndim-1);
    gamma=dmatrix(0,mdim-1,0,ndim-1);
    lineqsys=dmatrix(0,ndim-1,0,ndim-1);
    les_inverse=dmatrix(0,ndim-1,0,ndim-1);
    perm_indx=ivector(0,ndim-1);
    delta=dvector(0,ndimtot-1); /* ==0*/
    delta_norm=dvector(0,ndimtot-1); /*==0*/
    fxi1=dvector(0,mdim-1);
    fxi2=dvector(0,mdim-1);
    force_xi=dvector(0,mdim-1);
    p=dvector(0,ndim-1);
    q=dvector(0,ndim-1);
    

    /* Clear vector delta */
    for(k=0;k<ndimtot;k++) delta[k]=0;

    /* calculate the first vector */
    F=(*calc_forces)(xi,fxi1);
    
    printf("%d %f %f %f %f %f %f %d \n",
	       m,F,xi[0],xi[1],xi[2],xi[3],xi[4],fcalls);
    fflush(stdout);
    if (F<NOTHING) return;	/* If F is less than nothing, */
				/* what is there to do?*/
    (void) copy_vector(fxi1,force_xi,mdim);
    do { /*outer loop, includes recalculating gamma*/
	m=0;
	/* Init gamma */
    	if (i=gamma_init(gamma, d, xi, fxi1, ndim, mdim)) {
	  write_pot_table( &pair_pot, tempfile ); /*emergency writeout*/    
	  sprintf(errmsg, "F does not depend on xi[%d], fit impossible!\n",
		    idx[i-1]);
	    error(errmsg);
	}
    	(void) lineqsys_init(gamma,lineqsys,fxi1,p,ndim,mdim); /*init LES */
    	F3=F;
    	do {
            /* LU decomposition of LES matrix */
	    /* 1) Store original matrix (needed later) */
	    (void) copy_matrix(lineqsys,les_inverse,ndim,ndim);
            /* 2) LU decompose */
	    if (ludcmp_r(les_inverse,ndim,perm_indx,&perm_sig)) break;
            /* 3) Store original vector (needed later) */
	    (void) copy_vector(p,q,ndim);
	    /* 4) Backsubstitute vector in lu decomp. matrix */
            (void) lubksb_r(les_inverse,ndim,perm_indx,q); 
	    /* 5) Improve that solution by iterative improvement */
            (void) mprove_r(lineqsys,les_inverse,ndim,perm_indx,p,q); 
	    
	    /* get delta by multiplying q with the direction vectors */
            (void) matdotvec(d,q,delta,ndim,ndim);
            /* store delta */
            (void) copy_vector(delta, delta_norm, ndimtot);
	    /* normalize it */
/*	    (void) normalize_vector(delta_norm,ndimtot); */
            
	    F2=F;   /*shift F*/
	    /* minimize F(xi) along vector delta, return new F */
	    F=linmin_r(xi,delta,F,ndim,mdim,&xi1,&xi2,fxi1,fxi2);
    
            j=0;temp2=0.;
            for (i=0;i<ndim;i++) 
                if ((temp=fabs(p[i]*q[i]))>temp2) {j=i;temp2=temp;};
	    
	    /* Check for degeneracy in directions */
/* 	    for (i=0;i<ndim;i++) { */
/* 		if (i!=j) { */
/* 		    temp=0; */
/* 	            for (k=0;k<ndim;k++) temp+=d[k][i]*d[k][j]; */
/* 		    if (1-temp<=.0001) break;  /\*Emergency exit*\/ */
/* 	        } */
/* 	    } */
	    /*update gamma, but if fn returns 1, matrix will be sigular, 
	    break inner loop and restart with new matrix*/
            if (gamma_update(gamma,xi1,xi2,fxi1,fxi2,delta_norm,j,mdim,ndimtot,F)){
		sprintf(errmsg,"Matrix gamma singular after step %d,\nrestarting inner loop\n",m);
		warning(errmsg);
		break;
	    }
            /*set new d[j]*/
            for (i=0;i<ndim;i++) d[i][j]=delta_norm[idx[i]]; 

            (void) lineqsys_update(gamma,lineqsys,fxi1,p,j,ndim,mdim);
	    m++;           /*increment loop counter*/
	/* loop at least 7*ndim times, but at most INNERLOOPS or until no
	further improvement */ 
	    df=F2-F;  
        } while ((m<7*ndim+1 || (m<=INNERLOOPS && df>PRECISION) ) &&
			df<TOOBIG); 
	    /* Print the steps in current loop, F, a few values of xi, and
	       total number of fn calls */
        printf("%d %f %f %f %f %f %f %d \n",
	  m,F,xi[0],xi[1],xi[2],xi[3],xi[4],fcalls);
	fflush(stdout);
        n++;			/* increment outer loop counter */
        /* Output of potential at intermediate steps */
	/* printf("# %d\n",n);
	for (k=0;k<ndim;k++)
	    printf("%f %f\n",pair_pot.begin[0]+k*pair_pot.step[0],xi[k]);
	    printf("\n\n");*/
        /* End fit if break flagfile exists */
	ff = fopen(flagfile,"r");
	if (NULL != ff) {
	    printf("Fit terminated prematurely in presence of break flagfile %s!\n",
		   flagfile);
	    break;
	}
	if (tempfile != "\0" ) write_pot_table( &pair_pot, tempfile );
    /*End fit if whole series didn't improve F */
    } while (F3-F>PRECISION/10.);
	/* Free memory */
    free_dmatrix(d,0,ndim-1,0,ndim-1);
    free_dmatrix(gamma,0,mdim-1,0,ndim-1);
    free_dmatrix(lineqsys,0,ndim-1,0,ndim-1);
    free_dmatrix(les_inverse,0,ndim-1,0,ndim-1);
    free_ivector(perm_indx,0,ndim-1);
    free_dvector(delta,0,ndim-1);
    free_dvector(delta_norm,0,ndim-1);
    free_dvector(force_xi,0,mdim-1);
    free_dvector(fxi1,0,mdim-1);
    free_dvector(fxi2,0,mdim-1);
    free_dvector(p,0,ndim-1);
    free_dvector(q,0,ndim-1);

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

int gamma_init(real **gamma, real **d, real *xi, real *force_xi, int n,
         	int m){
    real *force;
    int i,j;			/* Auxiliary vars: Counters */
    real sum,temp;			/* Auxiliary var: Sum */
/*   Set direction vectors to coordinate directions d_ij=KroneckerDelta_ij */
    /*Initialize direction vectors*/
    for (i=0;i<n;i++) {
	    for (j=0;j<n;j++) if (i==j) d[i][j]=1.; else d[i][j]=0.;
    }
/* Initialize gamma by calculating numerical derivatives    */
    force=dvector(0,m-1);
    for (i=0;i<n;i++) {           /*initialize gamma*/
	xi[idx[i]]+=EPS;                /*increase xi[idx[i]]...*/
	sum = 0.;
	(void) (*calc_forces)(xi,force);
       	for (j=0;j<m;j++) {
	  temp =(force[j]-force_xi[j])/EPS;
	  gamma[j][i]=temp;
	  sum += SQR(temp);
	}
	xi[idx[i]]-=EPS;                /*...and decrease xi[idx[i]] again*/
/* scale gamma so that sum_j(gamma^2)=1                      */
	temp = sqrt(sum);
	if (temp>NOTHING) { 
	    for (j=0;j<m;j++) gamma[j][i] /= temp; /*normalize gamma*/
	    d[i][i]/=temp;	/* rescale d */
	}
	else
	    return i+1;		/* singular matrix, abort */

    }
    free_dvector(force,0,m-1);

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
	    temp=gamma[i][j]=((fa[i]-fb[i])/(a-b));
	    mu+=temp*fa[i];
    }
    mu/=fmin;
    for (i=0;i<m;i++) {
      temp=gamma[i][j]-=mu*fa[i];
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

    /* calculating vector p (lineqsys . q == P in LinEqSys) */

    for (i=0;i<n;i++) {
	p[i]=0.;
	for (j=0;j<m;j++) {
	    p[i]-=gamma[j][i]*deltaforce[j];
	}
    }
    /* calculating the linear equation system matrix gamma^t.gamma */
    for (i=0;i<n;i++) {
	for (k=0;k<n;k++) {
	    lineqsys[i][k]=0;
	    for (j=0;j<m;j++) {
		lineqsys[i][k]+=gamma[j][i]*gamma[j][k];
	    }
	}
    }
    return;
}

/*******************************************************************
*
* lineqsys_update: Update LinEqSys matrix line and column i, vector
*            p component i.
*
*******************************************************************/

void lineqsys_update(real **gamma, real **lineqsys, real *force_xi,
		     real *p, int i, int n, int m){
#ifdef _OPENMP
#pragma omp parallel
#endif 
  {
    int j,k;
    real temp;
/* #ifdef _OPENMP3 */
/* #pragma omp sections */
/* #endif */
/*  { */
/* #ifdef _OPENMP3 */
/* #pragma omp section */
/* #endif  */
/*     for (k=0;k<n;k++){ */
/*     	p[k]=0.; */
/*     	for (j=0;j<m;j++) p[k]-=gamma[j][k]*force_xi[j]; */
/*     } */
/* #ifdef _OPENMP3 */
/* #pragma omp section */
/* #endif  */
/*     /\* divide calculation of new line i and column i in 3 parts*\/  */
/*     for (k=0;k<i;k++) {           /\*Part 1: All elements ki and ik with k<i*\/ */
/*         lineqsys[i][k]=0.; */
/* 	lineqsys[k][i]=0.; */
/* 	for (j=0;j<m;j++) { */
/* 	    temp=gamma[j][i]*gamma[j][k]; */
/* 	    lineqsys[i][k]+=temp; */
/* 	    lineqsys[k][i]+=temp; */
/*         } */
/*     } */
/*     /\* divide calculation of new line i and column i in 3 parts*\/  */
#ifdef _OPENMP
#pragma omp for
#endif
    for (k=0;k<n;k++) {           /*Part 1: All elements ki and ik with k<i*/
      p[k]=0.;  
      lineqsys[i][k]=0.;
      lineqsys[k][i]=0.;
      for (j=0;j<m;j++) {
	p[k]-=gamma[j][k]*force_xi[j];
	temp=gamma[j][i]*gamma[j][k];
	lineqsys[i][k]+=temp;
	lineqsys[k][i]+=temp;
      }
    }
/* #ifdef _OPENMP3 */
/* #pragma omp section */
/* #endif  */
/*  { */
/*     lineqsys[i][i]=0.;             /\*Part 2: Element ii *\/ */
/*     for (j=0;j<m;j++) { */
/* 	lineqsys[i][i]+=SQR(gamma[j][i]); */
/*     } */
/*  } */
/* #ifdef _OPENMP3 */
/* #pragma omp section */
/* #endif  */
/*     for (k=i+1;k<n;k++){	   /\* Part 3: All elements ki and ik with k>i*\/ */
/*         lineqsys[i][k]=0; */
/* 	lineqsys[k][i]=0; */
/* 	for (j=0;j<m;j++) { */
/* 	    temp=gamma[j][i]*gamma[j][k]; */
/* 	    lineqsys[i][k]+=temp; */
/* 	    lineqsys[k][i]+=temp; */
/*         } */
/*     } */
/*  } */
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

