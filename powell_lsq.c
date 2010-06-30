/****************************************************************
 *
 * powell_lsq.c: Contains Powell optimization algorithm
 *	and some small subroutines
 *
 ****************************************************************
 *
 * Copyright 2002-2010 Peter Brommer, Daniel Schopf
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://www.itap.physik.uni-stuttgart.de/
 *
 ****************************************************************
 *
 *   This file is part of potfit.
 *
 *   potfit is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   potfit is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

/******************************************************************************
*
* This utility will optimize a parameter vector xi[N] to minimize the sum of
* squares U=sum_{i=1..M}(f_i(xi)-F_i)^2, where F_i is an array passed to the
* utility, and f_i(xi) is a function of the vector xi.
*
******************************************************************************/

#ifdef ACML
#include <acml.h>
#else /* ACML */
#include <mkl_lapack.h>
#endif /* ACML */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "potfit.h"
#include "utils.h"
#include "bracket.h"
#include "powell_lsq.h"
#define EPS .0001
#define PRECISION 1.E-7
/* Well, almost nothing */
#define NOTHING 1.E-12
#define INNERLOOPS 801
#define TOOBIG 10000


void powell_lsq(real *xi)
{
  char  uplo[1] = "U";		/* char used in dsysvx */
  char  fact[1] = "N";		/* char used in dsysvx */
  int   i, j, m = 0, n = 0;	/* Simple counting variables */
  real *force_xi;		/* calculated force, alt */
  real **d;			/* Direction vectors */
  real **gamma;			/* Matrix of derivatives */
  real **lineqsys;		/* Lin.Eq.Sys. Matrix */
  real **les_inverse;		/* LU decomp. of the lineqsys */
  real *delta;			/* Vector pointing into correct dir'n */
  real *delta_norm;		/* Normalized vector delta */
  real *fxi1, *fxi2;		/* two latest force vectors */
#ifndef ACML			/* work arrays not needed for ACML */
  real *work;			/* work array to be used by dsysvx */
  int  *iwork;
  int   worksize;		/* Size of work array (dsysvx) */
#endif /* ACML */
  int  *perm_indx;		/* Keeps track of LU pivoting */
  int   breakflag;		/* Breakflag */
  real  cond;			/* Condition number dsysvx */
  real *p, *q;			/* Vectors needed in Powell's algorithm */
/*  real  perm_sig;		|+ Signature of permutation in LU decomp +|*/
  real  F, F2, F3 = 0, df, xi1, xi2;	/* Fn values, changes, steps ... */
  real  temp, temp2;		/* as the name indicates: temporary vars */
#ifdef APOT
  int   itemp, itemp2;		/* the same for integer */
#endif
  real  ferror, berror;		/* forward/backward error estimates */
  static char errmsg[256];	/* Error message */
  FILE *ff;			/* Exit flagfile */

  d = mat_real(ndim, ndim);
  gamma = mat_real(mdim, ndim);
  lineqsys = mat_real(ndim, ndim);
  les_inverse = mat_real(ndim, ndim);
  perm_indx = vect_int(ndim);
  delta_norm = vect_real(ndimtot);
				 /*==0*/
  force_xi = vect_real(mdim);
  p = vect_real(ndim);
  q = vect_real(ndim);
  delta = vect_real(ndimtot);	/* ==0 */
  fxi1 = vect_real(mdim);
  fxi2 = vect_real(mdim);
#ifndef ACML			/* work arrays not needed */
  worksize = 64 * ndim;
  work = (real *)malloc(worksize * sizeof(real));
  iwork = (int *)malloc(ndim * sizeof(int));
#endif /* ACML */

  /* clear delta */
  for (i = 0; i < ndimtot; i++)
    delta[i] = 0.;

  /* calculate the first force */
  F = (*calc_forces) (xi, fxi1, 0);
#ifndef APOT
  printf("%d %f %f %f %f %f %f %d\n",
	 m, F, xi[0], xi[1], xi[2], xi[3], xi[4], fcalls);
  fflush(stdout);
#endif

  if (F < NOTHING) {
    printf("Error already too small to optimize, aborting ...\n");
    return;			/* If F is less than nothing, */
    /* what is there to do? */
  }
  (void)copy_vector(fxi1, force_xi, mdim);
#ifdef APOT
  printf("loops\t\terror_sum\tforce calculations\n");
  printf("%5d\t%17.6f\t%6d\n", m, F, fcalls);
#endif
  do {				/*outer loop, includes recalculating gamma */
    m = 0;
    /* Init gamma */
    if (i = gamma_init(gamma, d, xi, fxi1)) {
#ifdef EAM
#ifndef NORESCALE
      /* perhaps rescaling helps? - Last resort... */
      sprintf(errmsg, "F does not depend on xi[%d], trying to rescale!\n",
	      idx[i - 1]);
      warning(errmsg);
      rescale(&opt_pot, 1., 1);
      /* wake other threads and sync potentials */
      F = calc_forces(xi, fxi1, 2);
/*      printf("F=%f\n", F);*/
      i = gamma_init(gamma, d, xi, fxi1);
#endif /* NORESCALE */
#endif /* EAM */
      /* try again */
      if (i) {
	/* ok, now this is serious, better exit cleanly */
#ifndef APOT
	write_pot_table(&opt_pot, tempfile);	/*emergency writeout */
	sprintf(errmsg, "F does not depend on xi[%d], fit impossible!\n",
		idx[i - 1]);
#else
	for (n = 0; n < ndim; n++)
	  apot_table.values[apot_table.idxpot[n]][apot_table.idxparam[n]] =
	    xi[idx[n]];
	write_pot_table(&apot_table, tempfile);
	itemp = apot_table.idxpot[i - 1];
	itemp2 = apot_table.idxparam[i - 1];
#ifdef DIPOLE
	sprintf(errmsg,
		"F does not depend on the %d. parameter (%s) of the %d. potential.\nFit impossible!\n",
		itemp2 + 1, apot_table.param_name[itemp][itemp2],
		itemp + 1);
#else
	sprintf(errmsg,
		"F does not depend on the %d. parameter (%s) of the %d. potential (%s).\nFit impossible!\n",
		itemp2 + 1, apot_table.param_name[itemp][itemp2],
		itemp + 1, apot_table.names[itemp]);
#endif /* DIPOLE */
#endif
	warning(errmsg);
	break;
      }
    }
    (void)lineqsys_init(gamma, lineqsys, fxi1, p, ndim, mdim);	/*init LES */
    F3 = F;
    breakflag = 0;
    do {			/*inner loop - only calculate changed rows/lines in gamma */
      /* (a) solve linear equation */

      /* All in one driver routine */
      j = 1;			/* 1 rhs */
      /* Linear Equation Solution (lapack) */
#ifdef ACML
      dsysvx('N', 'U', ndim, j, &lineqsys[0][0], ndim,
	     &les_inverse[0][0], ndim, perm_indx, p, ndim, q, ndim,
	     &cond, &ferror, &berror, &i);
#else /* ACML */
      dsysvx(fact, uplo, &ndim, &j, &lineqsys[0][0], &ndim,
	     &les_inverse[0][0], &ndim, perm_indx, p, &ndim, q, &ndim,
	     &cond, &ferror, &berror, work, &worksize, iwork, &i);
#endif /* ACML */
#if defined DEBUG && !(defined APOT)
      printf("q0: %d %f %f %f %f %f %f %f %f\n", i, q[0], q[1], q[2],
	     q[3], q[4], q[5], q[6], q[7]);
#endif
      if (i > 0 && i <= ndim) {
	sprintf(errmsg, "Linear equation system singular after step %d i=%d",
		m, i);
	warning(errmsg);
	break;
      }
      /* (b) get delta by multiplying q with the direction vectors */
      for (i = 0; i < ndim; i++) {
	delta[idx[i]] = 0.;
	for (j = 0; j < ndim; j++)
	  delta[idx[i]] += d[i][j] * q[j];
#ifndef APOT
	if ((usemaxch) && (maxchange[idx[i]] > 0) &&
	    (fabs(delta[idx[i]]) > maxchange[idx[i]])) {
	  /* something seriously went wrong,
	     parameter idx[i] out of control */
	  sprintf(errmsg,
		  "Direction vector component %d out of range in step %d\n",
		  idx[i], m);
	  sprintf(errmsg, "%s(%g instead of %g).\n",
		  errmsg, fabs(delta[idx[i]]), maxchange[idx[i]]);
	  sprintf(errmsg, "%sRestarting inner loop\n", errmsg);
	  warning(errmsg);
	  breakflag = 1;
	}
#else
	if ((xi[idx[i]] + delta[idx[i]]) <
	    apot_table.pmin[apot_table.idxpot[i]][apot_table.idxparam[i]]) {
	  delta[idx[i]] =
	    apot_table.pmin[apot_table.idxpot[i]][apot_table.idxparam[i]] -
	    xi[idx[i]];
	}
	if ((xi[idx[i]] + delta[idx[i]]) >
	    apot_table.pmax[apot_table.idxpot[i]][apot_table.idxparam[i]]) {
	  delta[idx[i]] =
	    apot_table.pmax[apot_table.idxpot[i]][apot_table.idxparam[i]] -
	    xi[idx[i]];
	}
#endif
      }
      if (breakflag)
	break;
      /*     and store delta */
      copy_vector(delta, delta_norm, ndimtot);

      F2 = F;			/*shift F */

      /* (c) minimize F(xi) along vector delta, return new F */
      F = linmin(xi, delta, F, &xi1, &xi2, fxi1, fxi2);

#ifdef DEBUG
      printf("%f %6g %f %f %d\n", F, cond, ferror, berror, i);
#endif
      /* (d) if error estimate is too high after minimization
         in 5 directions: restart outer loop */
      if (ferror + berror > 1. && m > 5)
	break;

      /* (e) find optimal direction to replace */
      j = 0;
      temp2 = 0.;
      for (i = 0; i < ndim; i++)
	if ((temp = fabs(p[i] * q[i])) > temp2) {
	  j = i;
	  temp2 = temp;
	};

      /* (f) update gamma, but if fn returns 1, matrix will be sigular,
         break inner loop and restart with new matrix */
      if (gamma_update
	  (gamma, xi1, xi2, fxi1, fxi2, delta_norm, j, mdim, ndimtot, F)) {
	sprintf(errmsg,
		"Matrix gamma singular after step %d, restarting inner loop",
		m);
	warning(errmsg);
	break;
      }

      /* (g) set new direction vector */
      for (i = 0; i < ndim; i++)
	d[i][j] = delta_norm[idx[i]];

      /* (h) update linear equation system */
      lineqsys_update(gamma, lineqsys, fxi1, p, j, ndim, mdim);

      m++;			/*increment loop counter */
      df = F2 - F;

      /* loop at least ndim times, but at most INNERLOOPS or until no
         further improvement */
    } while ((m < ndim || (m <= INNERLOOPS && df > PRECISION))
	     && df < TOOBIG);
    /* inner loop */

    n++;			/* increment outer loop counter */

    /* Print the steps in current loop, F, a few values of xi, and
       total number of fn calls */
#ifdef APOT
    printf("%5d\t%17.6f\t%6d\n", m, F, fcalls);
#else
    printf("%d %f %f %f %f %f %f %d\n",
	   m, F, xi[0], xi[1], xi[2], xi[3], xi[4], fcalls);
#endif
    fflush(stdout);

    /* End fit if break flagfile exists */
    if (*flagfile != '\0') {
      ff = fopen(flagfile, "r");
      if (NULL != ff) {
	printf
	  ("Fit terminated prematurely in presence of break flagfile \"%s\"!\n",
	   flagfile);
	fclose(ff);
	remove(flagfile);
	break;
      }
    }
#ifdef EAM
#ifndef NORESCALE
    /* Check for rescaling... every fourth step */
    if ((n % 4) == 0) {
      temp = rescale(&opt_pot, 1., 0);
      /* Was rescaling necessary ? */
      if (temp != 0.) {
	/* wake other threads and sync potentials */
	F = calc_forces(xi, fxi1, 2);
      }
    }
#endif /* NORESCALE */
#endif /* EAM */
    /* write temp file  */
    if (*tempfile != '\0')
#ifndef APOT
      write_pot_table(&opt_pot, tempfile);	/*emergency writeout */
#else
      for (i = 0; i < ndim; i++)
	apot_table.values[apot_table.idxpot[i]][apot_table.idxparam[i]] =
	  xi[idx[i]];
    write_pot_table(&apot_table, tempfile);
#endif

    /*End fit if whole series didn't improve F */
  } while (((F3 - F > PRECISION / 10.) || (F3 - F < 0)) && (F3 - F > d_eps));
  /* outer loop */

  if (fabs(F3 - F) < PRECISION && F3 != F)
    printf("Precision reached: %10g\n", F3 - F);
  else if (F3 == F)
    printf("Could not find any further improvements, aborting!\n");
  else if ((fabs(F3 - F) > PRECISION && F3 != F && fabs(F3 - F) < d_eps))
    printf("Last improvement was smaller than d_eps (%f), aborting!\n",
	   d_eps);
  else
    printf("Precision not reached!\n");
#ifdef APOT
  for (i = 0; i < ndim; i++)
    apot_table.values[apot_table.idxpot[i]][apot_table.idxparam[i]] =
      xi[idx[i]];
#endif

  /* Free memory */
  free_vect_real(delta);
  free_vect_real(fxi1);
  free_vect_real(fxi2);
  free_mat_real(d);
  free_mat_real(gamma);
  free_mat_real(lineqsys);
  free_mat_real(les_inverse);
  free_vect_int(perm_indx);
  free_vect_real(delta_norm);
  free_vect_real(force_xi);
  free_vect_real(p);
  free_vect_real(q);
#ifndef ACML
  free(work);
  free(iwork);
#endif /* ACML */
  return;
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
  static real *force;
  int   i, j;			/* Auxiliary vars: Counters */
  real  sum, temp, scale, store;	/* Auxiliary var: Sum */
/*   Set direction vectors to coordinate directions d_ij=KroneckerDelta_ij */
  /*Initialize direction vectors */
  for (i = 0; i < ndim; i++) {
    for (j = 0; j < ndim; j++)
      d[i][j] = (i == j) ? 1. : 0.;
  }
/* Initialize gamma by calculating numerical derivatives    */
  if (force == NULL) {
    force = (real *)malloc(mdim * sizeof(real));
    if (force == NULL)
      error("Error in real vector allocation");
    for (i = 0; i < mdim; i++)
      force[i] = 0;
    reg_for_free(force, "force from init_gamma");
  }

  for (i = 0; i < ndim; i++) {	/*initialize gamma */
    store = xi[idx[i]];
#ifdef APOT
    scale = apot_table.pmax[apot_table.idxpot[i]][apot_table.idxparam[i]] -
      apot_table.pmin[apot_table.idxpot[i]][apot_table.idxparam[i]];
    xi[idx[i]] += (EPS * scale);
#else
    scale = 1.;
    xi[idx[i]] += EPS;		/*increase xi[idx[i]]... */
#endif
    sum = 0.;
    (void)(*calc_forces) (xi, force, 0);
    for (j = 0; j < mdim; j++) {
      temp = (force[j] - force_xi[j]) / (EPS * scale);
      gamma[j][i] = temp;
      sum += SQR(temp);
    }
    temp = sqrt(sum);
    xi[idx[i]] = store;		/*...and reset [idx[i]] again */
/* scale gamma so that sum_j(gamma^2)=1                      */
    if (temp > NOTHING) {
      for (j = 0; j < mdim; j++)
	gamma[j][i] /= temp;	/*normalize gamma */
      d[i][i] /= temp;		/* rescale d */
    } else
      return i + 1;		/* singular matrix, abort */
  }
  return 0;
}

/*******************************************************************
*
* gamma_update: Update column j of gamma ( to newly calculated
*           numerical derivatives (calculated from fa, fb
*           at a,b); normalize new vector.
*
*******************************************************************/

int gamma_update(real **gamma, real a, real b, real *fa, real *fb,
		 real *delta, int j, int m, int n, real fmin)
{
  int   i;
  real  temp;
  real  sum = 0.;
  real  mu = 0.;
  for (i = 0; i < m; i++) {
    temp = ((fa[i] - fb[i]) / (a - b));
    gamma[i][j] = temp;
    mu += temp * fa[i];
  }
  mu /= fmin;
  for (i = 0; i < m; i++) {
    temp = gamma[i][j] - mu * fa[i];
    gamma[i][j] = temp;
    sum += temp * temp;
  }
  temp = sqrt(sum);		/* normalization factor */
  if (temp > NOTHING) {
    for (i = 0; i < m; i++)
      gamma[i][j] /= temp;
    for (i = 0; i < n; i++)
      delta[i] /= temp;
  } else
    return 1;			/*Matrix will be singular: Restart! */
  return 0;
}

/*******************************************************************
*
* lineqsys_init: Initialize LinEqSys matrix, vector p in
*              lineqsys . q == p
*
*******************************************************************/

void lineqsys_init(real **gamma, real **lineqsys, real *deltaforce,
		   real *p, int n, int m)
{
  int   i, j, k;		/* Auxiliary vars: Counters */
/*   real  temp; */
  /* calculating vector p (lineqsys . q == P in LinEqSys) */

  for (i = 0; i < n; i++) {
    p[i] = 0.;
    for (j = 0; j < m; j++) {
      p[i] -= gamma[j][i] * deltaforce[j];
    }
  }
  /* calculating the linear equation system matrix gamma^t.gamma */
  for (i = 0; i < n; i++) {
    lineqsys[i][i] = 0;
    for (j = 0; j < m; j++)
      lineqsys[i][i] += SQR(gamma[j][i]);
    for (k = i + 1; k < n; k++) {
      lineqsys[i][k] = 0.;
      for (j = 0; j < m; j++) {
	lineqsys[i][k] += gamma[j][i] * gamma[j][k];
      }
      lineqsys[k][i] = lineqsys[i][k];
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
		     real *p, int i, int n, int m)
{
  int   k;
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
/*     int j,k; */
/*     real temp; */
#ifdef _OPENMP
#pragma omp for
#endif
    for (k = 0; k < n; k++) {
      int   j;
/*       real  temp; */
      p[k] = 0.;
      lineqsys[i][k] = 0.;
      for (j = 0; j < m; j++) {
	p[k] -= gamma[j][k] * force_xi[j];
	lineqsys[i][k] += gamma[j][i] * gamma[j][k];
      }
      lineqsys[k][i] = lineqsys[i][k];
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

void copy_matrix(real **a, real **b, int n, int m)
{
  int   i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      b[j][i] = a[j][i];
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

void copy_vector(real *a, real *b, int n)
{
  int   i;
  for (i = 0; i < n; i++)
    b[i] = a[i];
  return;
}

/******************************************************************
*
* matdotvec: Calculates the product of matrix a (n x m) with column
* vector x (dim m), Result y (dim n). (A . x = m)
*
******************************************************************/

void matdotvec(real **a, real *x, real *y, int n, int m)
{
  int   i, j;
  for (i = 0; i < n; i++) {
    y[idx[i]] = 0.;
    for (j = 0; j < m; j++)
      y[idx[i]] += a[i][j] * x[j];
  }
  return;
}

/*******************************************************************
*
* normalize_vector: Normalizes vector to |vec|^2=1, returns norm of old vector
*
*******************************************************************/

real normalize_vector(real *v, int n)
{
  int   j;
  real  temp, sum = 0.0;
  for (j = 0; j < n; j++)
    sum += SQR(v[j]);
  temp = sqrt(sum);
  for (j = 0; j < n; j++)
    v[j] /= temp;
  return temp;
}
