/****************************************************************
 *
 * powell_lsq.c: Contains Powell optimization algorithm
 *	and some small subroutines
 *
 ****************************************************************
 *
 * Copyright 2002-2015
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.sourceforge.net/
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

/****************************************************************
 *
 *  This utility will optimize a parameter vector xi[N] to
 *  minimize the sum of squares U=sum_{i=1..M}(f_i(xi)-F_i)^2,
 *  where F_i is an array passed to the utility, and f_i(xi)
 *  is a function of the vector xi.
 *
 ****************************************************************/

#include "potfit.h"

#ifdef ACML
#include <acml.h>
#else /* ACML */
#include <mkl_lapack.h>
#endif /* ACML */

#include "bracket.h"
#include "forces.h"
#include "optimize.h"
#include "potential_input.h"
#include "potential_output.h"
#include "utils.h"

#define EPS 0.001
#define PRECISION 1.E-7
#define NOTHING 1.E-12		/* Well, almost nothing */
#define INNERLOOPS 801
#define TOOBIG 10000

/* powell least squares [powell_lsq.c] */
void  powell_lsq(double *);
int   gamma_init(double **, double **, double *, double *);
int   gamma_update(double **, double, double, double *, double *, double *, int, int, int, double);
void  lineqsys_init(double **, double **, double *, double *, int, int);
void  lineqsys_update(double **, double **, double *, double *, int, int, int);
void  copy_matrix(double **, double **, int, int);
void  copy_vector(double *, double *, int);
void  matdotvec(double **, double *, double *, int, int);
double normalize_vector(double *, int);

void run_powell_lsq(double *xi)
{
#ifndef ACML
  char  uplo[1] = "U";		/* char used in dsysvx */
  char  fact[1] = "N";		/* char used in dsysvx */
#endif /* ACML */
  int   i, j, m = 0, n = 0;	/* Simple counting variables */
  double *force_xi;		/* calculated force, alt */
  double **d;			/* Direction vectors */
  double **gamma;		/* Matrix of derivatives */
  double **lineqsys;		/* Lin.Eq.Sys. Matrix */
  double **les_inverse;		/* LU decomp. of the lineqsys */
  double *delta;		/* Vector pointing into correct dir'n */
  double *delta_norm;		/* Normalized vector delta */
  double *fxi1, *fxi2;		/* two latest force vectors */
#ifndef ACML			/* work arrays not needed for ACML */
  double *work;			/* work array to be used by dsysvx */
  int  *iwork;
  int   worksize;		/* Size of work array (dsysvx) */
#endif /* ACML */
  int  *perm_indx;		/* Keeps track of LU pivoting */
  int   breakflag;		/* Breakflag */
  double cond = 0.0;		/* Condition number dsysvx */
  double *p, *q;		/* Vectors needed in Powell's algorithm */
  double F, F2, F3 = 0, df, xi1, xi2;	/* Fn values, changes, steps ... */
  double temp, temp2;		/* as the name indicates: temporary vars */
#ifdef APOT
  int   itemp, itemp2;		/* the same for integer */
#endif /* APOT */
  double ferror = 0.0;
  double berror = 0.0;		/* forward/backward error estimates */
  FILE *ff;			/* Exit flagfile */

  d = mat_double(g_calc.ndim, g_calc.ndim);
  gamma = mat_double(g_calc.mdim, g_calc.ndim);
  lineqsys = mat_double(g_calc.ndim, g_calc.ndim);
  les_inverse = mat_double(g_calc.ndim, g_calc.ndim);
  perm_indx = vect_int(g_calc.ndim);
  delta_norm = vect_double(g_calc.ndimtot);
				 /*==0*/
  force_xi = vect_double(g_calc.mdim);
  p = vect_double(g_calc.ndim);
  q = vect_double(g_calc.ndim);
  delta = vect_double(g_calc.ndimtot);	/* ==0 */
  fxi1 = vect_double(g_calc.mdim);
  fxi2 = vect_double(g_calc.mdim);
#ifndef ACML			/* work arrays not needed */
  worksize = 64 * g_calc.ndim;
  work = (double *)malloc(worksize * sizeof(double));
  iwork = (int *)malloc(g_calc.ndim * sizeof(int));
#endif /* ACML */

  /* clear delta */
  for (i = 0; i < g_calc.ndimtot; i++)
    delta[i] = 0.0;

  /* calculate the first force */
  F = g_calc_forces(xi, fxi1, 0);
#ifndef APOT
  printf("%d %f %f %f %f %f %f %d\n", m, F, xi[0], xi[1], xi[2], xi[3], xi[4], g_calc.fcalls);
  fflush(stdout);
#endif /* APOT */

  if (F < NOTHING) {
    printf("Error already too small to optimize, aborting ...\n");
    return;			/* If F is less than nothing, */
    /* what is there to do? */
  }

  copy_vector(fxi1, force_xi, g_calc.mdim);
#if defined(APOT)
  printf("loops\t\terror_sum\tforce calculations\n");
  printf("%5d\t%17.6f\t%6d\n", m, F, g_calc.fcalls);
#endif /* APOT */

  do {				/*outer loop, includes recalculating gamma */
    m = 0;

    /* Init gamma */
    i = gamma_init(gamma, d, xi, fxi1);
    if (0 != i) {
#ifdef RESCALE
#if defined EAM || defined ADP || defined MEAM
      /* perhaps rescaling helps? - Last resort... */
      warning("F does not depend on xi[%d], trying to rescale!\n", g_todo.idx[i - 1]);
      rescale(&g_pot.opt_pot, 1.0, 1);
      /* wake other threads and sync potentials */
      F = (*g_calc_forces)(xi, fxi1, 2);
      i = gamma_init(gamma, d, xi, fxi1);
#endif /* EAM */
#endif /* RESCALE */

      /* try again */
      if (0 != i) {
	/* ok, now this is serious, better exit cleanly */
#ifndef APOT
	write_pot_table_potfit(g_files.tempfile);	/*emergency writeout */
	warning("F does not depend on xi[%d], fit impossible!\n", g_todo.idx[i - 1]);
#else
	update_apot_table(xi);
        write_pot_table_potfit(g_files.tempfile);
	itemp = g_pot.apot_table.idxpot[i - 1];
        itemp2 = g_pot.apot_table.idxparam[i - 1];
	warning("F does not depend on the %d. parameter (%s) of the %d. potential.\n",
                itemp2 + 1, g_pot.apot_table.param_name[itemp][itemp2], itemp + 1);
	warning("Fit impossible!\n");
#endif /* APOT */
	break;
      }
    }
    lineqsys_init(gamma, lineqsys, fxi1, p, g_calc.ndim, g_calc.mdim);	/*init LES */
    F3 = F;
    breakflag = 0;

    /*inner loop - only calculate changed rows/lines in gamma */
    do {
      /* (a) solve linear equation */

      /* All in one driver routine */
      j = 1;			/* 1 rhs */

      /* Linear Equation Solution (lapack) */
#ifdef ACML
dsysvx('N', 'U', g_calc.ndim, j, &lineqsys[0][0], g_calc.ndim, &les_inverse[0][0], g_calc.ndim,
       perm_indx, p, g_calc.ndim, q, g_calc.ndim, &cond, &ferror, &berror, &i);
#else
      dsysvx(fact, uplo, &g_calc.ndim, &j, &lineqsys[0][0], &g_calc.ndim, &les_inverse[0][0],
             &g_calc.ndim, perm_indx, p, &g_calc.ndim, q, &g_calc.ndim, &cond, &ferror, &berror, work, &worksize, iwork, &i);
#endif /* ACML */
#if defined DEBUG && !(defined APOT)
      printf("q0: %d %f %f %f %f %f %f %f %f\n", i, q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7]);
#endif /* DEBUG && !APOT */
      if (i > 0 && i <= g_calc.ndim) {
	warning("Linear equation system singular after step %d i=%d\n", m, i);
	break;
      }
      /* (b) get delta by multiplying q with the direction vectors */
      for (i = 0; i < g_calc.ndim; i++) {
        delta[g_todo.idx[i]] = 0.0;
        for (j = 0; j < g_calc.ndim; j++)
          delta[g_todo.idx[i]] += d[i][j] * q[j];
#ifndef APOT
        if ((g_param.usemaxch) && (g_todo.maxchange[g_todo.idx[i]] > 0)
          && (fabs(delta[g_todo.idx[i]]) > g_todo.maxchange[g_todo.idx[i]])) {
	  /* something seriously went wrong,
	     parameter idx[i] out of control */
	  warning("Direction vector component %d out of range in step %d\n", g_todo.idx[i], m);
	  warning("(%g instead of %g).\n", fabs(delta[g_todo.idx[i]]), g_todo.maxchange[g_todo.idx[i]]);
	  warning("Restarting inner loop\n");
	  breakflag = 1;
	}
#else
if ((xi[g_todo.idx[i]] + delta[g_todo.idx[i]]) < g_pot.apot_table.pmin[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]]) {
  delta[g_todo.idx[i]] = g_pot.apot_table.pmin[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]] - xi[g_todo.idx[i]];
	}
	if ((xi[g_todo.idx[i]] + delta[g_todo.idx[i]]) > g_pot.apot_table.pmax[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]]) {
          delta[g_todo.idx[i]] = g_pot.apot_table.pmax[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]] - xi[g_todo.idx[i]];
	}
#endif /* !APOT */
      }
      if (breakflag)
	break;
      /*     and store delta */
      copy_vector(delta, delta_norm, g_calc.ndimtot);

      F2 = F;			/*shift F */

      /* (c) minimize F(xi) along vector delta, return new F */
      F = linmin(xi, delta, F, &xi1, &xi2, fxi1, fxi2);

#ifdef DEBUG
      printf("%f %6g %f %f %d\n", F, cond, ferror, berror, i);
#endif /* DEBUG */

      /* (d) if error estimate is too high after minimization
         in 5 directions: restart outer loop */
      if (ferror + berror > 1.0 && m > 5)
	break;

      /* (e) find optimal direction to replace */
      j = 0;
      temp2 = 0.0;
      for (i = 0; i < g_calc.ndim; i++)
	if ((temp = fabs(p[i] * q[i])) > temp2) {
	  j = i;
	  temp2 = temp;
	};

      /* (f) update gamma, but if fn returns 1, matrix will be sigular,
         break inner loop and restart with new matrix */
      if (gamma_update(gamma, xi1, xi2, fxi1, fxi2, delta_norm, j, g_calc.mdim, g_calc.ndimtot, F)) {
	warning("Matrix gamma singular after step %d, restarting inner loop\n", m);
	break;
      }

      /* (g) set new direction vector */
      for (i = 0; i < g_calc.ndim; i++)
	d[i][j] = delta_norm[g_todo.idx[i]];

      /* (h) update linear equation system */
      lineqsys_update(gamma, lineqsys, fxi1, p, j, g_calc.ndim, g_calc.mdim);

      m++;			/*increment loop counter */
      df = F2 - F;

      /* loop at least ndim times, but at most INNERLOOPS or until no
         further improvement */
    } while ((m < g_calc.ndim || (m <= INNERLOOPS && df > PRECISION))
      && df < TOOBIG);
    /* inner loop */

    n++;			/* increment outer loop counter */

    /* Print the steps in current loop, F, a few values of xi, and
       total number of fn calls */
#ifdef APOT
printf("%5d\t%17.6f\t%6d\n", m, F, g_calc.fcalls);
#else
    printf("%d %f %f %f %f %f %f %d\n", m, F, xi[0], xi[1], xi[2], xi[3], xi[4], g_calc.fcalls);
#endif /* APOT */
    fflush(stdout);

    /* End fit if break flagfile exists */
    if (*g_files.flagfile != '\0') {
      ff = fopen(g_files.flagfile, "r");
      if (NULL != ff) {
        printf("Fit terminated prematurely in presence of break flagfile \"%s\"!\n", g_files.flagfile);
	fclose(ff);
        remove(g_files.flagfile);
	break;
      }
    }

    /* WARNING: This rescaling is not necessary for EAM. Causes more problems. */
#ifdef RESCALE
#if defined xEAM || defined xMEAM
    /* Check for rescaling... every fourth step */
    if ((n % 4) == 0) {
      temp = rescale(&opt_pot, 1.0, 0);
      /* Was rescaling necessary ? */
      if (temp != 0.0) {
	/* wake other threads and sync potentials */
	F = calc_forces(xi, fxi1, 2);
      }
    }
#endif /* xEAM || xMEAM */
#endif /* RESCALE */

    /* write temp file  */
    if (*g_files.tempfile != '\0') {
#ifndef APOT
      write_pot_table_potfit(g_files.tempfile);	/*emergency writeout */
#else
      update_apot_table(xi);
      write_pot_table_potfit(g_files.tempfile);
#endif /* APOT */
    }

    /*End fit if whole series didn't improve F */
  } while (((F3 - F > PRECISION / 10.0) || (F3 - F < 0)) && (F3 - F > g_calc.d_eps));
  /* outer loop */

  if (fabs(F3 - F) < PRECISION && F3 != F)
    printf("Precision reached: %10g\n", F3 - F);
  else if (F3 == F)
    printf("Could not find any further improvements, aborting!\n");
  else if ((fabs(F3 - F) > PRECISION && F3 != F && fabs(F3 - F) < g_calc.d_eps))
    printf("Last improvement was smaller than d_eps (%f), aborting!\n", g_calc.d_eps);
  else
    printf("Precision not reached!\n");
#ifdef APOT
  update_apot_table(xi);
#endif /* APOT */

  /* Free memory */
  free_vect_double(delta);
  free_vect_double(fxi1);
  free_vect_double(fxi2);
  free_mat_double(d);
  free_mat_double(gamma);
  free_mat_double(lineqsys);
  free_mat_double(les_inverse);
  free_vect_int(perm_indx);
  free_vect_double(delta_norm);
  free_vect_double(force_xi);
  free_vect_double(p);
  free_vect_double(q);
#ifndef ACML
  free(work);
  free(iwork);
#endif /* ACML */
  return;
}


/****************************************************************
 *
 * gamma_init: (Re-)Initialize gamma[j][i] (Gradient Matrix) after
 *            (Re-)Start or whenever necessary by calculating numerical
 *            gradients in coordinate directions. Includes re-setting the
 *            direction vectors to coordinate directions.
 *
 ****************************************************************/

int gamma_init(double **gamma, double **d, double *xi, double *force_xi)
{
  static double *force;
  int   i, j;			/* Auxiliary vars: Counters */
  double sum, temp, scale, store;	/* Auxiliary var: Sum */
/*   Set direction vectors to coordinate directions d_ij=KroneckerDelta_ij */
  /*Initialize direction vectors */
  for (i = 0; i < g_calc.ndim; i++) {
    for (j = 0; j < g_calc.ndim; j++)
      d[i][j] = (i == j) ? 1.0 : 0.0;
  }
/* Initialize gamma by calculating numerical derivatives    */
  if (force == NULL) {
    force = (double *)malloc(g_calc.mdim * sizeof(double));
    if (force == NULL)
      error(1, "Error in double vector allocation");
    for (i = 0; i < g_calc.mdim; i++)
      force[i] = 0;
    reg_for_free(force, "force from init_gamma");
  }

  for (i = 0; i < g_calc.ndim; i++) {	/*initialize gamma */
    store = xi[g_todo.idx[i]];
#ifdef APOT
    scale =
      g_pot.apot_table.pmax[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]] -
      g_pot.apot_table.pmin[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]];
    xi[g_todo.idx[i]] += (EPS * scale);
#else
    scale = 1.0;
    xi[g_todo.idx[i]] += EPS;		/*increase xi[idx[i]]... */
#endif /* APOT */
    sum = 0.0;
    (*g_calc_forces)(xi, force, 0);
    for (j = 0; j < g_calc.mdim; j++) {
      temp = (force[j] - force_xi[j]) / (EPS * scale);
      gamma[j][i] = temp;
      sum += dsquare(temp);
    }
    temp = sqrt(sum);
    xi[g_todo.idx[i]] = store;		/*...and reset [idx[i]] again */
/* scale gamma so that sum_j(gamma^2)=1                      */
    if (temp > NOTHING) {
      for (j = 0; j < g_calc.mdim; j++)
	gamma[j][i] /= temp;	/*normalize gamma */
      d[i][i] /= temp;		/* rescale d */
    } else
      return i + 1;		/* singular matrix, abort */
  }
  return 0;
}

/****************************************************************
 *
 * gamma_update: Update column j of gamma ( to newly calculated
 *           numerical derivatives (calculated from fa, fb
 *           at a,b); normalize new vector.
 *
 ****************************************************************/

int gamma_update(double **gamma, double a, double b, double *fa, double *fb,
  double *delta, int j, int m, int n, double fmin)
{
  int   i;
  double temp;
  double sum = 0.0;
  double mu = 0.0;
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
    return 1;			/* Matrix will be singular: Restart! */
  return 0;
}

/****************************************************************
 *
 * lineqsys_init: Initialize LinEqSys matrix, vector p in
 *              lineqsys . q == p
 *
 ****************************************************************/

void lineqsys_init(double **gamma, double **lineqsys, double *deltaforce, double *p, int n, int m)
{
  int   i, j, k;		/* Auxiliary vars: Counters */
/*   double  temp; */
  /* calculating vector p (lineqsys . q == P in LinEqSys) */

  for (i = 0; i < n; i++) {
    p[i] = 0.0;
    for (j = 0; j < m; j++) {
      p[i] -= gamma[j][i] * deltaforce[j];
    }
  }
  /* calculating the linear equation system matrix gamma^t.gamma */
  for (i = 0; i < n; i++) {
    lineqsys[i][i] = 0;
    for (j = 0; j < m; j++)
      lineqsys[i][i] += dsquare(gamma[j][i]);
    for (k = i + 1; k < n; k++) {
      lineqsys[i][k] = 0.0;
      for (j = 0; j < m; j++) {
	lineqsys[i][k] += gamma[j][i] * gamma[j][k];
      }
      lineqsys[k][i] = lineqsys[i][k];
    }
  }
  return;
}

/****************************************************************
 *
 * lineqsys_update: Update LinEqSys matrix row and column i, vector
 *            p.
 *
 ****************************************************************/

void lineqsys_update(double **gamma, double **lineqsys, double *force_xi, double *p, int i, int n, int m)
{
  int   j, k;
  for (k = 0; k < n; k++) {
    p[k] = 0.0;
    lineqsys[i][k] = 0.0;
    for (j = 0; j < m; j++) {
      p[k] -= gamma[j][k] * force_xi[j];
      lineqsys[i][k] += gamma[j][i] * gamma[j][k];
    }
    lineqsys[k][i] = lineqsys[i][k];
  }
  return;
}



/****************************************************************
 *
 *  copy_matrix: Copies data from Matrix a into matrix b
 *              (matrix dimension n x m)
 *
 ****************************************************************/

void copy_matrix(double **a, double **b, int n, int m)
{
  int   i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      b[j][i] = a[j][i];
    }
  }
  return;
}

/****************************************************************
 *
 * copy_vector: Copies data from vector a into vector b (both dim n)
 *
 ****************************************************************/

void copy_vector(double *a, double *b, int n)
{
  int   i;
  for (i = 0; i < n; i++)
    b[i] = a[i];
}

/****************************************************************
 *
 * matdotvec: Calculates the product of matrix a (n x m) with column
 * vector x (dim m), Result y (dim n). (A . x = m)
 *
 ****************************************************************/

void matdotvec(double **a, double *x, double *y, int n, int m)
{
  int   i, j;
  for (i = 0; i < n; i++) {
    y[g_todo.idx[i]] = 0.0;
    for (j = 0; j < m; j++)
      y[g_todo.idx[i]] += a[i][j] * x[j];
  }
}

/****************************************************************
 *
 * normalize_vector: Normalizes vector to |vec|^2=1, returns norm of old vector
 *
 ****************************************************************/

double normalize_vector(double *v, int n)
{
  int   j;
  double temp, sum = 0.0;
  for (j = 0; j < n; j++)
    sum += dsquare(v[j]);
  temp = sqrt(sum);
  for (j = 0; j < n; j++)
    v[j] /= temp;
  return temp;
}
