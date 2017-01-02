/****************************************************************
 *
 * powell_lsq.c: Contains Powell optimization algorithm
 *               and some small subroutines
 *
 ****************************************************************
 *
 * Copyright 2002-2017 - the potfit development team
 *
 * https://www.potfit.net/
 *
 ****************************************************************
 *
 * This file is part of potfit.
 *
 * potfit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * potfit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with potfit; if not, see <http://www.gnu.org/licenses/>.
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

#if defined(ACML)
#include <acml.h>
#else
#include <mkl_lapack.h>
#endif  // ACML

#include "bracket.h"
#include "force.h"
#include "memory.h"
#include "optimize.h"
#include "potential_input.h"
#include "potential_output.h"
#include "rescale.h"
#include "utils.h"

#define EPS 0.001
#define PRECISION 1.E-7
#define VERY_SMALL 1.E-12
#define INNERLOOPS 801
#define TOOBIG 10000

int gamma_init(double**, double**, double*, double*);
int gamma_update(double**, double, double, double*, double*, double*, int, int,
                 int, double);
void lineqsys_init(double**, double**, double*, double*, int, int);
void lineqsys_update(double**, double**, double*, double*, int, int, int);
double normalize_vector(double*, int);

double** mat_double(int rowdim, int coldim)
{
  double** matrix = NULL;

  /* matrix: array of array of pointers */
  /* matrix: pointer to rows */
  matrix = (double**)Malloc(rowdim * sizeof(double*));

  /* matrix[0]: pointer to elements */
  matrix[0] = (double*)Malloc(rowdim * coldim * sizeof(double));

  for (int i = 1; i < rowdim; i++)
    matrix[i] = matrix[i - 1] + coldim;

  return matrix;
}

/****************************************************************
 *
 * run_powell_lsq
 *
 ****************************************************************/

void run_powell_lsq(double* xi)
{
  int n = 0;
#if !defined(ACML)
  char uplo[1] = "U"; /* char used in dsysvx */
  char fact[1] = "N"; /* char used in dsysvx */
#endif                // ACML
  int breakflag;
  double cond = 0.0;
  double F1 = 0.0;
  double F2 = 0.0;
  double F3 = 0.0;
  double df = 0.0;
  double xi1 = 0.0;
  double xi2 = 0.0;
  double ferror = 0.0;
  double berror = 0.0; /* forward/backward error estimates */

  /* Direction vectors */
  double** d = mat_double(g_calc.ndim, g_calc.ndim);

  /* Matrix of derivatives */
  double** gamma = mat_double(g_calc.mdim, g_calc.ndim);

  /* Lin.Eq.Sys. Matrix */
  double** lineqsys = mat_double(g_calc.ndim, g_calc.ndim);

  /* LU decomp. of the lineqsys */
  double** les_inverse = mat_double(g_calc.ndim, g_calc.ndim);

  /* Keeps track of LU pivoting */
  int* perm_indx = (int*)Malloc(g_calc.ndim * sizeof(int));

  /* Normalized vector delta */
  double* delta_norm = (double*)Malloc(g_calc.ndimtot * sizeof(double));

  /* calculated forces */
  double* forces_1 = (double*)Malloc(g_calc.mdim * sizeof(double));
  double* forces_2 = (double*)Malloc(g_calc.mdim * sizeof(double));

  /* Vectors needed in Powell's algorithm */
  double* p = (double*)Malloc(g_calc.ndim * sizeof(double));
  double* q = (double*)Malloc(g_calc.ndim * sizeof(double));

  /* Vector pointing into correct dir'n */
  double* delta = (double*)Malloc(g_calc.ndimtot * sizeof(double)); /* ==0 */

#if !defined(ACML) /* work arrays not needed */
  int worksize = 64 * g_calc.ndim;
  /* work array to be used by dsysvx */
  double* work = (double*)Malloc(worksize * sizeof(double));
  int* iwork = (int*)Malloc(g_calc.ndim * sizeof(int));
#endif  // ACML

  /* calculate the first force */
  F1 = calc_forces(xi, forces_1, 0);

  if (F1 < VERY_SMALL) {
    printf("Error already too small to optimize, aborting ...\n");
    return;
  }

  memcpy(forces_2, forces_1, g_calc.mdim * sizeof(double));

#if defined(APOT)
  printf("loops\t\terror_sum\tforce calculations\n");
  printf("%5d\t%17.6f\t%6d\n", 0, F1, g_calc.fcalls);
#else
  printf("%d %f %f %f %f %f %f %d\n", 0, F1, xi[0], xi[1], xi[2], xi[3], xi[4],
         g_calc.fcalls);
#endif  // APOT
  fflush(stdout);

  do {
    /*outer loop, includes recalculating gamma */
    int m = 0;

    /* Init gamma */
    int i = gamma_init(gamma, d, xi, forces_2);

    if (i != 0) {
#if defined(RESCALE) && (defined(EAM) || defined(ADP) || defined(MEAM))
      /* perhaps rescaling helps? - Last resort... */
      warning("F does not depend on xi[%d], trying to rescale!\n",
              g_pot.opt_pot.idx[i - 1]);

      rescale(&g_pot.opt_pot, 1.0, 1);

      /* wake other threads and sync potentials */
      F1 = calc_forces(xi, forces_1, 2);

      i = gamma_init(gamma, d, xi, forces_1);
#endif  // RESCALE && ( EAM || ADP || MEAM )

      /* try again */
      if (i != 0) {
/* ok, now this is serious, better exit cleanly */
#if !defined(APOT)
        write_pot_table_potfit(g_files.tempfile); /*emergency writeout */
        warning("F does not depend on xi[%d], fit impossible!\n",
                g_pot.opt_pot.idx[i - 1]);
#else
        update_apot_table(xi);
        write_pot_table_potfit(g_files.tempfile);
        warning(
            "F does not depend on the %d. parameter (%s) of the %d. "
            "potential.\n",
            g_pot.apot_table.idxparam[i - 1] + 1,
            g_pot.apot_table.param_name[g_pot.apot_table.idxpot[i - 1]]
                                       [g_pot.apot_table.idxparam[i - 1]],
            g_pot.apot_table.idxpot[i - 1] + 1);
        warning("Fit impossible!\n");
#endif  // APOT
        break;
      }
    }

    /*init LES */
    lineqsys_init(gamma, lineqsys, forces_1, p, g_calc.ndim, g_calc.mdim);

    F3 = F1;

    breakflag = 0;

    /*inner loop - only calculate changed rows/lines in gamma */
    do {
      /* (a) solve linear equation */

      /* All in one driver routine */
      int j = 1; /* 1 rhs */

/* Linear Equation Solution (lapack) */
#if defined(ACML)
      dsysvx('N', 'U', g_calc.ndim, j, &lineqsys[0][0], g_calc.ndim,
             &les_inverse[0][0], g_calc.ndim, perm_indx, p, g_calc.ndim, q,
             g_calc.ndim, &cond, &ferror, &berror, &i);
#else
      dsysvx(fact, uplo, &g_calc.ndim, &j, &lineqsys[0][0], &g_calc.ndim,
             &les_inverse[0][0], &g_calc.ndim, perm_indx, p, &g_calc.ndim, q,
             &g_calc.ndim, &cond, &ferror, &berror, work, &worksize, iwork, &i);
#endif  // ACML

#if defined(DEBUG) && !(defined APOT)
      printf("q0: %d %f %f %f %f %f %f %f %f\n", i, q[0], q[1], q[2], q[3],
             q[4], q[5], q[6], q[7]);
#endif  // DEBUG && !APOT

      if (i > 0 && i <= g_calc.ndim) {
        warning("Linear equation system singular after step %d i=%d\n", m, i);
        break;
      }

      /* (b) get delta by multiplying q with the direction vectors */
      for (i = 0; i < g_calc.ndim; i++) {
        delta[g_pot.opt_pot.idx[i]] = 0.0;

        for (j = 0; j < g_calc.ndim; j++)
          delta[g_pot.opt_pot.idx[i]] += d[i][j] * q[j];

#if !defined(APOT)
        if ((g_param.usemaxch) &&
            (g_calc.maxchange[g_pot.opt_pot.idx[i]] > 0) &&
            (fabs(delta[g_pot.opt_pot.idx[i]]) >
             g_calc.maxchange[g_pot.opt_pot.idx[i]])) {
          /* something seriously went wrong,
             parameter idx[i] out of control */
          warning("Direction vector component %d out of range in step %d\n",
                  g_pot.opt_pot.idx[i], m);
          warning("(%g instead of %g).\n", fabs(delta[g_pot.opt_pot.idx[i]]),
                  g_calc.maxchange[g_pot.opt_pot.idx[i]]);
          warning("Restarting inner loop\n");
          breakflag = 1;
        }
#else
        if ((xi[g_pot.opt_pot.idx[i]] + delta[g_pot.opt_pot.idx[i]]) <
            g_pot.apot_table.pmin[g_pot.apot_table.idxpot[i]]
                                 [g_pot.apot_table.idxparam[i]]) {
          delta[g_pot.opt_pot.idx[i]] =
              g_pot.apot_table.pmin[g_pot.apot_table.idxpot[i]]
                                   [g_pot.apot_table.idxparam[i]] -
              xi[g_pot.opt_pot.idx[i]];
        }
        if ((xi[g_pot.opt_pot.idx[i]] + delta[g_pot.opt_pot.idx[i]]) >
            g_pot.apot_table.pmax[g_pot.apot_table.idxpot[i]]
                                 [g_pot.apot_table.idxparam[i]]) {
          delta[g_pot.opt_pot.idx[i]] =
              g_pot.apot_table.pmax[g_pot.apot_table.idxpot[i]]
                                   [g_pot.apot_table.idxparam[i]] -
              xi[g_pot.opt_pot.idx[i]];
        }
#endif  // !APOT
      }

      if (breakflag)
        break;

      /*     and store delta */
      memcpy(delta_norm, delta, g_calc.ndimtot * sizeof(double));

      F2 = F1; /*shift F */

      /* (c) minimize F(xi) along vector delta, return new F */
      F1 = linmin(xi, delta, F1, &xi1, &xi2, forces_1, forces_2);

#if defined(DEBUG)
      printf("%f %6g %f %f %d\n", F1, cond, ferror, berror, i);
#endif  // DEBUG

      /* (d) if error estimate is too high after minimization
         in 5 directions: restart outer loop */
      if (ferror + berror > 1.0 && m > 5)
        break;

      /* (e) find optimal direction to replace */
      j = 0;
      double temp2 = 0.0;

      for (i = 0; i < g_calc.ndim; i++) {
        double temp = fabs(p[i] * q[i]);
        if (temp > temp2) {
          j = i;
          temp2 = temp;
        }
      }

      /* (f) update gamma, but if fn returns 1, matrix will be sigular,
         break inner loop and restart with new matrix */
      if (gamma_update(gamma, xi1, xi2, forces_1, forces_2, delta_norm, j,
                       g_calc.mdim, g_calc.ndimtot, F1)) {
        warning("Matrix gamma singular after step %d, restarting inner loop\n",
                m);
        break;
      }

      /* (g) set new direction vector */
      for (i = 0; i < g_calc.ndim; i++)
        d[i][j] = delta_norm[g_pot.opt_pot.idx[i]];

      /* (h) update linear equation system */
      lineqsys_update(gamma, lineqsys, forces_1, p, j, g_calc.ndim,
                      g_calc.mdim);

      m++; /*increment loop counter */
      df = F2 - F1;

      /* loop at least ndim times, but at most INNERLOOPS or until no
         further improvement */
    } while ((m < g_calc.ndim || (m <= INNERLOOPS && df > PRECISION)) &&
             df < TOOBIG);
    /* inner loop */

    n++; /* increment outer loop counter */

/* Print the steps in current loop, F, a few values of xi, and total number of
 * fn calls */
#if defined(APOT)
    printf("%5d\t%17.6f\t%6d\n", m, F1, g_calc.fcalls);
#else
    printf("%d %f %f %f %f %f %f %d\n", m, F1, xi[0], xi[1], xi[2], xi[3],
           xi[4], g_calc.fcalls);
#endif  // APOT
    fflush(stdout);

    /* End fit if break flagfile exists */
    if (*g_files.flagfile != '\0') {
      FILE* ff = fopen(g_files.flagfile, "r");
      if (ff != NULL) {
        printf(
            "Fit terminated prematurely in presence of break flagfile "
            "\"%s\"!\n",
            g_files.flagfile);
        fclose(ff);
        remove(g_files.flagfile);
        break;
      }
    }

/* WARNING: This rescaling is not necessary for EAM. Causes more problems. */
#if defined(RESCALE) && (defined(xEAM) || defined(xMEAM))
    /* Check for rescaling... every fourth step */
    if ((n % 4) == 0) {
      temp = rescale(&opt_pot, 1.0, 0);

      /* Was rescaling necessary ? */
      if (temp != 0.0) {
        /* wake other threads and sync potentials */
        F1 = calc_forces(xi, forces_1, 2);
      }
    }
#endif  // RESCALE && ( xEAM || xMEAM )

    /* write temp file  */
    if (*g_files.tempfile != '\0') {
#if defined(APOT)
      update_apot_table(xi);
#endif  // APOT
      write_pot_table_potfit(g_files.tempfile);
    }

    /*End fit if whole series didn't improve F */
  } while (((F3 - F1 > PRECISION / 10.0) || (F3 - F1 < 0)) &&
           (F3 - F1 > g_calc.d_eps));
  /* outer loop */

  if (fabs(F3 - F1) < PRECISION && F3 != F1)
    printf("Precision reached: %10g\n", F3 - F1);
  else if (F3 == F1)
    printf("Could not find any further improvements, aborting!\n");
  else if ((fabs(F3 - F1) > PRECISION && F3 != F1 &&
            fabs(F3 - F1) < g_calc.d_eps))
    printf("Last improvement was smaller than d_eps (%f), aborting!\n",
           g_calc.d_eps);
  else
    printf("Precision not reached!\n");

#if defined(APOT)
  update_apot_table(xi);
#endif  // APOT
}

/****************************************************************
 *
 * gamma_init: (Re-)Initialize gamma[j][i] (Gradient Matrix) after
 *            (Re-)Start or whenever necessary by calculating numerical
 *            gradients in coordinate directions. Includes re-setting the
 *            direction vectors to coordinate directions.
 *
 ****************************************************************/

int gamma_init(double** gamma, double** d, double* xi, double* force_xi)
{
  static double* force;

  double sum, temp, scale, store; /* Auxiliary var: Sum */

  /* Set direction vectors to coordinate directions d_ij=KroneckerDelta_ij */
  for (int i = 0; i < g_calc.ndim; i++)
    for (int j = 0; j < g_calc.ndim; j++)
      d[i][j] = (i == j) ? 1.0 : 0.0;

  /* Initialize gamma by calculating numerical derivatives */
  if (force == NULL)
    force = (double*)Malloc(g_calc.mdim * sizeof(double));

  /*initialize gamma */
  for (int i = 0; i < g_calc.ndim; i++) {
    store = xi[g_pot.opt_pot.idx[i]];
#if defined(APOT)
    scale =
        g_pot.apot_table
            .pmax[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]] -
        g_pot.apot_table
            .pmin[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]];
    xi[g_pot.opt_pot.idx[i]] += (EPS * scale);
#else
    scale = 1.0;
    xi[g_pot.opt_pot.idx[i]] += EPS; /*increase xi[idx[i]]... */
#endif  // APOT

    sum = 0.0;

    calc_forces(xi, force, 0);

    for (int j = 0; j < g_calc.mdim; j++) {
      temp = (force[j] - force_xi[j]) / (EPS * scale);
      gamma[j][i] = temp;
      sum += dsquare(temp);
    }

    temp = sqrt(sum);

    xi[g_pot.opt_pot.idx[i]] = store; /*...and reset [idx[i]] again */

    /* scale gamma so that sum_j(gamma^2)=1                      */
    if (temp > VERY_SMALL) {
      for (int j = 0; j < g_calc.mdim; j++)
        gamma[j][i] /= temp; /*normalize gamma */
      d[i][i] /= temp;       /* rescale d */
    } else
      return i + 1; /* singular matrix, abort */
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

int gamma_update(double** gamma, double a, double b, double* fa, double* fb,
                 double* delta, int j, int m, int n, double fmin)
{
  double temp;
  double sum = 0.0;
  double mu = 0.0;

  for (int i = 0; i < m; i++) {
    temp = ((fa[i] - fb[i]) / (a - b));
    gamma[i][j] = temp;
    mu += temp * fa[i];
  }

  mu /= fmin;

  for (int i = 0; i < m; i++) {
    temp = gamma[i][j] - mu * fa[i];
    gamma[i][j] = temp;
    sum += temp * temp;
  }

  temp = sqrt(sum); /* normalization factor */

  if (temp > VERY_SMALL) {
    for (int i = 0; i < m; i++)
      gamma[i][j] /= temp;
    for (int i = 0; i < n; i++)
      delta[i] /= temp;
  } else
    return 1; /* Matrix will be singular: Restart! */

  return 0;
}

/****************************************************************
 *
 * lineqsys_init: Initialize LinEqSys matrix, vector p in
 *              lineqsys . q == p
 *
 ****************************************************************/

void lineqsys_init(double** gamma, double** lineqsys, double* deltaforce,
                   double* p, int n, int m)
{
  /* calculating vector p (lineqsys . q == P in LinEqSys) */
  for (int i = 0; i < n; i++) {
    p[i] = 0.0;
    for (int j = 0; j < m; j++) {
      p[i] -= gamma[j][i] * deltaforce[j];
    }
  }

  /* calculating the linear equation system matrix gamma^t.gamma */
  for (int i = 0; i < n; i++) {
    lineqsys[i][i] = 0;
    for (int j = 0; j < m; j++)
      lineqsys[i][i] += dsquare(gamma[j][i]);
    for (int k = i + 1; k < n; k++) {
      lineqsys[i][k] = 0.0;
      for (int j = 0; j < m; j++) {
        lineqsys[i][k] += gamma[j][i] * gamma[j][k];
      }
      lineqsys[k][i] = lineqsys[i][k];
    }
  }
}

/****************************************************************
 *
 * lineqsys_update: Update LinEqSys matrix row and column i, vector
 *            p.
 *
 ****************************************************************/

void lineqsys_update(double** gamma, double** lineqsys, double* force_xi,
                     double* p, int i, int n, int m)
{
  for (int k = 0; k < n; k++) {
    p[k] = 0.0;

    lineqsys[i][k] = 0.0;

    for (int j = 0; j < m; j++) {
      p[k] -= gamma[j][k] * force_xi[j];
      lineqsys[i][k] += gamma[j][i] * gamma[j][k];
    }
    lineqsys[k][i] = lineqsys[i][k];
  }
}
