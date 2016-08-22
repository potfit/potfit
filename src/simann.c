/****************************************************************
 *
 * simann.c: Contains all routines used for simulated annealing.
 *
 *****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
 *
 * http://potfit.sourceforge.net/
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
 *****************************************************************/

#include "potfit.h"

#if !defined(EVO)

#include <ctype.h>

#include "force.h"
#include "memory.h"
#include "optimize.h"
#include "potential_input.h"
#include "potential_output.h"
#include "random.h"
#include "rescale.h"
#include "utils.h"

#define EPS 0.1
#define NEPS 4
#define NSTEP 20
#define NTEMP (3 * g_calc.ndim)
#define STEPVAR 2.0
#define TEMPVAR 0.85
#define KMAX 1000

#define ONE_OVER_SQRT_2_PI 0.39894228040143267794
#define GAUSS(a) (ONE_OVER_SQRT_2_PI * (exp(-((a) * (a)) / 2.0)))

/****************************************************************
 *
 * void randomize_parameter
 *      const int n:    index of parameter to change
 * 	double *xi: 	pointer to all parameters
 * 	double *v: 	pointer to displacement vector
 *
 * Function to generate random parameters for analytic/tabulated
 * potentials.
 *
 * APOT:
 *      We create up to 10 new random parameters until we find one inside
 *      the predefined range, specified by the user.
 * TAB:
 *      Displaces equidistant sampling points of a function.
 *      Displacement is given by gaussian.
 *
 ****************************************************************/

void randomize_parameter(const int n, double* xi, double* v)
{
#if defined(APOT)
  const double min =
      g_pot.apot_table
          .pmin[g_pot.apot_table.idxpot[n]][g_pot.apot_table.idxparam[n]];
  const double max =
      g_pot.apot_table
          .pmax[g_pot.apot_table.idxpot[n]][g_pot.apot_table.idxparam[n]];
  double temp = 0.0;

  if (v[n] > max - min)
    v[n] = max - min;

  // try 10 times, then just use the last value

  for (int i = 0; i < 10; i++) {
    temp = xi[g_pot.opt_pot.idx[n]] + ((2.0 * eqdist() - 1.0) * v[n]);

    if (temp >= min && temp <= max)
      break;
  }

  xi[g_pot.opt_pot.idx[n]] = temp;
#else
  int pot_index = 0;
  const double width = fabs(normdist());
  const double height = normdist() * v[n];

  /* find potential to which the parameter belongs */
  while (g_pot.opt_pot.last[pot_index] < g_pot.opt_pot.idx[n])
    pot_index++;

  for (int i = 0; i <= 4.0 * width; i++) {
    /* using idx avoids moving fixed points */
    if ((n + i <= g_calc.ndim) &&
        (g_pot.opt_pot.idx[n + i] <= g_pot.opt_pot.last[pot_index])) {
      xi[g_pot.opt_pot.idx[n + i]] += GAUSS((double)i / width) * height;
    }
    if ((n - i >= 0) &&
        (g_pot.opt_pot.idx[n - i] >= g_pot.opt_pot.first[pot_index])) {
      xi[g_pot.opt_pot.idx[n - i]] += GAUSS((double)i / width) * height;
    }
  }
#endif  // APOT
}

/****************************************************************
 *
 * double get_annealing_temperature
 *      const double* xi:
 *      double* xi_new:
 *      double* forces:
 *      const double* displacements:
 *      double F:
 *
 * Function to determine annealing starting temperature.
 *
 ****************************************************************/

double get_annealing_temperature(const double* xi, double* xi_new,
                                 double* forces, double* displacements,
                                 double F)
{
  double T = -1.0;

  /* check for automatic temperature */

  if (tolower(g_param.anneal_temp[0]) == 'a') {
    double chi = 0.8;
    double dF = 0.0;

    int u = 10 * g_calc.ndim;
    int m1 = 0;

    printf("Determining optimal starting temperature T ...\n");

    for (int e = 0; e < u; e++) {
      memcpy(xi_new, xi, g_calc.ndimtot * sizeof(double));

      randomize_parameter((int)(eqdist() * g_calc.ndim), xi_new, displacements);

      double F_new = calc_forces(xi_new, forces, 0);

      if (F_new <= F) {
        m1++;
      } else {
        dF += F_new - F;
      }
    }

    printf("Performed %d trial steps, %d of them were downhill.\n", u, m1);

    u -= m1;
    dF /= u;
    T = dF / log(u / (u * chi + (1 - chi) * m1));

    if (isnan(T) || isinf(T))
      error(1, "Simann failed because T was %f, please set it manually.\n", T);
    if (T < 0)
      T = -T;

    printf("Setting T=%f\n\n", T);
  } else {
    T = atof(g_param.anneal_temp);

    if (T < 0)
      error(1, "The value for anneal_temp (%f) is invalid!\n", T);
  }

  return T;
}

/****************************************************************
 *
 * store_pot_data
 *      pot_data_t* pot_data:
 *
 ****************************************************************/

typedef struct {
  double* begin;
  double* end;
  double* step;
  double* invstep;
  double* xcoord;
} pot_data_t;

// Need to save xcoord of this F potential because we use the
// optimum potential in the future, and the current potential
// could be rescaled differently from the optimum

void store_pot_data(pot_data_t* pot_data)
{
#if !defined(APOT)
  if (pot_data->begin == NULL) {
    pot_data->begin = (double*)Malloc(g_param.ntypes * sizeof(double));
    pot_data->end = (double*)Malloc(g_param.ntypes * sizeof(double));
    pot_data->step = (double*)Malloc(g_param.ntypes * sizeof(double));
    pot_data->invstep = (double*)Malloc(g_param.ntypes * sizeof(double));
    pot_data->xcoord = (double*)Malloc(g_calc.ndimtot * sizeof(double));
  }

  int k = 0;

  for (int i = g_calc.paircol + g_param.ntypes;
       i < g_calc.paircol + 2 * g_param.ntypes; i++) {
    pot_data->begin[k] = g_pot.opt_pot.begin[i];
    pot_data->end[k] = g_pot.opt_pot.end[i];
    pot_data->step[k] = g_pot.opt_pot.step[i];
    pot_data->invstep[k] = g_pot.opt_pot.invstep[i];

    // Loop through each spline knot of F
    for (int j = g_pot.opt_pot.first[i]; j <= g_pot.opt_pot.last[i]; j++)
      pot_data->xcoord[j] = g_pot.opt_pot.xcoord[j];
    k++;
  }
#endif  // !APOT
}

/****************************************************************
 *
 * restore_pot_data
 *      pot_data_t* pot_data:
 *
 ****************************************************************/

// Need to put back xcoord of optimum F potential

void restore_pot_data(const pot_data_t* pot_data)
{
  int k = 0;

  for (int i = g_calc.paircol + g_param.ntypes;
       i < g_calc.paircol + 2 * g_param.ntypes; ++i) {
    g_pot.opt_pot.begin[i] = pot_data->begin[k];
    g_pot.opt_pot.end[i] = pot_data->end[k];
    g_pot.opt_pot.step[i] = pot_data->step[k];
    g_pot.opt_pot.invstep[i] = pot_data->invstep[k];

    // Loop through each spline knot of F
    for (int j = g_pot.opt_pot.first[i]; j <= g_pot.opt_pot.last[i]; ++j)
      g_pot.opt_pot.xcoord[j] = pot_data->xcoord[j];
    ++k;
  }
}

/****************************************************************
 *
 * run_simulated_annealing
 * 	double *xi: 	pointer to all parameters
 *
 * Anneals a vector xi to minimize a function F(xi).
 * Algorithm according to Corana et al.
 *
 ****************************************************************/

void run_simulated_annealing(double* const xi)
{
  int loop_counter = 0;
  int loop_again = 0;

#if defined(RESCALE) && !defined(APOT) && \
    (defined(EAM) || defined(ADP) || defined(MEAM))
  int do_rescale = 1;
#endif  // RESCALE && !APOT && (EAM || ADP || MEAM)

  pot_data_t pot_data;

  pot_data.begin = NULL;
  pot_data.end = NULL;
  pot_data.step = NULL;
  pot_data.invstep = NULL;
  pot_data.xcoord = NULL;

  double F = 0;
  double F_opt = 0;
  double F_new = 0;

  /* backlog of previous F values */
  double* F_old = (double*)Malloc(NEPS * sizeof(double));

  /* displacement vectors */
  double* v = (double*)Malloc(g_calc.ndim * sizeof(double));

  /* optimal value */
  double* xi_opt = (double*)Malloc(g_calc.ndimtot * sizeof(double));
  double* xi_new = (double*)Malloc(g_calc.ndimtot * sizeof(double));

  /* latest force vector */
  double* forces = (double*)Malloc(g_calc.mdim * sizeof(double));

  /* number of accepted changes in dir */
  int* naccept = (int*)Malloc(g_calc.ndim * sizeof(int));

  /* init displacement vector */
  for (int i = 0; i < g_calc.ndim; i++)
    v[i] = 0.1;

  memcpy(xi_new, xi, g_calc.ndimtot * sizeof(double));
  memcpy(xi_opt, xi, g_calc.ndimtot * sizeof(double));

  F = calc_forces(xi, forces, 0);

  F_opt = F;

  /* Temperature */
  double T = get_annealing_temperature(xi, xi_new, forces, v, F);

  /* don't anneal if starttemp equal zero */
  if (T == 0.0)
    return;

  store_pot_data(&pot_data);

  printf("  k\tT        \t  m\tF          \tF_opt\n");
  printf("%3d\t%f\t%3d\t%f\t%f\n", 0, T, 0, F, F_opt);
  fflush(stdout);

  for (int n = 0; n < NEPS; n++)
    F_old[n] = F;

  /* annealing loop */
  do {
    for (int m = 0; m < NTEMP; m++) {
      for (int j = 0; j < NSTEP; j++) {
        for (int h = 0; h < g_calc.ndim; h++) {
          /* Step #1 */
          memcpy(xi_new, xi, g_calc.ndimtot * sizeof(double));

          randomize_parameter(h, xi_new, v);

          F_new = calc_forces(xi_new, forces, 0);

          /* accept new point */
          if (F_new <= F) {
#if defined(APOT)
            xi[g_pot.opt_pot.idx[h]] = xi_new[g_pot.opt_pot.idx[h]];
#else
            memcpy(xi, xi_new, g_calc.ndimtot * sizeof(double));
#endif  // APOT
            F = F_new;

            naccept[h]++;

            if (F_new < F_opt) {
              memcpy(xi_opt, xi_new, g_calc.ndimtot * sizeof(double));

              store_pot_data(&pot_data);

              F_opt = F_new;

              if (*g_files.tempfile != '\0') {
#if defined(APOT)
                update_apot_table(xi);
#endif
                write_pot_table_potfit(g_files.tempfile);
              }
            }
          } else if (eqdist() < (exp((F - F_new) / T))) {
            memcpy(xi, xi_new, g_calc.ndimtot * sizeof(double));
            F = F_new;
            naccept[h]++;
          }
        }  // loop over parameters
      }    // steps per temperature

      /* Step adjustment */
      for (int n = 0; n < g_calc.ndim; n++) {
        if (naccept[n] > (0.6 * NSTEP))
          v[n] *= (1 + STEPVAR * ((double)naccept[n] / NSTEP - 0.6) / 0.4);
        else if (naccept[n] < (0.4 * NSTEP))
          v[n] /= (1 + STEPVAR * (0.4 - (double)naccept[n] / NSTEP) / 0.4);
        naccept[n] = 0;
      }

      printf("%3d\t%f\t%3d\t%f\t%f\n", loop_counter, T, m + 1, F, F_opt);
      fflush(stdout);

      /* End annealing if break flagfile exists */
      if (g_files.flagfile && *g_files.flagfile != '\0') {
        FILE* ff = fopen(g_files.flagfile, "r");
        if (NULL != ff) {
          printf("Annealing terminated in presence of break flagfile \"%s\"!\n",
                 g_files.flagfile);
          printf("Temperature was %f, returning optimum configuration\n", T);

          for (int n = 0; n < g_calc.ndimtot; n++)
            xi[n] = xi_opt[n];

          F = F_opt;
          loop_counter = KMAX + 1;
          fclose(ff);
          remove(g_files.flagfile);
          break;
        }
      }

#if defined(RESCALE) && !defined(APOT) && \
    (defined(EAM) || defined(ADP) || defined(MEAM))
      /* Check for rescaling... every tenth step */
      if (((m + 1) % 10 == 0) && (do_rescale == 1)) {
        /* Was rescaling necessary ? */
        if (rescale(&g_pot.opt_pot, 1.0, 0) != 0.0) {
          /* wake other threads and sync potentials */
          printf("F before rescale = %f\n", F);
          F = calc_forces(xi, forces, 2);
          printf("F after rescale = %f\n", F);
        }
      }
#endif  // RESCALE && !APOT && ( EAM || ADP || MEAM )
    }

    /*Temp adjustment */
    T *= TEMPVAR;
    loop_counter++;

    for (int i = 0; i < NEPS - 1; i++)
      F_old[i] = F_old[i + 1];

    F_old[NEPS - 1] = F;

    loop_again = 0;

    for (int n = 0; n < NEPS - 1; n++) {
      if (fabs(F - F_old[n]) > (EPS * F * 0.01)) {
        loop_again = 1;
        break;
      }
    }

    if (!loop_again && ((F - F_opt) > (EPS * F * 0.01))) {
      for (int n = 0; n < g_calc.ndimtot; n++)
        xi[n] = xi_opt[n];

      F = F_opt;

#if defined(MEAM) && !defined(APOT)
      restore_pot_data(&pot_data);

      /* wake other threads and sync potentials */
      F = calc_forces(xi, forces, 2);

#if defined(RESCALE)
      // Turn off rescaling
      do_rescale = 0;
#endif
#endif  // MEAM && !APOT

      loop_again = 1;
    }
  } while (loop_counter < KMAX && loop_again);

  memcpy(xi, xi_opt, g_calc.ndimtot * sizeof(double));

#if defined(MEAM) && !defined(APOT)
  restore_pot_data(&pot_data);

  // wake other threads and sync potentials
  F = calc_forces(xi, forces, 2);
#endif  // MEAM && !APOT
  printf("Finished annealing, starting powell minimization ...\n");

  if (*g_files.tempfile != '\0') {
#if defined(APOT)
    update_apot_table(xi_opt);
#endif  // APOT
    write_pot_table_potfit(g_files.tempfile);
  }
}

#endif  // !EVO
