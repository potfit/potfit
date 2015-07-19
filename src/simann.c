/****************************************************************
 *
 * simann.c: Contains all routines used for simulated annealing.
 *
 *****************************************************************
 *
 * Copyright 2002-2014
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.sourceforge.net/
 *
 *****************************************************************
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
//  *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 *****************************************************************/

#include "potfit.h"

#if !defined(EVO)

#include <ctype.h>

#include "forces.h"
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
  const double min = g_pot.apot_table.pmin[g_pot.apot_table.idxpot[n]][g_pot.apot_table.idxparam[n]];
  const double max = g_pot.apot_table.pmax[g_pot.apot_table.idxpot[n]][g_pot.apot_table.idxparam[n]];
  double temp = 0.0;

  if (v[n] > max - min)
    v[n] = max - min;

  // try 10 times, then just use the last value

  for (int i = 0; i < 10; i++)
  {
    temp = xi[g_todo.idx[n]] + ((2.0 * eqdist() - 1.0) * v[n]);

    if (temp >= min && temp <= max)
      break;
  }

  xi[g_todo.idx[n]] = temp;
#else
  int pot_index = 0;
  const double width = fabs(normdist());
  const double height = normdist() * v[n];

  /* find potential to which the parameter belongs */
  while (g_pot.opt_pot.last[pot_index] < g_todo.idx[n])
    pot_index++;

  for (int i = 0; i <= 4.0 * width; i++)
  {
    /* using idx avoids moving fixed points */
    if ((n + i <= g_calc.ndim) && (g_todo.idx[n + i] <= g_pot.opt_pot.last[pot_index]))
    {
      xi[g_todo.idx[n + i]] += GAUSS((double)i / width) * height;
    }
    if ((n - i >= 0) && (g_todo.idx[n - i] >= g_pot.opt_pot.first[pot_index]))
    {
      xi[g_todo.idx[n - i]] += GAUSS((double)i / width) * height;
    }
  }
#endif // APOT
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

double get_annealing_temperature(const double* xi, double* xi_new, double* forces, double* displacements, double F)
{
  double T = -1.0;

  /* check for automatic temperature */

  if (tolower(g_param.anneal_temp[0]) == 'a')
  {
    double chi = 0.8;
    double dF = 0.0;

    int u = 10 * g_calc.ndim;
    int m1 = 0;

    printf("Determining optimal starting temperature T ...\n");

    for (int e = 0; e < u; e++)
    {
      memcpy(xi_new, xi, g_calc.ndimtot * sizeof(double));

      randomize_parameter((int)(eqdist() * g_calc.ndim), xi_new, displacements);

      double F_new = (*g_calc_forces)(xi_new, forces, 0);

      if (F_new <= F)
      {
        m1++;
      }
      else
      {
        dF += F_new - F;
      }
    }

    printf("Performed %d trial steps, %d of them were downhill.\n", u, m1);

    u -= m1;
    dF /= u;
    T = dF / log(u / (u * chi + (1 - chi) * m1));

    if (isnan(T) || isinf(T))
      error(1, "Simann failed because T was %f, please set it manually.", T);
    if (T < 0)
      T = -T;

    printf("Setting T=%f\n\n", T);
  }
  else
  {
    T = atof(g_param.anneal_temp);

    if (T < 0)
      error(1, "The value for anneal_temp (%f) is invalid!\n", T);
  }

  return T;
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

void run_simulated_annealing(double* xi)
{
  int loop_counter = 0;
  int loopagain = 0; /* loop flag */

#if defined(RESCALE) && !defined(APOT) && (defined(EAM) || defined(ADP) || defined(MEAM))
  int do_rescale = 1;      /* rescaling flag */
#endif // RESCALE && !APOT && (EAM || ADP || MEAM)

  double F = 0;
  double F_opt = 0;
  double F_new = 0;

  /* backlog of previous F values */
  double F_old[NEPS];

  /* displacement vectors */
  double* v = (double*)Malloc(g_calc.ndim * sizeof(double));

  /* optimal value */
  double* xi_opt = (double*)Malloc(g_calc.ndimtot * sizeof(double));
  double* xi_new = (double*)Malloc(g_calc.ndimtot * sizeof(double));

  /* latest force vector */
  double* forces = (double*)Malloc(g_calc.mdim * sizeof(double));

  /* number of accepted changes in dir */
  int* naccept = (int*)Malloc(g_calc.ndim * sizeof(int));

#if !defined(APOT)
  // Optimum potential x-coord arrays
  int col, col2;
  double* optbegin, *optend, *optstep, *optinvstep, *optxcoord;
  optbegin = (double*)Malloc(g_param.ntypes * sizeof(double));
  optend = (double*)Malloc(g_param.ntypes * sizeof(double));
  optstep = (double*)Malloc(g_param.ntypes * sizeof(double));
  optinvstep = (double*)Malloc(g_param.ntypes * sizeof(double));
  optxcoord = (double*)Malloc(g_calc.ndimtot * sizeof(double));
#endif // !APOT

  /* init displacement vector */
  for (int i = 0; i < g_calc.ndim; i++)
    v[i] = 0.1;

  memcpy(xi_new, xi, g_calc.ndimtot * sizeof(double));
  memcpy(xi_opt, xi, g_calc.ndimtot * sizeof(double));

  F = (*g_calc_forces)(xi, forces, 0);

  F_opt = F;

  /* Temperature */
  double T = get_annealing_temperature(xi, xi_new, forces, v, F);

  /* don't anneal if starttemp equal zero */
  if (T == 0.0)
    return;

#if !defined(APOT)
  // Need to save xcoord of this F potential because we use the
  // optimum potential in the future, and the current potential
  // could be rescaled differently from the optimum
  col2 = 0;
  for (col = g_calc.paircol + g_param.ntypes; col < g_calc.paircol + 2 * g_param.ntypes;
       ++col)
  {
    optbegin[col2] = g_pot.opt_pot.begin[col];
    optend[col2] = g_pot.opt_pot.end[col];
    optstep[col2] = g_pot.opt_pot.step[col];
    optinvstep[col2] = g_pot.opt_pot.invstep[col];

    // Loop through each spline knot of F
    for (int n = g_pot.opt_pot.first[col]; n <= g_pot.opt_pot.last[col]; ++n)
      optxcoord[n] = g_pot.opt_pot.xcoord[n];
    ++col2;
  }
#endif // !APOT

  printf("  k\tT        \t  m\tF          \tF_opt\n");
  printf("%3d\t%f\t%3d\t%f\t%f\n", 0, T, 0, F, F_opt);
  fflush(stdout);

  for (int n = 0; n < NEPS; n++)
    F_old[n] = F;

  /* annealing loop */
  do
  {
    for (int m = 0; m < NTEMP; m++)
    {
      for (int j = 0; j < NSTEP; j++)
      {
        for (int h = 0; h < g_calc.ndim; h++)
        {
          /* Step #1 */
          memcpy(xi_new, xi, g_calc.ndimtot * sizeof(double));

          randomize_parameter(h, xi_new, v);

          F_new = (*g_calc_forces)(xi_new, forces, 0);

          /* accept new point */
          if (F_new <= F)
          {
#if defined(APOT)
            xi[g_todo.idx[h]] = xi_new[g_todo.idx[h]];
#else
            memcpy(xi, xi_new, g_calc.ndimtot * sizeof(double));
#endif // APOT
            F = F_new;

            naccept[h]++;

            if (F_new < F_opt)
            {
              memcpy(xi_opt, xi_new, g_calc.ndimtot * sizeof(double));

#if !defined(APOT)
              // Need to save xcoord of this F potential because we use the
              // optimum potential in the future, and the current potential
              // could be rescaled differently from the optimum
              col2 = 0;
              for (col = g_calc.paircol + g_param.ntypes;
                   col < g_calc.paircol + 2 * g_param.ntypes; ++col)
              {
                optbegin[col2] = g_pot.opt_pot.begin[col];
                optend[col2] = g_pot.opt_pot.end[col];
                optstep[col2] = g_pot.opt_pot.step[col];
                optinvstep[col2] = g_pot.opt_pot.invstep[col];

                // Loop through each spline knot of F
                for (int n = g_pot.opt_pot.first[col]; n <= g_pot.opt_pot.last[col]; ++n)
                  optxcoord[n] = g_pot.opt_pot.xcoord[n];

                ++col2;
              }
#endif // APOT

              F_opt = F_new;

              if (*g_files.tempfile != '\0')
              {
#if defined(APOT)
                update_apot_table(xi);
#endif
                write_pot_table_potfit(g_files.tempfile);
              }
            }
          }
          else if (eqdist() < (exp((F - F_new) / T)))
          {
            memcpy(xi, xi_new, g_calc.ndimtot * sizeof(double));
            F = F_new;
            naccept[h]++;
          }
        } // loop over parameters
      } // steps per temperature

      /* Step adjustment */
      for (int n = 0; n < g_calc.ndim; n++)
      {
        if (naccept[n] > (0.6 * NSTEP))
          v[n] *= (1 + STEPVAR * ((double)naccept[n] / NSTEP - 0.6) / 0.4);
        else if (naccept[n] < (0.4 * NSTEP))
          v[n] /= (1 + STEPVAR * (0.4 - (double)naccept[n] / NSTEP) / 0.4);
        naccept[n] = 0;
      }

      printf("%3d\t%f\t%3d\t%f\t%f\n", loop_counter, T, m + 1, F, F_opt);
      fflush(stdout);

      /* End annealing if break flagfile exists */
      if (*g_files.flagfile != '\0')
      {
        FILE* ff = fopen(g_files.flagfile, "r");
        if (NULL != ff)
        {
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
#if defined(RESCALE)
#if !defined(APOT) && (defined EAM || defined ADP || defined MEAM)
      /* Check for rescaling... every tenth step */
      if (((m + 1) % 10 == 0) && (do_rescale == 1))
      {
        /* Was rescaling necessary ? */
        if (rescale(&g_pot.opt_pot, 1.0, 0) != 0.0)
        {
          /* wake other threads and sync potentials */
          printf("F before rescale = %f\n", F);
          F = (*g_calc_forces)(xi, forces, 2);
          printf("F after rescale = %f\n", F);
        }
      }
#endif /* !APOT && ( EAM || ADP || MEAM ) */
#endif /* RESCALE */
    }

    /*Temp adjustment */
    T *= TEMPVAR;
    loop_counter++;

    for (int i = 0; i < NEPS - 1; i++)
      F_old[i] = F_old[i + 1];
    F_old[NEPS - 1] = F;

    loopagain = 0;

    for (int n = 0; n < NEPS - 1; n++)
    {
      if (fabs(F - F_old[n]) > (EPS * F * 0.01))
        loopagain = 1;
    }

    if (!loopagain && ((F - F_opt) > (EPS * F * 0.01)))
    {
      for (int n = 0; n < g_calc.ndimtot; n++)
        xi[n] = xi_opt[n];
      F = F_opt;
#if defined MEAM && !defined APOT
      // Need to put back xcoord of optimum F potential
      col2 = 0;
      for (col = g_calc.paircol + g_param.ntypes;
           col < g_calc.paircol + 2 * g_param.ntypes; ++col)
      {
        g_pot.opt_pot.begin[col] = optbegin[col2];
        g_pot.opt_pot.end[col] = optend[col2];
        g_pot.opt_pot.step[col] = optstep[col2];
        g_pot.opt_pot.invstep[col] = optinvstep[col2];

        // Loop through each spline knot of F
        for (int n = g_pot.opt_pot.first[col]; n <= g_pot.opt_pot.last[col]; ++n)
          g_pot.opt_pot.xcoord[n] = optxcoord[n];

        ++col2;
      }

      /* wake other threads and sync potentials */
      F = (*g_calc_forces)(xi, forces, 2);

#if defined(RESCALE)
      // Turn off rescaling
      do_rescale = 0;
#endif
#endif /* MEAM && !APOT */

      loopagain = 1;
    }
  } while (loop_counter < KMAX && loopagain);

  memcpy(xi, xi_opt, g_calc.ndimtot * sizeof(double));

#if defined(MEAM) && !defined(APOT)
  // Need to put back xcoord of optimum F potential
  col2 = 0;
  for (col = g_calc.paircol + g_param.ntypes; col < g_calc.paircol + 2 * g_param.ntypes;
       ++col)
  {
    g_pot.opt_pot.begin[col] = optbegin[col2];
    g_pot.opt_pot.end[col] = optend[col2];
    g_pot.opt_pot.step[col] = optstep[col2];
    g_pot.opt_pot.invstep[col] = optinvstep[col2];

    // Loop through each spline knot of F
    for (int n = g_pot.opt_pot.first[col]; n <= g_pot.opt_pot.last[col]; ++n)
      g_pot.opt_pot.xcoord[n] = optxcoord[n];

    ++col2;
  }

  // wake other threads and sync potentials
  F = (*g_calc_forces)(xi, forces, 2);
#endif /* MEAM && !APOT */
  printf("Finished annealing, starting powell minimization ...\n");

  if (*g_files.tempfile != '\0')
  {
#if defined(APOT)
    update_apot_table(xi_opt);
#endif // APOT
    write_pot_table_potfit(g_files.tempfile);
  }
}

#endif  // !EVO
