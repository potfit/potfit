/****************************************************************
 *
 * diff_evo.c: Implementation of the differential evolution
 *	algorithm for global optimization
 *
 ****************************************************************
 *
 * Copyright 2002-2018 - the potfit development team
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

#include "potfit.h"

#if defined(EVO)

#include "force.h"
#include "memory.h"
#include "optimize.h"
#include "potential_output.h"
#include "random.h"
#include "rescale.h"
#include "utils.h"

#define D (g_calc.ndimtot + 2)
#define NP 15 * D  // number of total population

#define JR 0.6  // jumping rate for opposite algorithm

// boundary values for self-adapting parameters
#define F_LOWER 0.1 /* lower value for F */
#define F_UPPER 0.9 /* upper value for F */
#define TAU_1 0.1   /* probability for changing F */
#define TAU_2 0.1   /* probability for changing CR */

#if defined(APOT) || defined(KIM)
void opposite_check(double** population, double* cost, int do_init);
void quicksort(double* cost, int start, int end, double** population);
int partition(double* cost, int start, int end, int index, double** population);
void swap_population(double* population_1, double* population_2);
#endif  // APOT

/****************************************************************
 *
 *  initialize population with random numbers
 *
 ****************************************************************/

void init_population(double** pop, double* xi, double* cost)
{
  // copy initial population into all populations
  for (int i = 0; i < NP; i++) {
    for (int j = 0; j < (D - 2); j++)
      pop[i][j] = xi[j];
    pop[i][D - 2] = F_LOWER + eqdist() * F_UPPER;
    pop[i][D - 1] = eqdist();
  }

  // create random populations (except for first one)
  for (int i = 1; i < NP; i++) {
    for (int j = 0; j < g_calc.ndim; j++) {
      double val = xi[g_pot.opt_pot.idx[j]];
#if defined(APOT) || defined(KIM)
      double min = g_pot.apot_table.pmin[g_pot.apot_table.idxpot[j]][g_pot.apot_table.idxparam[j]];
      double max = g_pot.apot_table.pmax[g_pot.apot_table.idxpot[j]][g_pot.apot_table.idxparam[j]];

      // initialize with normal distribution
      double temp = normdist() / 3.0;

      if (fabs(temp) > 1)
        temp /= fabs(temp);

      if (temp > 0)
        pop[i][g_pot.opt_pot.idx[j]] = val + temp * (max - val);
      else
        pop[i][g_pot.opt_pot.idx[j]] = val + temp * (val - min);
#else
      pop[i][g_pot.opt_pot.idx[j]] = 10.0 * val * (2 * eqdist() - 1);
#endif  // APOT || KIM
    }
  }

  double* forces = (double*)malloc(g_calc.mdim * sizeof(double));

  for (int i = 0; i < NP; i++)
    cost[i] = calc_forces(pop[i], forces, 0);

#if defined(APOT) || defined(KIM)
  opposite_check(pop, cost, 1);
#endif  // APOT || KIM

  free(forces);
}

#if defined(APOT) || defined(KIM)

/****************************************************************
 *
 *  create and check opposite population
 *
 ****************************************************************/

void opposite_check(double** population, double* cost, int do_init)
{
  printf("opposite_check size = %u\n", g_calc.mdim);
  double* fxi = (double*)malloc(g_calc.mdim * sizeof(double));
  double min = 0.0;
  double max = 0.0;
  double minp[g_calc.ndim];
  double maxp[g_calc.ndim];

  static double* tot_cost;  // cost of two populations
  static double** tot_P;    // two populations

  // allocate memory if not done yet
  if (tot_P == NULL) {
    tot_P = (double**)Malloc(2 * NP * sizeof(double*));

    for (int i = 0; i < 2 * NP; i++)
      tot_P[i] = (double*)Malloc(D * sizeof(double));
  }

  if (tot_cost == NULL)
    tot_cost = (double*)Malloc(2 * NP * sizeof(double));

  if (!do_init) {
    for (int i = 0; i < g_calc.ndim; i++) {
      minp[i] = 10e30;
      maxp[i] = -10e30;
    }

    for (int i = 0; i < NP; i++) {
      for (int j = 0; j < g_calc.ndim; j++) {
        if (population[i][g_pot.opt_pot.idx[j]] < minp[j])
          minp[j] = population[i][g_pot.opt_pot.idx[j]];
        if (population[i][g_pot.opt_pot.idx[j]] > maxp[j])
          maxp[j] = population[i][g_pot.opt_pot.idx[j]];
      }
    }
  }

  // generate opposite population
  for (int i = 0; i < NP; i++) {
    memcpy(tot_P[i], population[i], D * sizeof(double));
    memcpy(tot_P[NP + i], population[i], D * sizeof(double));
  }

  for (int i = 0; i < NP; i++) {
    for (int j = 0; j < g_calc.ndim; j++) {
      if (do_init) {
        min =
            g_pot.apot_table
                .pmin[g_pot.apot_table.idxpot[j]][g_pot.apot_table.idxparam[j]];
        max =
            g_pot.apot_table
                .pmax[g_pot.apot_table.idxpot[j]][g_pot.apot_table.idxparam[j]];
      } else {
        min = minp[j];
        max = maxp[j];
      }
      tot_P[i + NP][g_pot.opt_pot.idx[j]] =
          min + max - tot_P[i][g_pot.opt_pot.idx[j]];
    }
    tot_P[i + NP][D - 2] = tot_P[i][D - 2];
    tot_P[i + NP][D - 1] = tot_P[i][D - 1];
  }

  // calculate cost of opposite population
  for (int i = 0; i < NP; i++)
    tot_cost[i] = cost[i];

  for (int i = NP; i < 2 * NP; i++)
    tot_cost[i] = calc_forces(tot_P[i], fxi, 0);

  // evaluate the NP best individuals from both populations
  // sort with quicksort and return NP best indivuals

  quicksort(tot_cost, 0, 2 * NP - 1, tot_P);

  for (int i = 0; i < NP; i++) {
    memcpy(population[i], tot_P[i], D * sizeof(double));
    cost[i] = tot_cost[i];
  }

  free(fxi);
}

/****************************************************************
 *
 *  quicksort algorithm for opposition-based diff_evo
 *
 ****************************************************************/

void quicksort(double* cost, int low, int high, double** population)
{
  if (low < high) {
    int index = (low + high) / 2;
    int newIndex = partition(cost, low, high, index, population);
    quicksort(cost, low, newIndex - 1, population);
    quicksort(cost, newIndex + 1, high, population);
  }
}

/****************************************************************
 *
 *  partition
 *
 ****************************************************************/

int partition(double* cost, int low, int high, int index, double** population)
{
  const double ind_val = cost[index];
  double temp = 0.0;

  SWAP(cost[index], cost[high], temp);
  swap_population(population[index], population[high]);

  int current_index = low;

  for (int i = low; i < high; i++) {
    if (cost[i] <= ind_val && i != current_index) {
      SWAP(cost[i], cost[current_index], temp);
      swap_population(population[i], population[current_index]);
      current_index++;
    }
  }

  SWAP(cost[current_index], cost[high], temp);
  swap_population(population[current_index], population[high]);

  return current_index;
}

/****************************************************************
 *
 *  swap_population
 *
 ****************************************************************/

void swap_population(double* pop_1, double* pop_2)
{
  int size = (g_calc.ndimtot + 2) * sizeof(double);
  static double* temp = NULL;

  if (temp == NULL)
    temp = (double*)Malloc((g_calc.ndimtot + 2) * sizeof(double));

  memcpy(temp, pop_1, size);
  memcpy(pop_1, pop_2, size);
  memcpy(pop_2, temp, size);
}

#endif  // APOT || KIM

/****************************************************************
 *
 *  differential evolution
 *
 ****************************************************************/

void run_differential_evolution(double* xi)
{
  int a, b;                /* store randomly picked numbers */
                           //  int 	c; 			/* additional vector */
                           //  int 	d; 			/* additional vector */
                           //  int 	e; 			/* additional vector */
  int count = 0;           /* counter for loops */
  double cost_sum = 0.0;   /* average sum of squares for all configurations */
  double crit = 1000.0;    /* treshold for stopping criterion */
  double min_cost = 10e10; /* current minimum for all configurations */
  double max_cost = 0.0;   /* current maximum for all configurations */
#if defined(APOT) || defined(KIM)
  int jump_steps = 0;
  double jump_rate = JR;
#endif  // APOT

  if (g_param.evo_threshold == 0.0)
    return;

  // vector for force calculation
  double* forces = (double*)Malloc(g_calc.mdim * sizeof(double));

  // vector with new configuration
  double* trial = (double*)Malloc(D * sizeof(double));
  memcpy(trial, xi, (D - 2) * sizeof(double));

  // allocate memory for all configurations
  double** pop_1 = (double**)Malloc(NP * sizeof(double*));
  double** pop_2 = (double**)Malloc(NP * sizeof(double*));
  double* best = (double*)Malloc(D * sizeof(double));
  double* cost = (double*)Malloc(NP * sizeof(double));

  for (int i = 0; i < NP; i++) {
    pop_1[i] = (double*)Malloc(D * sizeof(double));
    pop_2[i] = (double*)Malloc(D * sizeof(double));
  }

  printf("Initializing population ... ");
  fflush(stdout);

  init_population(pop_1, xi, cost);

  for (int i = 0; i < NP; i++) {
    if (cost[i] < min_cost) {
      min_cost = cost[i];
      memcpy(best, pop_1[i], D * sizeof(double));
    }
    if (cost[i] > max_cost)
      max_cost = cost[i];

    cost_sum += cost[i];
  }

  printf("done\n");

  crit = max_cost - min_cost;

  printf("Loops\t\tOptimum\t\tAverage error sum\t\tMax-Min\n");
  printf("%5d\t\t%15f\t%20f\t\t%.2e\n", count, min_cost, cost_sum / (NP), crit);
  fflush(stdout);

  // main differential evolution loop
  while (crit >= g_param.evo_threshold && min_cost >= g_param.evo_threshold) {
    max_cost = 0.0;

    // randomly create new populations
    for (int i = 0; i < NP; i++) {
      // generate random numbers
      do
        a = (int)floor(eqdist() * NP);
      while (a == i);

      do
        b = (int)floor(eqdist() * NP);
      while (b == i || b == a);

      //       do
      //         c = (int)floor(eqdist() * NP);
      //       while (c == i || c == a || c == b);
      //
      //       do
      //         d = (int)floor(eqdist() * NP);
      //       while (d == i || d == a || d == b || d == c);
      //
      //       do
      //         e = (int)floor(eqdist() * NP);
      //       while (e == i || e == a || e == b || e == c || e == d);

      int j = (int)floor(eqdist() * g_calc.ndim);

      // self-adaptive parameters
      if (eqdist() < TAU_1)
        trial[D - 2] = F_LOWER + eqdist() * F_UPPER;
      else
        trial[D - 2] = pop_1[i][D - 2];

      if (eqdist() < TAU_2)
        trial[D - 1] = eqdist();
      else
        trial[D - 1] = pop_1[i][D - 1];

      double temp = 0.0;

      // create trail vectors with different methods
      for (int k = 1; k <= g_calc.ndim; k++) {
        if (eqdist() < trial[D - 1] || k == j) {
          /* DE/rand/1/exp */
          //           temp = pop_1[c][g_pot.opt_pot.idx[j]] + trial[D - 2] *
          //           (pop_1[a][g_pot.opt_pot.idx[j]] -
          //           pop_1[b][g_pot.opt_pot.idx[j]]);
          /* DE/best/1/exp */
          temp = best[g_pot.opt_pot.idx[j]] +
                 trial[D - 2] * (pop_1[a][g_pot.opt_pot.idx[j]] -
                                 pop_1[b][g_pot.opt_pot.idx[j]]);
/* DE/rand/2/exp */
//           temp = pop_1[e][j] + trial[D-2] * (pop_1[a][j] + pop_1[b][j] -
//           pop_1[c][j] -
//           pop_1[d][j]);
/* DE/best/2/exp */
//           temp = best[j] + trial[D-2] * (pop_1[a][j] + pop_1[b][j] -
//           pop_1[c][j] -
//           pop_1[d][j]);
/* DE/rand-to-best/1/exp */
//           temp = pop_1[c][j] + (1 - trial[D-2]) * (best[j] - pop_1[c][j]) +
//           trial[D-2]
//           * (pop_1[a][j] - pop_1[b][j]);
/* DE/rand-to-best/2/exp */
//           temp = pop_1[e][j] + (1 - trial[D-2]) * (best[j] - pop_1[e][j]) +
//           trial[D-2]
//           * (pop_1[a][j] + pop_1[b][j] - pop_1[c][j] - pop_1[d][j]);
#if defined(APOT) || defined(KIM)
          double pmin = g_pot.apot_table.pmin[g_pot.apot_table.idxpot[j]]
                                             [g_pot.apot_table.idxparam[j]];
          double pmax = g_pot.apot_table.pmax[g_pot.apot_table.idxpot[j]]
                                             [g_pot.apot_table.idxparam[j]];

          if (temp > pmax) {
            trial[g_pot.opt_pot.idx[j]] = pmax;
          } else if (temp < pmin) {
            trial[g_pot.opt_pot.idx[j]] = pmin;
          } else
            trial[g_pot.opt_pot.idx[j]] = temp;
#else
          trial[g_pot.opt_pot.idx[j]] = temp;
#endif  // APOT || KIM
        } else {
          trial[g_pot.opt_pot.idx[j]] = pop_1[i][g_pot.opt_pot.idx[j]];
        }

        j = (j + 1) % g_calc.ndim;
      }

      double force = calc_forces(trial, forces, 0);

      if (force < min_cost) {
        memcpy(best, trial, D * sizeof(double));

        if (g_files.tempfile && strlen(g_files.tempfile)) {
          for (j = 0; j < g_calc.ndim; j++)
#if defined(APOT) || defined(KIM)
            g_pot.apot_table.values[g_pot.apot_table.idxpot[j]]
                                   [g_pot.apot_table.idxparam[j]] =
                trial[g_pot.opt_pot.idx[j]];
#else
            xi[g_pot.opt_pot.idx[j]] = trial[g_pot.opt_pot.idx[j]];
#endif  // APOT || KIM
          write_pot_table_potfit(g_files.tempfile);
        }
        min_cost = force;
      }

      if (force <= cost[i]) {
        memcpy(pop_2[i], trial, D * sizeof(double));

        cost[i] = force;

        if (force > max_cost)
          max_cost = force;
      } else {
        memcpy(pop_2[i], pop_1[i], D * sizeof(double));

        if (cost[i] > max_cost)
          max_cost = cost[i];
      }
    }

#if defined(APOT) || defined(KIM)
    if (eqdist() < jump_rate) {
      opposite_check(pop_2, cost, 0);

      jump_steps++;

      if (jump_steps > 10) {
        jump_rate *= 0.9;
        jump_steps = 0;
      }
    }
#endif  // APOT || KIM

    cost_sum = 0.0;

    for (int i = 0; i < NP; i++)
      cost_sum += cost[i];

    printf("%5d\t\t%15f\t%20f\t\t%.2e\n", count + 1, min_cost, cost_sum / (NP),
           max_cost - min_cost);
    fflush(stdout);

    for (int i = 0; i < NP; i++)
      memcpy(pop_1[i], pop_2[i], D * sizeof(double));
    count++;

    /* End optimization if break flagfile exists */
    if (g_files.flagfile && strlen(g_files.flagfile)) {
      FILE* ff = fopen(g_files.flagfile, "r");

      if (NULL != ff) {
        printf("\nEvolutionary algorithm terminated ");
        printf("in presence of break flagfile \"%s\"!\n\n", g_files.flagfile);
        fclose(ff);
        remove(g_files.flagfile);
        break;
      }
    }

    crit = max_cost - min_cost;
  }

  printf("Finished differential evolution.\n");
  fflush(stdout);

  memcpy(xi, best, g_calc.ndimtot * sizeof(double));
}

#endif  // EVO
