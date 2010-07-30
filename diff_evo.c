/****************************************************************
 *
 * diff_evo.c: Implementation of the differential evolution
 *	algorithm for global optimization
 *
 ****************************************************************
 *
 * Copyright 2009-2010 Daniel Schopf
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

#ifdef EVO

#include <math.h>
#include "potfit.h"
#include "utils.h"

#ifdef APOT
#define D ndim
#else
#define D ndimtot
#endif /* APOT */

/* parameters for the differential evolution algorithm */
#define NP 25*D			/* number of total population */
#define CR 0.5			/* crossover constant, in [0,1] */
#define F 0.2			/* coupling constant with others, in [0,1] */

#define MAX_LOOPS 1e6		/* max number of loops performed */
#define MAX_UNCHANGED 100	/* abort after number of unchanged steps */

#ifdef APOT

real *calc_vect(real *x)
{
  int   i, j, k = 0, n = 0;
  static real *vect;

  if (vect == NULL) {
    vect = (real *)malloc(ndimtot * sizeof(real));
    for (i = 0; i < ndimtot; i++)
      vect[i] = 0;
  }

  for (i = 0; i < apot_table.number; i++) {
    vect[k++] = 0;
    vect[k++] = 0;
    for (j = 0; j < apot_table.n_par[i]; j++) {
      if (!invar_pot[i] && !apot_table.invar_par[i][j])
	vect[k++] = x[n++];
      else
	vect[k++] = apot_table.values[i][j];
    }
  }

#ifdef PAIR
  if (enable_cp) {
    for (i = 0; i < ntypes; i++)
      vect[k++] = x[n++];
  }
#endif /* PAIR */

  if (have_globals) {
    for (i = 0; i < apot_table.globals; i++) {
      vect[k++] = x[n++];
    }
  }

  return vect;
}

#endif /* APOT */

void init_population(real **pop, real *xi, int size, real scale)
{
  int   i, j;
  real  temp, max, min, val;

  for (i = 0; i < size; i++)
#ifdef APOT
    pop[0][i] = xi[idx[i]];
#else
    pop[0][i] = xi[i];
#endif /* APOT */
  for (i = 1; i < NP; i++) {
    for (j = 0; j < size; j++) {
#ifdef APOT
      val = xi[idx[j]];
      min = apot_table.pmin[apot_table.idxpot[j]][apot_table.idxparam[j]];
      max = apot_table.pmax[apot_table.idxpot[j]][apot_table.idxparam[j]];
#else
      val = xi[j];
      min = .9 * val;
      max = 1.1 * val;
#endif /* APOT */
      /* force normal distribution to [-1:1] or less */
      temp = normdist() / (3 * scale);
      if (temp > 0)
	pop[i][j] = val + temp * (max - val);
      else
	pop[i][j] = val + temp * (val - min);
    }
  }
}

void diff_evo(real *xi)
{
  int   a, b, c, d, e, i, j, k;
  int   count = 0, last_changed = 0, finished = 0, restart = 0;
  real  force, min = 10e10, avg = 0, temp = 0, sum = 0, tmpsum = 0, pmin =
    0, pmax = 0;
  real *cost, *fxi, *trial, *best, *opt;
  real **x1, **x2;

  if (evo_width == 0)
    return;

  /* vector for force calculation */
  fxi = vect_real(mdim);

  /* vector with new configuration */
  trial = (real *)malloc(D * sizeof(real));

  /* all configurations */
  x1 = (real **)malloc(NP * sizeof(real *));
  x2 = (real **)malloc(NP * sizeof(real *));
  cost = (real *)malloc(NP * sizeof(real));
  best = (real *)malloc(NP * sizeof(real));
  if (x1 == NULL || x2 == NULL || trial == NULL || cost == NULL)
    error("Could not allocate memory for population vector!\n");
  for (i = 0; i < NP; i++) {
    x1[i] = (real *)malloc(D * sizeof(real));
    x2[i] = (real *)malloc(D * sizeof(real));
    if (x1[i] == NULL || x2[i] == NULL)
      error("Could not allocate memory for population vector!\n");
    for (j = 0; j < D; j++) {
      x1[i][j] = 0;
      x2[i][j] = 0;
    }
  }

  init_population(x1, xi, D, 1. / evo_width);
  for (i = 0; i < NP; i++) {
#ifdef APOT
    opt = calc_vect(x1[i]);
#else
    opt = x1[i];
#endif /* APOT */
    cost[i] = (*calc_forces) (opt, fxi, 0);
    if (cost[i] < min) {
      min = cost[i];
      for (j = 0; j < D; j++)
	best[j] = x1[i][j];
    }
  }
  avg = 0;
  for (i = 0; i < NP; i++)
    avg += cost[i];

#ifdef DEBUG
  printf("Starting Differential Evolution with the following parameters:\n");
  printf("D=%d, NP=%d, CR=%f, F=%f\n", D, NP, CR, F);
#endif /* DEBUG */

  printf("Loops\t\tOptimum\t\tAverage error sum\n");
  printf("%5d\t\t%f\t%f\n", count, min, avg / (NP));
  fflush(stdout);

  /* main differential evolution loop */
  while (count < MAX_LOOPS && last_changed < MAX_UNCHANGED && !finished
    && restart < 4) {
    sum = 0;
    /* randomly create new populations */
    for (i = 0; i < NP; i++) {
      tmpsum = 0;
      do
	a = (int)(dsfmt_genrand_close_open(&dsfmt) * NP);
      while (a == i);
      do
	b = (int)(dsfmt_genrand_close_open(&dsfmt) * NP);
      while (b == i || b == a);
      do
	c = (int)(dsfmt_genrand_close_open(&dsfmt) * NP);
      while (c == i || c == a || c == b);
      do
	d = (int)(dsfmt_genrand_close_open(&dsfmt) * NP);
      while (d == i || d == a || d == b || d == c);
      do
	e = (int)(dsfmt_genrand_close_open(&dsfmt) * NP);
      while (e == i || e == a || e == b || e == c || e == d);
      j = (int)(dsfmt_genrand_close_open(&dsfmt) * D);
      for (k = 1; k <= D; k++) {
	if (dsfmt_genrand_close_open(&dsfmt) < CR || k == D) {
	  /* DE/rand/1/exp */
/*          temp = x1[c][j] + F * (x1[a][j] - x1[b][j]);*/
	  /* DE/best/1/exp */
/*          temp = best[j] + F * (x1[a][j] - x1[b][j]);*/
	  /* DE/rand/2/exp */
/*          temp = x1[e][j] + F * (x1[a][j] + x1[b][j] - x1[c][j] - x1[d][j]);*/
	  /* DE/best/2/exp */
/*          temp = best[j] + F * (x1[a][j] + x1[b][j] - x1[c][j] - x1[d][j]);*/
	  /* DE/rand-to-best/1/exp */
/*          temp = x1[c][j] + (1 - F) * (best[j] - x1[c][j]) +*/
/*            F * (x1[a][j] - x1[b][j]);*/
	  /* DE/rand-to-best/2/exp */
	  temp = x1[e][j] + (1 - F) * (best[j] - x1[e][j]) +
	    F * (x1[a][j] + x1[b][j] - x1[c][j] - x1[d][j]);
#ifdef APOT
	  pmin =
	    apot_table.pmin[apot_table.idxpot[j]][apot_table.idxparam[j]];
	  pmax =
	    apot_table.pmax[apot_table.idxpot[j]][apot_table.idxparam[j]];
	  if (temp > pmax || temp < pmin) {
	    trial[j] = x1[(int)(dsfmt_genrand_close_open(&dsfmt) * D)][j];
	  } else
	    trial[j] = temp;
#else
	  trial[j] = temp;
#endif /* APOT */
	  tmpsum += fabs(x1[i][j] - temp);
	} else {
	  trial[j] = x1[i][j];
	}
	j = (j + 1) % D;
      }
#ifdef APOT
      opt = calc_vect(trial);
#else
      opt = trial;
#endif /* APOT */
      force = (*calc_forces) (opt, fxi, 0);
      if (force < min) {
	last_changed = 0;
	for (j = 0; j < D; j++)
	  best[j] = trial[j];
	if (*tempfile != '\0') {
	  for (j = 0; j < ndim; j++)
#ifdef APOT
	    apot_table.values[apot_table.idxpot[j]][apot_table.idxparam[j]] =
	      trial[j];
	  write_pot_table(&apot_table, tempfile);
#else
	    xi[j] = trial[j];
	  write_pot_table(&opt_pot, tempfile);
#endif /* APOT */
	}
	min = force;
      }
      if (force <= cost[i]) {
	if (force < cost[i])
	  sum += tmpsum;
	for (j = 0; j < D; j++)
	  x2[i][j] = trial[j];
	cost[i] = force;
      } else
	for (j = 0; j < D; j++)
	  x2[i][j] = x1[i][j];
    }
    avg = 0;
    for (i = 0; i < NP; i++)
      avg += cost[i];
#ifdef APOT
    printf("%5d\t\t%f\t%f\n", count + 1, min, avg / (NP));
#else
    printf("%5d\t\t%f\t%f\n", count + 1, min, avg / (NP));
#endif /* APOT */
    fflush(stdout);
    for (i = 0; i < NP; i++)
      for (j = 0; j < D; j++)
	x1[i][j] = x2[i][j];
    count++;
    last_changed++;
    if (last_changed == MAX_UNCHANGED && restart > 2) {
      printf
	("\nCould not find any improvements in the last %d steps.\n",
	MAX_UNCHANGED);
      printf("Aborting evolution algorithm ...\n\n");
    }
    if ((avg / (NP) - min) < 1e-10) {
      printf("Average cost equals minimum cost, nothing more to improve\n");
      finished = 1;
    }
    if (last_changed == MAX_UNCHANGED && restart < 3) {
      restart++;
      printf("\nCould not find any improvements in the last %d steps.\n",
	MAX_UNCHANGED);
      printf("Restarting algorithm. (%d tries left)\n\n", 3 - restart);
      init_population(x1, xi, D, log(restart * exp(10)) / evo_width);
      for (i = 0; i < NP; i++) {
#ifdef APOT
	opt = calc_vect(x1[i]);
#else
	opt = x1[i];
#endif /* APOT */
	cost[i] = (*calc_forces) (opt, fxi, 0);
	if (cost[i] < min) {
	  min = cost[i];
	  for (j = 0; j < D; j++)
	    best[j] = x1[i][j];
	}
      }
      last_changed = 0;
    }
  }
#ifdef APOT
  opt = calc_vect(best);
#else
  opt = best;
#endif /* APOT */
  for (j = 0; j < ndimtot; j++)
    xi[j] = opt[j];

  /* clean up */
  for (i = 0; i < NP; i++) {
    free(x1[i]);
    free(x2[i]);
  }
  free(x1);
  free(x2);
  free(trial);
  free(cost);
  free(best);
  free(fxi);
}

#endif /* EVO */
