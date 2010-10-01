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

#include "potfit.h"
#include "utils.h"

#define D (ndimtot+2)
#define NP 10*D			/* number of total population */

/* boundary values for self-adapting parameters */
#define F_LOWER 0.1		/* lower value for F */
#define F_UPPER 0.9		/* upper value for F */
#define TAU_1 0.1		/* probability for changing F */
#define TAU_2 0.1		/* probability for changing CR */

/****************************************************************
 *
 *  initialize population with random numbers
 *
 ****************************************************************/

void init_population(real **pop, real *xi)
{
  int   i, j;
  real  temp, max, min, val;

  for (i = 0; i < NP; i++) {
    for (j = 0; j < (D - 2); j++)
      pop[i][j] = xi[j];
    pop[i][D - 2] = F_LOWER + dsfmt_genrand_close_open(&dsfmt) * F_UPPER;
    pop[i][D - 1] = dsfmt_genrand_close_open(&dsfmt);
  }
  for (i = 1; i < NP; i++) {
    for (j = 0; j < ndim; j++) {
#ifdef APOT
      val = xi[idx[j]];
      min = apot_table.pmin[apot_table.idxpot[j]][apot_table.idxparam[j]];
      max = apot_table.pmax[apot_table.idxpot[j]][apot_table.idxparam[j]];
#else /* APOT */
      val = xi[idx[j]];
      min = .75 * val;
      max = 1.25 * val;
#endif /* APOT */
      /* scale normal distribution to [-1:1] or less */
      temp = normdist() / 3.;
      if (fabs(temp) > 1)
	temp /= fabs(temp);
      if (temp > 0)
	pop[i][idx[j]] = val + temp * (max - val);
      else
	pop[i][idx[j]] = val + temp * (val - min);
    }
  }
}

/****************************************************************
 *
 *  differential evolution
 *
 ****************************************************************/

void diff_evo(real *xi)
{
  int   a, b, c;		/* store randomly picked numbers */
/*  int 	d, e;			|+ enable this line for more vectors +|*/
  int   count = 0;		/* counter for loops */
  int   i, j, k;		/* counters */
  real  avg = 0.;		/* average sum of squares for all configurations */
  real  crit = 1000.;		/* treshold for stopping criterion */
  real  force = 0.;		/* holds the current sum of squares */
  real  min = 10e10;		/* current minimum for all configurations */
  real  max = 0.;		/* current maximum for all configurations */
  real  temp = 0.;		/* temp storage */
#ifdef APOT
  real  pmin = 0.;		/* lower bound for parameter */
  real  pmax = 0.;		/* upper bound for parameter */
#endif
  real *best;			/* best configuration */
  real *cost;			/* cost values for all configurations */
  real *fxi;			/* force vector */
  real *opt;			/* used for force calculation */
  real *trial;			/* current trial configuration */
  real **x1;			/* current population */
  real **x2;			/* next generation */
  FILE *ff;			/* exit flagfile */

  if (evo_threshold == 0)
    return;

  /* vector for force calculation */
  fxi = vect_real(mdim);

  /* vector with new configuration */
  trial = (real *)malloc(D * sizeof(real));

  /* allocate memory for all configurations */
  x1 = (real **)malloc(NP * sizeof(real *));
  x2 = (real **)malloc(NP * sizeof(real *));
  best = (real *)malloc(NP * sizeof(real));
  cost = (real *)malloc(NP * sizeof(real));
  if (x1 == NULL || x2 == NULL || trial == NULL || cost == NULL
    || best == NULL)
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

  printf("Initializing population ... ");
  fflush(stdout);

  init_population(x1, xi);
  for (i = 0; i < NP; i++) {
    cost[i] = (*calc_forces) (x1[i], fxi, 0);
    if (cost[i] < min) {
      min = cost[i];
      for (j = 0; j < D; j++)
	best[j] = x1[i][j];
    }
    if (cost[i] > max)
      max = cost[i];
  }
  for (i = 0; i < NP; i++)
    avg += cost[i];
  printf("done\n");

#ifdef DEBUG
  printf("Starting Differential Evolution with the following parameters:\n");
  printf("D=%d, NP=%d, CR=%f, F=%f\n", D, NP, CR, F);
#endif /* DEBUG */

  crit = max - min;

  printf("Loops\t\tOptimum\t\tAverage error sum\t\tMax-Min\n");
  printf("%5d\t\t%12f\t%20f\t\t%.2e\n", count, min, avg / (NP), crit);
  fflush(stdout);

  /* main differential evolution loop */
  while (crit >= evo_threshold) {
    max = 0;
    /* randomly create new populations */
    for (i = 0; i < NP; i++) {
      /* generate random numbers */
      do
	a = (int)floor(dsfmt_genrand_close_open(&dsfmt) * NP);
      while (a == i);
      do
	b = (int)floor(dsfmt_genrand_close_open(&dsfmt) * NP);
      while (b == i || b == a);
/*      do*/
/*        c = (int)floor(dsfmt_genrand_close_open(&dsfmt) * NP);*/
/*      while (c == i || c == a || c == b);*/
/*      do*/
/*        d = (int)floor(dsfmt_genrand_close_open(&dsfmt) * NP);*/
/*      while (d == i || d == a || d == b || d == c);*/
/*      do*/
/*        e = (int)floor(dsfmt_genrand_close_open(&dsfmt) * NP);*/
/*      while (e == i || e == a || e == b || e == c || e == d);*/
      j = (int)floor(dsfmt_genrand_close_open(&dsfmt) * ndim);

      /* self-adaptive parameters */
      if (dsfmt_genrand_close_open(&dsfmt) < TAU_1)
	trial[D - 2] = F_LOWER + dsfmt_genrand_close_open(&dsfmt) * F_UPPER;
      else
	trial[D - 2] = x1[i][D - 2];
      if (dsfmt_genrand_close_open(&dsfmt) < TAU_2)
	trial[D - 1] = dsfmt_genrand_close_open(&dsfmt);
      else
	trial[D - 1] = x1[i][D - 1];

      /* create trail vectors with different methods */
      for (k = 1; k <= ndim; k++) {
	if (dsfmt_genrand_close_open(&dsfmt) < trial[D - 1] || k == j) {
	  /* DE/rand/1/exp */
/*          temp = x1[c][idx[j]] + trial[D - 2] * (x1[a][idx[j]] - x1[b][idx[j]]);*/
	  /* DE/best/1/exp */
	  temp =
	    best[idx[j]] + trial[D - 2] * (x1[a][idx[j]] - x1[b][idx[j]]);
	  /* DE/rand/2/exp */
/*          temp = x1[e][j] + trial[D-2] * (x1[a][j] + x1[b][j] - x1[c][j] - x1[d][j]);*/
	  /* DE/best/2/exp */
/*          temp = best[j] + trial[D-2] * (x1[a][j] + x1[b][j] - x1[c][j] - x1[d][j]);*/
	  /* DE/rand-to-best/1/exp */
/*          temp = x1[c][j] + (1 - trial[D-2]) * (best[j] - x1[c][j]) +*/
/*            trial[D-2] * (x1[a][j] - x1[b][j]);*/
	  /* DE/rand-to-best/2/exp */
/*          temp = x1[e][j] + (1 - trial[D-2]) * (best[j] - x1[e][j]) +*/
/*            trial[D-2] * (x1[a][j] + x1[b][j] - x1[c][j] - x1[d][j]);*/
#ifdef APOT
	  pmin =
	    apot_table.pmin[apot_table.idxpot[j]][apot_table.idxparam[j]];
	  pmax =
	    apot_table.pmax[apot_table.idxpot[j]][apot_table.idxparam[j]];
	  if (temp > pmax) {
	    trial[idx[j]] = pmax;
	  } else if (temp < pmin) {
	    trial[idx[j]] = pmin;
	  } else
	    trial[idx[j]] = temp;
#else
	  trial[idx[j]] = temp;
#endif
	} else {
	  trial[idx[j]] = x1[i][idx[j]];
	}
	j = (j + 1) % ndim;
      }

      force = (*calc_forces) (trial, fxi, 0);
      if (force < min) {
	for (j = 0; j < D; j++)
	  best[j] = trial[j];
	if (*tempfile != '\0') {
	  for (j = 0; j < ndim; j++)
#ifdef APOT
	    apot_table.values[apot_table.idxpot[j]][apot_table.idxparam[j]] =
	      trial[idx[j]];
	  write_pot_table(&apot_table, tempfile);
#else
	    xi[idx[j]] = trial[idx[j]];
	  write_pot_table(&opt_pot, tempfile);
#endif /* APOT */
	}
	min = force;
      }
      if (force <= cost[i]) {
	for (j = 0; j < D; j++)
	  x2[i][j] = trial[j];
	cost[i] = force;
	if (force > max)
	  max = force;
      } else {
	for (j = 0; j < D; j++)
	  x2[i][j] = x1[i][j];
	if (cost[i] > max)
	  max = cost[i];
      }
    }
    avg = 0;
    for (i = 0; i < NP; i++)
      avg += cost[i];
    printf("%5d\t\t%12f\t%20f\t\t%.2e\n", count + 1, min, avg / (NP),
      max - min);
    fflush(stdout);
    for (i = 0; i < NP; i++)
      for (j = 0; j < D; j++)
	x1[i][j] = x2[i][j];
    count++;

    /* End optimization if break flagfile exists */
    if (*flagfile != '\0') {
      ff = fopen(flagfile, "r");
      if (NULL != ff) {
	printf("\nEvolutionary algorithm terminated ");
	printf("in presence of break flagfile \"%s\"!\n\n", flagfile);
	fclose(ff);
	remove(flagfile);
	break;
      }
    }

    crit = max - min;
  }

  for (j = 0; j < ndimtot; j++)
    xi[j] = best[j];

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
