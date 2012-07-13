/****************************************************************
 *
 * diff_evo.c: Implementation of the differential evolution
 *	algorithm for global optimization
 *
 ****************************************************************
 *
 * Copyright 2009-2012
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.itap.physik.uni-stuttgart.de/
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

#include "optimize.h"
#include "utils.h"

#define D (ndimtot+2)
#define NP 15*D			/* number of total population */

#define JR 0.6			/* jumping rate for opposite algorithm */

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

void init_population(double **pop, double *xi, double *cost)
{
  int   i, j;
  double temp, max, min, val;
  double fxi[mdim];

  for (i = 0; i < NP; i++) {
    for (j = 0; j < (D - 2); j++)
      pop[i][j] = xi[j];
    pop[i][D - 2] = F_LOWER + eqdist() * F_UPPER;
    pop[i][D - 1] = eqdist();
  }
  for (i = 1; i < NP; i++) {
    for (j = 0; j < ndim; j++) {
      val = xi[idx[j]];
#ifdef APOT
      min = apot_table.pmin[apot_table.idxpot[j]][apot_table.idxparam[j]];
      max = apot_table.pmax[apot_table.idxpot[j]][apot_table.idxparam[j]];
      /* initialize with normal distribution */
      temp = normdist() / 3.;
      if (fabs(temp) > 1)
	temp /= fabs(temp);
      if (temp > 0)
	pop[i][idx[j]] = val + temp * (max - val);
      else
	pop[i][idx[j]] = val + temp * (val - min);
#else /* APOT */
      min = -10. * val;
      max = 10. * val;
      /* initialize with uniform distribution in [-1:1] */
      temp = dsfmt_genrand_close_open(&dsfmt);
/*      pop[i][idx[j]] = temp * 100.;*/
      pop[i][idx[j]] = val + temp * (max - min);
#endif /* APOT */
    }
  }
  for (i = 0; i < NP; i++)
    cost[i] = (*calc_forces) (pop[i], fxi, 0);
#ifdef APOT
  opposite_check(pop, cost, 1);
#endif /* APOT */
}

#ifdef APOT

/****************************************************************
 *
 *  create and check opposite population
 *
 ****************************************************************/

void opposite_check(double **P, double *costP, int init)
{
  int   i, j;
  double fxi[mdim];
  double max, min;
  double minp[ndim], maxp[ndim];
  static double *tot_cost;	/* cost of two populations */
  static double **tot_P;	/* two populations */

  /* allocate memory if not done yet */
  if (tot_P == NULL) {
    tot_P = (double **)malloc(2 * NP * sizeof(double *));
    if (tot_P == NULL)
      error(1, "Could not allocate memory for opposition vector!\n");
    for (i = 0; i < 2 * NP; i++) {
      tot_P[i] = (double *)malloc(D * sizeof(double));
      for (j = 0; j < D; j++)
	tot_P[i][j] = 0.;
    }
  }
  if (tot_cost == NULL)
    tot_cost = (double *)malloc(2 * NP * sizeof(double));
  for (i = 0; i < 2 * NP; i++)
    tot_cost[i] = 0.;

  if (!init) {
    for (i = 0; i < ndim; i++) {
      minp[i] = 10e30;
      maxp[i] = -10e30;
    }
    for (i = 0; i < NP; i++) {
      for (j = 0; j < ndim; j++) {
	if (P[i][idx[j]] < minp[j])
	  minp[j] = P[i][idx[j]];
	if (P[i][idx[j]] > maxp[j])
	  maxp[j] = P[i][idx[j]];
      }
    }
  }

  /* generate opposite population */
  for (i = 0; i < NP; i++)
    for (j = 0; j < D; j++)
      tot_P[i][j] = P[i][j];
  for (i = 0; i < NP; i++) {
    for (j = 0; j < ndim; j++) {
      if (init) {
	min = apot_table.pmin[apot_table.idxpot[j]][apot_table.idxparam[j]];
	max = apot_table.pmax[apot_table.idxpot[j]][apot_table.idxparam[j]];
      } else {
	min = minp[j];
	max = maxp[j];
      }
      tot_P[i + NP][idx[j]] = min + max - tot_P[i][idx[j]];
    }
    tot_P[i + NP][D - 2] = tot_P[i][D - 2];
    tot_P[i + NP][D - 1] = tot_P[i][D - 1];
  }

  /* calculate cost of opposite population */
  for (i = 0; i < NP; i++)
    tot_cost[i] = costP[i];
  for (i = NP; i < 2 * NP; i++)
    tot_cost[i] = (*calc_forces) (tot_P[i], fxi, 0);

  /* evaluate the NP best individuals from both populations */
  /* sort with quicksort and return NP best indivuals */
  quicksort(tot_cost, 0, 2 * NP - 1, tot_P);
  for (i = 0; i < NP; i++) {
    for (j = 0; j < D; j++)
      P[i][j] = tot_P[i][j];
    costP[i] = tot_cost[i];
  }
}

#endif /* APOT */

/****************************************************************
 *
 *  differential evolution
 *
 ****************************************************************/

void diff_evo(double *xi)
{
  int   a, b, c;		/* store randomly picked numbers */
/*  int   d, e; 			|+ enable this line for more vectors +|*/
  int   i, j, k;		/* counters */
  int   count = 0;		/* counter for loops */
  int   jsteps = 0;
  double avg = 0.;		/* average sum of squares for all configurations */
  double crit = 1000.;		/* treshold for stopping criterion */
  double force = 0.;		/* holds the current sum of squares */
  double jumprate = JR;
  double min = 10e10;		/* current minimum for all configurations */
  double max = 0.;		/* current maximum for all configurations */
  double temp = 0.;		/* temp storage */
#ifdef APOT
  double pmin = 0.;		/* lower bound for parameter */
  double pmax = 0.;		/* upper bound for parameter */
#endif /* APOT */
  double *best;			/* best configuration */
  double *cost;			/* cost values for all configurations */
  double *fxi;			/* force vector */
  double *trial;		/* current trial configuration */
  double **x1;			/* current population */
  double **x2;			/* next generation */
  FILE *ff;			/* exit flagfile */

  if (evo_threshold == 0.)
    return;

  /* vector for force calculation */
  fxi = vect_double(mdim);

  /* vector with new configuration */
  trial = (double *)malloc(D * sizeof(double));

  /* allocate memory for all configurations */
  x1 = (double **)malloc(NP * sizeof(double *));
  x2 = (double **)malloc(NP * sizeof(double *));
  best = (double *)malloc(NP * sizeof(double));
  cost = (double *)malloc(NP * sizeof(double));
  if (x1 == NULL || x2 == NULL || trial == NULL || cost == NULL || best == NULL)
    error(1, "Could not allocate memory for population vector!\n");
  for (i = 0; i < NP; i++) {
    x1[i] = (double *)malloc(D * sizeof(double));
    x2[i] = (double *)malloc(D * sizeof(double));
    if (x1[i] == NULL || x2[i] == NULL)
      error(1, "Could not allocate memory for population vector!\n");
    for (j = 0; j < D; j++) {
      x1[i][j] = 0;
      x2[i][j] = 0;
    }
  }

  printf("Initializing population ... ");
  fflush(stdout);

  init_population(x1, xi, cost);
  for (i = 0; i < NP; i++) {
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

  crit = max - min;

  printf("Loops\t\tOptimum\t\tAverage error sum\t\tMax-Min\n");
  printf("%5d\t\t%15f\t%20f\t\t%.2e\n", count, min, avg / (NP), crit);
  fflush(stdout);

  /* main differential evolution loop */
  while (crit >= evo_threshold && min >= evo_threshold) {
    max = 0.;
    /* randomly create new populations */
    for (i = 0; i < NP; i++) {
      /* generate random numbers */
      do
	a = (int)floor(eqdist() * NP);
      while (a == i);
      do
	b = (int)floor(eqdist() * NP);
      while (b == i || b == a);
/*      do*/
/*        c = (int)floor(eqdist() * NP);*/
/*      while (c == i || c == a || c == b);*/
/*      do*/
/*        d = (int)floor(eqdist() * NP);*/
/*      while (d == i || d == a || d == b || d == c);*/
/*      do*/
/*        e = (int)floor(eqdist() * NP);*/
/*      while (e == i || e == a || e == b || e == c || e == d);*/

      j = (int)floor(eqdist() * ndim);

      /* self-adaptive parameters */
      if (eqdist() < TAU_1)
	trial[D - 2] = F_LOWER + eqdist() * F_UPPER;
      else
	trial[D - 2] = x1[i][D - 2];
      if (eqdist() < TAU_2)
	trial[D - 1] = eqdist();
      else
	trial[D - 1] = x1[i][D - 1];

      /* create trail vectors with different methods */
      for (k = 1; k <= ndim; k++) {
	if (eqdist() < trial[D - 1] || k == j) {
	  /* DE/rand/1/exp */
/*          temp = x1[c][idx[j]] + trial[D - 2] * (x1[a][idx[j]] - x1[b][idx[j]]);*/
	  /* DE/best/1/exp */
	  temp = best[idx[j]] + trial[D - 2] * (x1[a][idx[j]] - x1[b][idx[j]]);
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
	  pmin = apot_table.pmin[apot_table.idxpot[j]][apot_table.idxparam[j]];
	  pmax = apot_table.pmax[apot_table.idxpot[j]][apot_table.idxparam[j]];
	  if (temp > pmax) {
	    trial[idx[j]] = pmax;
	  } else if (temp < pmin) {
	    trial[idx[j]] = pmin;
	  } else
	    trial[idx[j]] = temp;
#else
	  trial[idx[j]] = temp;
#endif /* APOT */
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
	    apot_table.values[apot_table.idxpot[j]][apot_table.idxparam[j]] = trial[idx[j]];
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
#ifdef APOT
    if (eqdist() < jumprate) {
      opposite_check(x2, cost, 0);
      jsteps++;
      if (jsteps > 10) {
	jumprate *= 0.9;
	jsteps = 0;
      }
    }
#endif /* APOT */
    avg = 0.;
    for (i = 0; i < NP; i++)
      avg += cost[i];
    printf("%5d\t\t%15f\t%20f\t\t%.2e\n", count + 1, min, avg / (NP), max - min);
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

  printf("Finished differential evolution.\n");
  fflush(stdout);

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
