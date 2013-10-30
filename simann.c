/****************************************************************
 *
 * simann.c: Contains all routines used for simulated annealing.
 *
 *****************************************************************
 *
 * Copyright 2002-2013
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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 *****************************************************************/

#ifndef EVO

#include <ctype.h>

#include "potfit.h"

#include "optimize.h"
#include "potential.h"
#include "utils.h"

#define EPS 0.1
#define NEPS 4
#define NSTEP 20
#define NTEMP (3*ndim)
#define STEPVAR 2.0
#define TEMPVAR 0.85
#define KMAX 1000
#define GAUSS(a) (1.0/sqrt(2*M_PI)*(exp(-((a)*(a))/2.)))

#ifdef APOT

/****************************************************************
 *
 * void randomize_parameter(int n, double *xi, double *v);
 * 	int n: 		index of the parameter
 * 	double *xi: 	pointer to all parameters
 * 	double *v: 	pointer to displacement vector
 *
 * Function to generate random parameters for analytic potentials.
 * We loop over a new random parameter until we find one inside
 * the predefined range, specified by the user.
 *
 ****************************************************************/

void randomize_parameter(int n, double *xi, double *v)
{
  double temp, rand;
  int   done = 0, count = 0;
  double min, max;

  min = apot_table.pmin[apot_table.idxpot[n]][apot_table.idxparam[n]];
  max = apot_table.pmax[apot_table.idxpot[n]][apot_table.idxparam[n]];

  if (v[n] > max - min)
    v[n] = max - min;

  do {
    temp = xi[idx[n]];
    rand = 2.0 * eqdist() - 1.;
    temp += (rand * v[n]);
    if (temp >= min && temp <= max)
      done = 1;
    count++;
  } while (!done);
  xi[idx[n]] = temp;
}

#else

/****************************************************************
 *
 * void makebump(double *x,double width, double height, int center);
 *      double *x: 	pointer to all parameters
 *      double width: 	width of the gaussian
 *      double height: 	height of the gaussian
 *      int center: 	index of the center point
 *
 * Displaces equidistant sampling points of a function.
 * Displacement is given by gaussian of given width and height.
 *
 ****************************************************************/

void makebump(double *x, double width, double height, int center)
{
  int   i, j = 0;

  /* find pot to which center belongs */
  while (opt_pot.last[j] < idx[center])
    j++;
  for (i = 0; i <= 4.0 * width; i++) {
    /* using idx avoids moving fixed points */
    if ((center + i <= ndim) && (idx[center + i] <= opt_pot.last[j])) {
      x[idx[center + i]] += GAUSS((double)i / width) * height;
    }
  }
  for (i = 1; i <= 4.0 * width; i++) {
    if ((center - i >= 0) && (idx[center - i] >= opt_pot.first[j])) {
      x[idx[center - i]] += GAUSS((double)i / width) * height;
    }
  }
  return;
}

#endif /* APOT */

/****************************************************************
 *
 * void anneal(double *x);
 * 	double *x: 	pointer to all parameters
 *
 * Anneals a vector xi to minimize a function F(xi).
 * Algorithm according to Corana et al.
 *
 ****************************************************************/

void anneal(double *xi)
{
  int   h = 0, j = 0, k = 0, n, m = 0;	/* counters */
  int   auto_T = 0;
  int   loopagain;		/* loop flag */
#ifndef APOT
  int   rescaleMe = 1;		/* rescaling flag */
#endif /* APOT */
  double T = -1.;		/* Temperature */
  double F, Fopt, F2;		/* Fn value */
  double *Fvar;			/* backlog of Fn vals */
  double *v;			/* step vector */
  double *xopt, *xi2;		/* optimal value */
  double *fxi1;			/* two latest force vectors */
#ifndef APOT
  double width, height;		/* gaussian bump size */
#endif /* APOT */
  FILE *ff;			/* exit flagfile */
  int  *naccept;		/* number of accepted changes in dir */

  /* check for automatic temperature */
  if (tolower(anneal_temp[0]) == 'a') {
    auto_T = 1;
  } else {
    T = atof(anneal_temp);
    if (T < 0)
      error(1, "The value for anneal_temp (%f) is invalid!\n", T);
  }

  if (T == 0. && auto_T != 1)
    return;			/* don't anneal if starttemp equal zero */

  Fvar = vect_double(KMAX + 5 + NEPS);	/* Backlog of old F values */
  v = vect_double(ndim);
  xopt = vect_double(ndimtot);
  xi2 = vect_double(ndimtot);
  fxi1 = vect_double(mdim);
  naccept = vect_int(ndim);
#ifndef APOT
  // Optimum potential x-coord arrays
  int   col, col2;
  double *optbegin, *optend, *optstep, *optinvstep, *optxcoord;
  optbegin = vect_double(ntypes);
  optend = vect_double(ntypes);
  optstep = vect_double(ntypes);
  optinvstep = vect_double(ntypes);
  optxcoord = vect_double(ndimtot);
#endif /* APOT */

  /* init step vector and optimum vector */
  for (n = 0; n < ndim; n++) {
    v[n] = .1;
    naccept[n] = 0;
  }
  for (n = 0; n < ndimtot; n++) {
    xi2[n] = xi[n];
    xopt[n] = xi[n];
  }
  F = (*calc_forces) (xi, fxi1, 0);
  Fopt = F;
#ifndef APOT
  // Need to save xcoord of this F potential because we use the
  // optimum potential in the future, and the current potential
  // could be rescaled differently from the optimum
  col2 = 0;
  for (col = paircol + ntypes; col < paircol + 2 * ntypes; ++col) {
    optbegin[col2] = opt_pot.begin[col];
    optend[col2] = opt_pot.end[col];
    optstep[col2] = opt_pot.step[col];
    optinvstep[col2] = opt_pot.invstep[col];

    // Loop through each spline knot of F
    for (n = opt_pot.first[col]; n <= opt_pot.last[col]; ++n)
      optxcoord[n] = opt_pot.xcoord[n];
    ++col2;
  }
#endif /* APOT */
  /* determine optimum temperature for annealing */
  if (auto_T) {
    int   e = 0;
    int   u = 10 * ndim;
    int   m1 = 0;
    double dF = 0.;
    double chi = .8;

    printf("Determining optimal starting temperature T ...\n");
    for (e = 0; e < u; e++) {
      for (n = 0; n < ndimtot; n++)
	xi2[n] = xi[n];
      h = (int)(eqdist() * ndim);
#ifdef APOT
      randomize_parameter(h, xi2, v);
#else
      /* Create a gaussian bump,
         width & hight distributed normally */
      width = fabs(normdist());
      height = normdist() * v[h];
      makebump(xi2, width, height, h);
#endif /* APOT */
      F2 = (*calc_forces) (xi2, fxi1, 0);
      if (F2 <= F) {
	m1++;
      } else {
	dF += F2 - F;
      }
    }
    printf("Did %d steps, %d were accepted\n", u, m1);
    u -= m1;
    dF /= u;

    T = dF / log(u / (u * chi + (1 - chi) * m1));
    if (isnan(T) || isinf(T))
      error(1, "Simann failed because T was %f, please set it manually.", T);
    if (T < 0)
      T = -T;
    printf("Setting T=%f\n\n", T);
  }

  printf("  k\tT        \t  m\tF          \tFopt\n");
  printf("%3d\t%f\t%3d\t%f\t%f\n", 0, T, 0, F, Fopt);
  fflush(stdout);
  for (n = 0; n <= NEPS; n++)
    Fvar[n] = F;

  /* annealing loop */
  do {
    for (m = 0; m < NTEMP; m++) {
      for (j = 0; j < NSTEP; j++) {
	for (h = 0; h < ndim; h++) {
	  /* Step #1 */
	  for (n = 0; n < ndimtot; n++) {
	    xi2[n] = xi[n];
	  }
#ifdef APOT
	  randomize_parameter(h, xi2, v);
#else
	  /* Create a gaussian bump,
	     width & hight distributed normally */
	  width = fabs(normdist());
	  height = normdist() * v[h];
	  makebump(xi2, width, height, h);
#endif /* APOT */
	  F2 = (*calc_forces) (xi2, fxi1, 0);
	  if (F2 <= F) {	/* accept new point */
#ifdef APOT
	    xi[idx[h]] = xi2[idx[h]];
#else
	    for (n = 0; n < ndimtot; n++)
	      xi[n] = xi2[n];
#endif /* APOT */
	    F = F2;
	    naccept[h]++;
	    if (F2 < Fopt) {
	      for (n = 0; n < ndimtot; n++)
		xopt[n] = xi2[n];

#ifndef APOT
	      // Need to save xcoord of this F potential because we use the
	      // optimum potential in the future, and the current potential
	      // could be rescaled differently from the optimum
	      col2 = 0;
	      for (col = paircol + ntypes; col < paircol + 2 * ntypes; ++col) {
		optbegin[col2] = opt_pot.begin[col];
		optend[col2] = opt_pot.end[col];
		optstep[col2] = opt_pot.step[col];
		optinvstep[col2] = opt_pot.invstep[col];

		// Loop through each spline knot of F
		for (n = opt_pot.first[col]; n <= opt_pot.last[col]; ++n)
		  optxcoord[n] = opt_pot.xcoord[n];

		++col2;
	      }
#endif /* APOT */
	      Fopt = F2;
	      if (*tempfile != '\0') {
#ifndef APOT
		write_pot_table(&opt_pot, tempfile);
#else
		update_apot_table(xi);
		write_pot_table(&apot_table, tempfile);
#endif /* APOT */
	      }
	    }
	  } else if (eqdist() < (exp((F - F2) / T))) {
	    for (n = 0; n < ndimtot; n++)
	      xi[n] = xi2[n];
	    F = F2;
	    naccept[h]++;
	  }
	}
      }

      /* Step adjustment */
      for (n = 0; n < ndim; n++) {
	if (naccept[n] > (0.6 * NSTEP))
	  v[n] *= (1 + STEPVAR * ((double)naccept[n] / NSTEP - 0.6) / 0.4);
	else if (naccept[n] < (0.4 * NSTEP))
	  v[n] /= (1 + STEPVAR * (0.4 - (double)naccept[n] / NSTEP) / 0.4);
	naccept[n] = 0;
      }

      printf("%3d\t%f\t%3d\t%f\t%f\n", k, T, m + 1, F, Fopt);
      fflush(stdout);

      /* End annealing if break flagfile exists */
      if (*flagfile != '\0') {
	ff = fopen(flagfile, "r");
	if (NULL != ff) {
	  printf("Annealing terminated in presence of break flagfile \"%s\"!\n", flagfile);
	  printf("Temperature was %f, returning optimum configuration\n", T);
	  for (n = 0; n < ndimtot; n++)
	    xi[n] = xopt[n];
	  F = Fopt;
	  k = KMAX + 1;
	  fclose(ff);
	  remove(flagfile);
	  break;
	}
      }
#if !defined APOT && ( defined EAM || defined ADP || defined MEAM ) && !defined NORESCALE
      /* Check for rescaling... every tenth step */
      if (((m + 1) % 10 == 0) && (rescaleMe == 1)) {
	/* Was rescaling necessary ? */
	if (rescale(&opt_pot, 1., 0) != 0.) {
	  /* wake other threads and sync potentials */
	  printf("F before rescale = %f\n", F);
	  F = (*calc_forces) (xi, fxi1, 2);
	  printf("F after rescale = %f\n", F);
	}
      }
#endif /* !APOT && ( EAM || ADP || MEAM ) && !NORESCALE */
    }

    /*Temp adjustment */
    T *= TEMPVAR;
    k++;
    Fvar[k + NEPS] = F;
    loopagain = 0;
    for (n = 1; n <= NEPS; n++) {
      if (fabs(F - Fvar[k - n + NEPS]) > (EPS * F * 0.01))
	loopagain = 1;
    }
    if (!loopagain && ((F - Fopt) > (EPS * F * 0.01))) {
      for (n = 0; n < ndimtot; n++)
	xi[n] = xopt[n];
      F = Fopt;
#if defined MEAM && !defined APOT
      // Need to put back xcoord of optimum F potential
      col2 = 0;
      for (col = paircol + ntypes; col < paircol + 2 * ntypes; ++col) {
	opt_pot.begin[col] = optbegin[col2];
	opt_pot.end[col] = optend[col2];
	opt_pot.step[col] = optstep[col2];
	opt_pot.invstep[col] = optinvstep[col2];

	// Loop through each spline knot of F
	for (n = opt_pot.first[col]; n <= opt_pot.last[col]; ++n)
	  opt_pot.xcoord[n] = optxcoord[n];

	++col2;
      }

      /* wake other threads and sync potentials */
      F = (*calc_forces) (xi, fxi1, 2);

      // Turn off rescaling
      rescaleMe = 0;
#endif /* MEAM && !APOT */

      loopagain = 1;
    }
  } while (k < KMAX && loopagain);
  for (n = 0; n < ndimtot; n++) {
    xi[n] = xopt[n];
  }

#if defined MEAM && !defined APOT
  // Need to put back xcoord of optimum F potential
  col2 = 0;
  for (col = paircol + ntypes; col < paircol + 2 * ntypes; ++col) {
    opt_pot.begin[col] = optbegin[col2];
    opt_pot.end[col] = optend[col2];
    opt_pot.step[col] = optstep[col2];
    opt_pot.invstep[col] = optinvstep[col2];

    // Loop through each spline knot of F
    for (n = opt_pot.first[col]; n <= opt_pot.last[col]; ++n)
      opt_pot.xcoord[n] = optxcoord[n];

    ++col2;
  }

  // wake other threads and sync potentials
  F = (*calc_forces) (xi, fxi1, 2);
#endif /* MEAM && !APOT */
  printf("Finished annealing, starting powell minimization ...\n");

  F = Fopt;
  if (*tempfile != '\0') {
#ifndef APOT
    write_pot_table(&opt_pot, tempfile);
#else
    update_apot_table(xopt);
    write_pot_table(&apot_table, tempfile);
#endif /* APOT */
  }

  free_vect_double(Fvar);
  free_vect_double(v);
  free_vect_double(xopt);
  free_vect_int(naccept);
  free_vect_double(xi2);
  free_vect_double(fxi1);
  return;
}

#endif /* !EVO */
