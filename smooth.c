/****************************************************************
*
*  smooth.c: Calculates a cutoff radius for a smooth cutoff
*
*****************************************************************/
/*
*   Copyright 2008 Daniel Schopf
*             Institute for Theoretical and Applied Physics
*             University of Stuttgart, D-70550 Stuttgart, Germany
*             http://www.itap.physik.uni-stuttgart.de/
*
*****************************************************************/
/*
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
*   along with potfit; if not, write to the Free Software
*   Foundation, Inc., 51 Franklin St, Fifth Floor,
*   Boston, MA  02110-1301  USA
*/
/****************************************************************
* $Revision: 1.7 $
* $Date: 2009/05/13 10:11:19 $
*****************************************************************/

#ifdef APOT

#include "potfit.h"

/******************************************************************************
*
* function for root bisecting used by the smooth algorithm
*
******************************************************************************/

real root_bisect(void (*function) (real, real *, real *), real *p,
		 real x1, real x2, real xacc)
{
  int   i = 0;
  real  f0, fl, fr, x0, xl, xr;

  xl = x1;
  xr = x2;

  function(xl, p, &fl);
  function(xr, p, &fr);
  x0 = xl - fl * (xr - xl) / (fr - fl);
  function(x0, p, &f0);

  while (f0 > xacc && i++ < 50) {
    if (fl * f0 < 0)
      xr = x0;
    else if (f0 * fr < 0)
      xl = x0;
    function(xl, p, &fl);
    function(xr, p, &fr);
    x0 = xl - fl * (xr - xl) / (fr - fl);
    function(x0, p, &f0);
  }
  return x0;
}

/******************************************************************************
*
* function to bracket root used by the smooth algorithm
*
******************************************************************************/

void bracket_root(void (*function) (real, real *, real *), real *p,
		  real x, real *x_lower, real *x_upper)
{
  real  size, f_x, f_temp;
  int   found = 1;

  size = .01 * x;
  function(x, p, &f_x);

  x_lower[0] = 0;
  x_upper[0] = 0;
  x_lower[1] = 0;
  x_upper[1] = 0;

  function(x - size, p, &f_temp);
  while (f_x * f_temp > 0) {
    size *= 1.05;
    function(x - size, p, &f_temp);
    if (size > (x / 2)) {
      found = 0;
      break;
    }
  }
  if (found) {
    x_lower[0] = x - size;
    x_upper[0] = x - size / 1.05;
  }

  found = 1;
  size = .01 * x;
  function(x + size, p, &f_temp);
  while (f_x * f_temp > 0) {
    size *= 1.05;
    function(x + size, p, &f_temp);
    if (size > (x / 2)) {
      found = 0;
      break;
    }
  }
  if (found) {
    x_lower[1] = x + size / 1.05;
    x_upper[1] = x + size;
  }

  return;
}

/******************************************************************************
*
* function for root finding used by the smooth algorithm
*
******************************************************************************/

int find_root(void (*function) (real, real *, real *),
	      real x, real *p, real *root)
{
  int   i = 0;
  real  x_lower[2], x_upper[2], x1, x2;

  bracket_root(function, p, x, x_lower, x_upper);

  if (x_lower[0] != 0 && x_upper[0] != 0) {
    x1 = root_bisect(function, p, x_lower[0], x_upper[0], 1e-6);
    i += 1;
  }
  if (x_lower[1] != 0 && x_upper[1] != 0) {
    x2 = root_bisect(function, p, x_lower[1], x_upper[1], 1e-6);
    i += 2;
  }

  root[0] = 0;
  root[1] = 0;
  if (i == 3) {
    root[0] = x1;
    root[1] = x2;
    return 2;
  } else if (i == 2) {
    root[0] = x2;
    return 1;
  } else if (i == 1) {
    root[0] = x1;
    return 0;
  } else if (i == 0) {
    return 0;
  }
  return -1;
}

/******************************************************************************
*
* fit a parabola to the given potential
*
******************************************************************************/

int fit_parab(void (*function) (real, real *, real *), real x, real *p,
	      real *par)
{
  real  f, g;
  real  temp;

  function(x, p, &f);
  function(x + 1e-6, p, &temp);
  function(x - 1e-6, p, &g);

  g = (temp - g) / 2e-6;

  if (f * g > 0)
    return 0;

  par[1] = g * g / (4 * f);
  par[2] = g - 2 * par[1] * x;
  par[3] = f + par[1] * x * x - g * x;
  par[0] = x - g / (2 * par[1]);

  return 1;
}

/******************************************************************************
*
* function for smooth cutoff radius
*
******************************************************************************/

real smooth(void (*function) (real, real *, real *),
	    real x, real *p, real xmin, real xmax, real *params)
{
  real  dist = 0, f, g, temp;
  real  a, b, c, d, x0 = 0, x1;
  int   i, nroot = 0, done = 0;
  char  msg[255];
  static real *roots;

  if (roots == NULL)
    roots = (real *)malloc(2 * sizeof(real));

  nroot = find_root(function, x, p, roots);

  if (nroot == 2) {
    dist = fabs(roots[1] - roots[0]) / 4;
    if (fabs(roots[1] - x) < fabs(roots[0] - x)) {
      temp = roots[1];
      roots[1] = roots[0];
      roots[0] = temp;
    }
  }

  if (nroot != 0) {
    for (i = 0; i < nroot; i++) {
      if (done)
	break;
      if (roots[i] > xmin && roots[i] < xmax && dist == 0) {
	x0 = roots[i] * 0.95;
      } else if (roots[i] > xmin && roots[i] < xmax && dist != 0) {
	x0 = roots[i] - dist;
      } else {
	x0 = 0;
      }

      if (x0 != 0 && fit_parab(function, x0, p, params)
	  && params[0] < (xmax / 1.2 * CUTOFF_MARGIN) && params[0] > xmin)
	done = 1;
      else
	x0 = 0;
    }
  } else {
    x0 = x * 0.95;
    if (!(x0 != 0 && fit_parab(function, x0, p, params)
	  && params[0] < (xmax / 1.2 * CUTOFF_MARGIN) && params[0] > xmin))
      x0 = 0;
  }

  if (x0 != 0) {
    return x0;
  } else {
    params[0] = x;
    params[1] = 0;
    params[2] = 0;
    params[3] = 0;
    return x;
  }
}

#endif /* APOT */
