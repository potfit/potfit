/****************************************************************
*
*  smooth.c: Calculates a cutoff radius for a smooth cutoff
*
*****************************************************************/
/*
*   Copyright 2008-2009 Daniel Schopf
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
* $Revision: 1.8 $
* $Date: 2009/05/15 08:58:39 $
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

  *x_lower = 0;
  *x_upper = 0;

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
    *x_lower = x - size;
    *x_upper = x - size / 1.05;
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
  real  x_lower, x_upper, x1;

  bracket_root(function, p, x, &x_lower, &x_upper);

  if (x_lower != 0 && x_upper != 0) {
    x1 = root_bisect(function, p, x_lower, x_upper, 1e-6);
    i = 1;
  }

  *root = 0;
  if (i == 1) {
    *root = x1;
    return 0;
  } else if (i == 0) {
    return 0;
  }
  return -1;
}

/******************************************************************************
*
* fit a cubic polynomial to the given potential
*
******************************************************************************/

int fit_cubic(void (*function) (real, real *, real *), real x1, real x0,
	      real *p, real *par)
{
  real  f, g;
  real  temp;

  function(x1, p, &f);
  function(x1 + 1e-6, p, &temp);
  function(x1 - 1e-6, p, &g);

  g = (temp - g) / 2e-6;

  /* *INDENT-OFF* */
  par[0] = (g * (x1 - x0) - 2 * f) / 
	  (x1 * x1 * x1 - x0 * x0 * x0 + 3 * x0 * x1 * (x0 - x1));
  par[1] = (f - par[0] * (x1 * x1 * x1 - 3 * x0 * x0 * x1 + 2 * x0 * x0 * x0)) / 
	  ((x1 - x0) * (x1 - x0));
  par[2] = -(3 * par[0] * x0 * x0 + 2 * par[1] * x0);
  par[3] = -(par[0] * x0 * x0 * x0 + par[1] * x0 * x0 + par[2] * x0);
  /* *INDENT-ON* */

  return 1;
}

/******************************************************************************
*
* function for smooth cutoff radius
*
******************************************************************************/

real smooth(void (*function) (real, real *, real *), real x, real *p,
	    real *params)
{
  real  dist = 0, f, g, temp;
  real  a, b, c, d, x0 = 0, x1;
  int   i, nroot = 0, done = 0;
  char  msg[255];
  real  root;

  nroot = find_root(function, x, p, &root);

  if (nroot == 1) {
    fit_cubic(function, root * 0.9, x, p, params);
    temp = root * 0.9;
  } else {
    fit_cubic(function, x * 0.85, x, p, params);
    temp = x * 0.85;
  }
  return temp;
}

#endif /* APOT */
