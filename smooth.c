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
* $Revision: 1.3 $
* $Date: 2008/11/03 11:46:21 $
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
  int   j = 0;
  real  dx, f, fmid, xmid, rtb;

  function(x1, p, &f);
  function(x2, p, &fmid);
  rtb = f < 0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
  do {
    j++;
    function(xmid = rtb + (dx *= .5), p, &fmid);
    if (fmid <= 0)
      rtb = xmid;
    if (fabs(dx) < xacc || fmid == 0)
      return rtb;
  }
  while (j < 100);
  return 0;
}

/******************************************************************************
*
* function to bracket root used by the smooth algorithm
*
******************************************************************************/

void bracket_range(void (*function) (real, real *, real *), real *p,
		   real x1, real x2, int n, real *xb1, real *xb2, int *nb)
{
  int   nbb, i;
  real  x, fp, fc, dx;

  nbb = 0;
  dx = (x2 - x1) / n;
  function(x = x1, p, &fp);
  for (i = 0; i < n; i++) {
    if (nbb > 9) {
      printf("Too many roots near cutoff radius found\n");
      exit(2);
    }
    function(x += dx, p, &fc);
    if (fc * fp < 0) {
      xb1[nbb] = x - dx;
      xb2[nbb++] = x;
      if (*nb == (nbb - 1))
	return;
    }
    fp = fc;
  }
  *nb = nbb;
}

/******************************************************************************
*
* function for root finding used by the smooth algorithm
*
******************************************************************************/

void find_root(void (*function) (real, real *, real *),
	       real x, real *p, real *root, real *dist, int *cut)
{
  int   n, i = 0;
  real *xb1, *xb2, x1, x2;
  n = 10;
  xb1 = (real *)malloc(n * sizeof(real));
  xb2 = (real *)malloc(n * sizeof(real));
  bracket_range(function, p, 0.8 * x, 1.2 * x, 200, xb1, xb2, &n);
  if (n > 1)
    while ((xb1[i] < x) && (xb2[i] < x)) {
      i++;
      if (i == n) {
	break;
      }
  } else {
    if (n == 1) {
      *cut = 1;
      *root = root_bisect(function, p, xb1[0], xb2[0], 1e-6);
      return;
    }
    *cut = 1;
    *root = 0;
    return;
  }
  x1 = root_bisect(function, p, xb1[--i], xb2[i], 1e-6);
  x2 = root_bisect(function, p, xb1[++i], xb2[i], 1e-6);
  *dist = fabs(x1 - x2) / 4;
  if (fabs(x1 - x) > fabs(x - x2)) {
    *root = x2;
  } else {
    *root = x1;
  }
  cut = 0;
  return;
}

/******************************************************************************
*
* function for smooth cutoff radius
*
******************************************************************************/

real smooth(void (*function) (real, real *, real *),
	    real x, real *p, real xmin, real xmax, real *params)
{
  real  dist, root, f, g, temp;
  real  a, b, c, x0, x1;
  int   cut = 0;
  char  msg[255];

  find_root(function, x, p, &root, &dist, &cut);
/*   printf("root=%f\n", root); */
  if (cut == 1 || root < xmin || root > xmax) {
    x0 = (root == 0 ? x : root) * 0.98;
  } else {
    x0 = root - dist;
  }

  /* we need the first derivative at x0 */
  function(x0, p, &f);
  function(x0 + 1e-6, p, &temp);
  function(x0 - 1e-6, p, &g);

  g = (temp - g) / 2e-6;
/*   printf("root=%f f=%f g=%f\n", x0, f, g); */

  if (f * g > 0) {
    sprintf(msg, "Could not truncate potential near the cutoff radius.\n");
    error(msg);
  }
  a = g * g / (4 * f);
  b = g - 2 * a * x0;
  c = f + a * x0 * x0 - g * x0;
  x1 = x0 - g / (2 * a);

  params[0] = x1;
  params[1] = a;
  params[2] = b;
  params[3] = c;

  return x0;
}

#endif /* APOT */
