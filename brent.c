/****************************************************************
 *
 * brent.c: Minmization of a multivariable function according to
 *	Brent's algorithm.
 *
 ****************************************************************
 *
 * Copyright 1996, 1997, 1998, 1999, 2000
 * 	Brian Gough
 * Copyright 2005-2012
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

#include "potfit.h"

#include "bracket.h"
#include "utils.h"

double brent(double ax, double bx, double cx, double fbx, double tol,
  double *xmin, double *xmin2, double *fxmin, double *fxmin2)
/* take bracket (a,b,c), f(b), tol, pointers to xmin, xmin2, vectors fxmin, fxmin2 */
{
  int   iter, j;
  double t2, tolerance;
  double midpoint;
  double x_left;
  double x_right;
  double z;
  double d = 0.;
  double e = 0.;
  double u, f_u;
  double v;
  double w;
  double f_v;
  double f_w;
  double f_z;
  double w_lower, w_upper;
  double *p_w, *p_z, *p_u, *p_temp;

  static double *vecu = NULL, *fxu = NULL;	/* Vector of location u */

  double p = 0, q = 0, r = 0;
  if (fxu == NULL)
    fxu = vect_double(mdim);
  if (vecu == NULL)
    vecu = vect_double(ndimtot);

  z = bx;
  f_z = fbx;
  x_left = (ax < cx ? ax : cx);
  x_right = (ax > cx ? ax : cx);

  v = w = x_left + CGOLD * (x_right - x_left);
  for (j = 0; j < ndimtot; j++)
    vecu[j] = xicom[j] + v * delcom[j];	/*set vecu */

  p_z = fxmin;
  p_w = fxmin2;
  p_u = fxu;
  f_v = f_w = (*calc_forces) (vecu, p_w, 0);

  for (iter = 1; iter <= ITMAX; iter++) {
    midpoint = 0.5 * (x_left + x_right);
    t2 = 2 * (tolerance = (tol * fabs(z) + ZEPS));
    w_lower = (x_left - z);
    w_upper = (x_right - z);
    if (fabs(z - midpoint) <= t2 - 0.5 * (x_right - x_left)) {
      *xmin = z;
      *xmin2 = w;
      /* Put correct values in pointers */
      for (j = 0; j < mdim; j++) {
	w_lower = p_z[j];	/* temporary storage */
	fxmin2[j] = p_w[j];
	fxmin[j] = w_lower;
      }
      return f_z;
    }
    if (fabs(e) > tolerance) {
      /* fit parabola */

      r = (z - w) * (f_z - f_v);
      q = (z - v) * (f_z - f_w);
      p = (z - v) * q - (z - w) * r;
      q = 2 * (q - r);

      if (q > 0) {
	p = -p;
      } else {
	q = -q;
      }

      r = e;
      e = d;

      if (fabs(p) < fabs(0.5 * q * r) && p > q * w_lower && p < q * w_upper) {

	d = p / q;
	u = z + d;

	if ((u - x_left) < t2 || (x_right - u) < t2) {
	  d = (z < midpoint) ? tolerance : -tolerance;
	}
      } else {
	e = (z < midpoint) ? x_right - z : -(z - x_left);
	d = CGOLD * e;
      }
    } else {
      e = (z < midpoint) ? x_right - z : -(z - x_left);
      d = CGOLD * e;
    }


    if (fabs(d) >= tolerance) {
      u = z + d;
    } else {
      u = z + ((d > 0) ? tolerance : -tolerance);
    }


    for (j = 0; j < ndimtot; j++)
      vecu[j] = xicom[j] + u * delcom[j];	/*set vecu */
    f_u = (*calc_forces) (vecu, p_u, 0);

    if (f_u > f_z) {
      if (u < z) {
	x_left = u;
	/* fertig */
      } else {
	x_right = u;
	/* done */
      }

      if (f_u <= f_w || w == z) {
	v = w;
	f_v = f_w;
	w = u;
	f_w = f_u;
	SWAP(p_w, p_u, p_temp);
	/* done */
      } else if (f_u <= f_v || v == z || v == w) {
	v = u;
	f_v = f_u;
	/* done */
      }
    } else if (f_u <= f_z) {
      if (u < z) {
	x_right = z;
      } else {
	x_left = z;
      }
      SHIFT(v, w, z, u);
      SHIFT(f_v, f_w, f_z, f_u);
      SWAP(p_w, p_z, p_temp);
      SWAP(p_z, p_u, p_temp);
      /* done */
    } else {
      error(1, "Problems in Brent minimization.\n");
    }

  }
  error(1, "Too many iterations in Brent minimization.\n");
  return 0.;
}
