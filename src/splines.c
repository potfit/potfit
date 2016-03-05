/****************************************************************
 *
 * spline.c: Contains all routines used for spline interpolation
 *	with equidistant or non-equidistant sampling points.
 *
 ****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
 *
 * http://potfit.sourceforge.net/
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

#include "memory.h"
#include "splines.h"

/****************************************************************
 *
 * spline_ed: initializes second derivatives used for spline interpolation
 *            (equidistant x[i])
 *
 ****************************************************************/

void spline_ed(double xstep, double* y, int n, double yp1, double ypn,
               double* y2)
{
  double qn = 0.0;
  double un = 0.0;
  static double* u = NULL;
  static int nmax = 0;

  if (n > nmax) {
    u = (double*)Realloc(u, (n - 1) * sizeof(double));
    nmax = n;
  }

  if (yp1 > 0.99e30) {
    y2[0] = 0.0;
    u[0] = 0.0;
  } else {
    y2[0] = -0.5;
    u[0] = (3.0 / (xstep)) * ((y[1] - y[0]) / (xstep)-yp1);
  }

  for (int i = 1; i < n - 1; i++) {
    /* sig=(x[i]-x[i-1])/(x[i+1]-x[i-1])=.5; */
    double p = 0.5 * y2[i - 1] + 2.0;
    y2[i] = (-0.5) / p;
    u[i] = (y[i + 1] - y[i]) / xstep - (y[i] - y[i - 1]) / (xstep);
    u[i] = (6.0 * u[i] / (2 * xstep) - 0.5 * u[i - 1]) / p;
  }

  if (ypn > 0.99e30) {
    qn = 0.0;
    un = 0.0;
  } else {
    qn = 0.5;
    un = (3.0 / (xstep)) * (ypn - (y[n - 1] - y[n - 2]) / (xstep));
  }

  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);

  for (int k = n - 2; k >= 0; --k)
    y2[k] = y2[k] * y2[k + 1] + u[k];
}

/****************************************************************
 *
 * splint_ed: interpolates the function with splines
 *            (equidistant x[i])
 *
 ****************************************************************/

double splint_ed(pot_table_t* pt, double* xi, int col, double r)
{
  /* check for distances shorter than minimal distance in table */
  double rr = r - pt->begin[col];

  if (rr < 0)
    error(1, "%f %f %d\nShort distance", r, pt->begin[col], col);

  /* indices into potential table */
  int k = (int)(rr * pt->invstep[col]);
  double b = (rr - k * pt->step[col]) * pt->invstep[col];
  k += pt->first[col];
  double a = 1.0 - b;
  double p1 = xi[k];
  double d21 = pt->d2tab[k++];
  double p2 = xi[k];
  double d22 = pt->d2tab[k];

  return a * p1 + b * p2 +
         ((a * a * a - a) * d21 + (b * b * b - b) * d22) /
             (6.0 * pt->invstep[col] * pt->invstep[col]);
}

/****************************************************************
 *
 * splint_comb_ed: calculates spline interpolation of a function (return value)
 *            and its gradiend (grad), equidistant x[i]
 *
 ****************************************************************/

double splint_comb_ed(pot_table_t* pt, double* xi, int col, double r,
                      double* grad)
{
  /* check for distances shorter than minimal distance in table */
  double rr = r - pt->begin[col];

  if (rr < 0)
    error(1, "short distance! in splint_comb_ed");

  /* indices into potential table */
  int k = (int)(rr * pt->invstep[col]);
  double b = (rr - k * pt->step[col]) * pt->invstep[col];
  k += pt->first[col];
  /* This fixes some problems, but causes a lot more ... */
  /*  if (rr = (pt->end[col] - pt->begin[col])) {*/
  /*    return xi[k];*/
  /*  }*/
  double a = 1.0 - b;
  double p1 = xi[k];
  double d21 = pt->d2tab[k++];
  double p2 = xi[k];
  double d22 = pt->d2tab[k];
  *grad = (p2 - p1) * pt->invstep[col] +
          ((3 * (b * b) - 1) * d22 - (3 * (a * a) - 1) * d21) /
              (6.0 * pt->invstep[col]);

  return a * p1 + b * p2 +
         ((a * a * a - a) * d21 + (b * b * b - b) * d22) /
             (6.0 * pt->invstep[col] * pt->invstep[col]);
}

/****************************************************************
 *
 * splint_grad_ed: calculates the first derivative from spline interpolation
 *            (equidistant x[i])
 *
 ****************************************************************/

double splint_grad_ed(pot_table_t* pt, double* xi, int col, double r)
{
  /* check for distances shorter than minimal distance in table */
  double rr = r - pt->begin[col];

  if (rr < 0)
    error(1, "short distance! in splint_grad_ed");

  /* indices into potential table */
  int k = (int)(rr * pt->invstep[col]);
  double b = (rr - k * pt->step[col]) * pt->invstep[col];
  k += pt->first[col];
  /* Check if we are at the last index */
  if (k >= pt->last[col]) {
    k--;
    b += 1.0;
  }
  double a = 1.0 - b;
  double p1 = xi[k];
  double d21 = pt->d2tab[k++];
  double p2 = xi[k];
  double d22 = pt->d2tab[k];

  return (p2 - p1) * pt->invstep[col] +
         ((3 * (b * b) - 1) * d22 - (3 * (a * a) - 1) * d21) /
             (6.0 * pt->invstep[col]);
}

/****************************************************************
 *
 * splint_dir: interpolates the function with splines
 *            (equidistant AND NON-eq.dist x[i])
 *            with known index position
 *
 ****************************************************************/

double splint_dir(pot_table_t* pt, double* xi, int k, double b, double step)
{
  /* indices into potential table */
  double a = 1.0 - b;
  double p1 = xi[k];
  double d21 = pt->d2tab[k++];
  double p2 = xi[k];
  double d22 = pt->d2tab[k];

  return a * p1 + b * p2 +
         ((a * a * a - a) * d21 + (b * b * b - b) * d22) * (step * step) / 6.0;
}

/****************************************************************
 *
 * splint_comb_dir: calculates spline interpolation of a function
 *            (return value)
 *            and its gradiend (grad), equidistant and non-eqd x[i]
 *            with known index position
 *
 ****************************************************************/

double splint_comb_dir(pot_table_t* pt, double* xi, int k, double b,
                       double step, double* grad)
{
  /* indices into potential table */
  double a = 1.0 - b;
  double p1 = xi[k];
  double d21 = pt->d2tab[k++];
  double p2 = xi[k];
  double d22 = pt->d2tab[k];
  *grad = (p2 - p1) / step +
          ((3 * (b * b) - 1) * d22 - (3 * (a * a) - 1) * d21) * step / 6.0;

  return a * p1 + b * p2 +
         ((a * a * a - a) * d21 + (b * b * b - b) * d22) * (step * step) / 6.0;
}

/****************************************************************
 *
 * splint_grad_dir: calculates the first derivative
 *            from spline interpolation (equidistant and NON-eqd. x[i])
 *            with known index position
 *
 ****************************************************************/

double splint_grad_dir(pot_table_t* pt, double* xi, int k, double b,
                       double step)
{
  /* indices into potential table */
  double a = 1.0 - b;
  double p1 = xi[k];
  double d21 = pt->d2tab[k++];
  double p2 = xi[k];
  double d22 = pt->d2tab[k];

  return (p2 - p1) / step +
         ((3 * (b * b) - 1) * d22 - (3 * (a * a) - 1) * d21) * step / 6.0;
}

/****************************************************************
 *
 * spline_ne  : initializes second derivatives used for spline interpolation
 *            (nonequidistant x[i])
 *
 ****************************************************************/

void spline_ne(double* x, double* y, int n, double yp1, double ypn, double* y2)
{
  double qn = 0.0;
  double un = 0.0;
  static double* u = NULL;
  static int nmax = 0;

  if (n > nmax) {
    u = (double*)Realloc(u, (n - 1) * sizeof(double));
    nmax = n;
  }

  if (yp1 > 0.99e30) {
    y2[0] = 0.0;
    u[0] = 0.0;
  } else {
    y2[0] = -0.5;
    u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }

  for (int i = 1; i < n - 1; i++) {
    double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    double p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
           (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
  }

  if (ypn > 0.99e30) {
    qn = 0.0;
    un = 0.0;
  } else {
    qn = 0.5;
    un = (3.0 / (x[n - 1] - x[n - 2])) *
         (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
  }

  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);

  for (int k = n - 2; k >= 0; k--)
    y2[k] = y2[k] * y2[k + 1] + u[k];
}

/****************************************************************
 *
 * splint_ne: interpolates the function with splines
 *            (nonequidistant x[i])
 *
 ****************************************************************/

double splint_ne(pot_table_t* pt, double* xi, int col, double r)
{
  int klo = pt->first[col];
  int khi = pt->last[col];

  /* Find index by bisection */
  while (khi - klo > 1) {
    int k = (khi + klo) >> 1;
    if (pt->xcoord[k] > r)
      khi = k;
    else
      klo = k;
  }

  double h = pt->xcoord[khi] - pt->xcoord[klo];

  double b = (r - pt->xcoord[klo]) / h;
  double a = (1.0 - b);

  return a * xi[klo] + b * xi[khi] +
         ((a * a * a - a) * pt->d2tab[klo] + (b * b * b - b) * pt->d2tab[khi]) *
             (h * h) / 6.0;
}

/******************************************************************************
 *
 * splint_ne_lin: interpolates the function with splines,
 *                linear extrapolation (nonequidistant x[i])
 *
 *****************************************************************************/

double splint_ne_lin(pot_table_t* pt, double* xi, int col, double r)
{
  int klo = pt->first[col];
  int khi = pt->last[col];

  /* Find index by bisection */
  while (khi - klo > 1) {
    int k = (khi + klo) >> 1;
    if (pt->xcoord[k] > r)
      khi = k;
    else
      klo = k;
  }

  double h = pt->xcoord[khi] - pt->xcoord[klo];

  /*   if (h == 0.0) error("Bad xa input to routine splint"); */
  double b = (r - pt->xcoord[klo]) / h;
  double a = (1.0 - b);

  if (r < pt->begin[col]) {
    b = 0.0;
    a = 1.0;
    double grad = (xi[khi] - xi[klo]) / h +
                  ((3 * (b * b) - 1) * pt->d2tab[khi] -
                   (3 * (a * a) - 1) * pt->d2tab[klo]) *
                      h / 6.0;
    return xi[klo] + grad * (r - pt->xcoord[klo]);
  } else if (r > pt->end[col]) {
    b = 1.0;
    a = 0.0;
    double grad = (xi[khi] - xi[klo]) / h +
                  ((3 * (b * b) - 1) * pt->d2tab[khi] -
                   (3 * (a * a) - 1) * pt->d2tab[klo]) *
                      h / 6.0;
    return xi[khi] + grad * (r - pt->xcoord[khi]);
  }

  return a * xi[klo] + b * xi[khi] +
         ((a * a * a - a) * pt->d2tab[klo] + (b * b * b - b) * pt->d2tab[khi]) *
             (h * h) / 6.0;
}

/****************************************************************
 *
 * splint_comb_ne: calculates spline interpolation of a function (return value)
 *            and its gradiend (grad), non-equidistant x[i]
 *
 ****************************************************************/

double splint_comb_ne(pot_table_t* pt, double* xi, int col, double r,
                      double* grad)
{
  int klo = pt->first[col];
  int khi = pt->last[col];

  /* Find index by bisection */
  while (khi - klo > 1) {
    int k = (khi + klo) >> 1;
    if (pt->xcoord[k] > r)
      khi = k;
    else
      klo = k;
  }

  double h = pt->xcoord[khi] - pt->xcoord[klo];

  double b = (r - pt->xcoord[klo]) / h;
  double a = (1.0 - b);

  *grad = (xi[khi] - xi[klo]) / h +
          ((3 * (b * b) - 1) * pt->d2tab[khi] -
           (3 * (a * a) - 1) * pt->d2tab[klo]) *
              h / 6.0;

  return a * xi[klo] + b * xi[khi] +
         ((a * a * a - a) * pt->d2tab[klo] + (b * b * b - b) * pt->d2tab[khi]) *
             (h * h) / 6.0;
}

/****************************************************************
 *
 * splint_grad_ne: calculates the first derivative from spline interpolation
 *            (nonequidistant x[i])
 *
 ****************************************************************/

double splint_grad_ne(pot_table_t* pt, double* xi, int col, double r)
{
  int klo = pt->first[col];
  int khi = pt->last[col];

  /* Find index by bisection */
  while (khi - klo > 1) {
    int k = (khi + klo) >> 1;
    if (pt->xcoord[k] > r)
      khi = k;
    else
      klo = k;
  }

  double h = pt->xcoord[khi] - pt->xcoord[klo];

  double b = (r - pt->xcoord[klo]) / h;
  double a = (1.0 - b);

  return (xi[khi] - xi[klo]) / h +
         ((3 * (b * b) - 1) * pt->d2tab[khi] -
          (3 * (a * a) - 1) * pt->d2tab[klo]) *
             h / 6.0;
}
