/****************************************************************
 *
 * potential_input.c: Routines for reading a potential table
 *
 ****************************************************************
 *
 * Copyright 2002-2013
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.sourceforge.net/
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

#ifdef PARABOLA

/****************************************************************
 *
 *  Evaluate value from parabola through three points.
 *  Extrapolates for all k.
 *
 ****************************************************************/

double parab_ed(pot_table_t *pt, double *xi, int col, double r)
{
  double rr, istep, chi, p0, p1, p2, dv, d2v;
  int   k;

  /* renorm to beginning of table */
  rr = r - pt->begin[col];

  /* indices into potential table */
  istep = pt->invstep[col];
  k = 0;
  chi = rr * istep;
  k = pt->first[col];

  /* intermediate values */
  p0 = xi[k++];
  p1 = xi[k++];
  p2 = xi[k];
  dv = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* return the potential value */
  return p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;
}

/****************************************************************
 *
 *  Evaluate value from parabola through three points.
 *  Extrapolates for all k. Nonequidistant points.
 *
 ****************************************************************/

double parab_ne(pot_table_t *pt, double *xi, int col, double r)
{
  double x0, x1, x2, chi0, chi1, chi2, p0, p1, p2;
  int   k;

  /* renorm to beginning of table */
  /* rr = r - pt->begin[col]; */
  k = pt->first[col];
  x0 = pt->xcoord[k];
  p0 = xi[k++];
  x1 = pt->xcoord[k];
  p1 = xi[k++];
  x2 = pt->xcoord[k];
  p2 = xi[k];

  /* indices into potential table */
  chi0 = (r - x0) / (x2 - x1);
  chi1 = (r - x1) / (x2 - x0);
  chi2 = (r - x2) / (x1 - x0);

  /* intermediate values */
  /* dv  = p1 - p0; */
  /* d2v = p2 - 2 * p1 + p0; */

  /* return the potential value */
  return chi1 * chi2 * p0 - chi0 * chi2 * p1 + chi0 * chi1 * p2;

}

/****************************************************************
 *
 *  Evaluate deritvative from parabola through three points.
 *  Extrapolates for all k.
 *
 ****************************************************************/

double parab_grad_ed(pot_table_t *pt, double *xi, int col, double r)
{
  double rr, istep, chi, p0, p1, p2, dv, d2v;
  int   k;

  /* renorm to beginning of table */
  rr = r - pt->begin[col];

  /* indices into potential table */
  istep = pt->invstep[col];
  k = 0;
  chi = rr * istep;
  k = pt->first[col];

  /* intermediate values */
  p0 = xi[k++];
  p1 = xi[k++];
  p2 = xi[k];
  dv = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* return the derivative */
  return istep * (dv + (chi - 0.5) * d2v);
}

/****************************************************************
 *
 *  Evaluate deritvative from parabola through three points.
 *  Extrapolates for all k.
 *
 ****************************************************************/

double parab_grad_ne(pot_table_t *pt, double *xi, int col, double r)
{
  double h0, h1, h2, x0, x1, x2, chi0, chi1, chi2, p0, p1, p2;
  int   k;

  /* renorm to beginning of table */
  k = pt->first[col];
  x0 = pt->xcoord[k];
  p0 = xi[k++];
  x1 = pt->xcoord[k];
  p1 = xi[k++];
  x2 = pt->xcoord[k];
  p2 = xi[k];

  h0 = x2 - x1;
  h1 = x2 - x0;
  h2 = x1 - x0;

  chi0 = (r - x0) / h0;
  chi1 = (r - x1) / h1;
  chi2 = (r - x2) / h2;

  /* return the potential value */
  return (chi2 / h1 + chi1 / h2) * p0 - (chi0 / h2 + chi2 / h0) * p1 + (chi0 / h1 + chi1 / h0) * p2;

}

/****************************************************************
 *
 *  Evaluate value and deritvative from parabola through three points.
 *  Extrapolates for all k.
 *
 ****************************************************************/

double parab_comb_ed(pot_table_t *pt, double *xi, int col, double r, double *grad)
{
  double rr, istep, chi, p0, p1, p2, dv, d2v;
  int   k;

  /* renorm to beginning of table */
  rr = r - pt->begin[col];

  /* indices into potential table */
  istep = pt->invstep[col];
  k = 0;
  chi = rr * istep;
  k = pt->first[col];

  /* intermediate values */
  p0 = xi[k++];
  p1 = xi[k++];
  p2 = xi[k];
  dv = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* set the derivative */
  *grad = istep * (dv + (chi - 0.5) * d2v);
  /* return the potential value */
  return p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;
}

/***************************************************************
 *
 *  Evaluate value and deritvative from parabola through three points.
 *  Extrapolates for all k.
 *
 ***************************************************************/

double parab_comb_ne(pot_table_t *pt, double *xi, int col, double r, double *grad)
{
  double h0, h1, h2, x0, x1, x2, chi0, chi1, chi2, p0, p1, p2;
  int   k;

  /* renorm to beginning of table */
  k = pt->first[col];
  x0 = pt->xcoord[k];
  p0 = xi[k++];
  x1 = pt->xcoord[k];
  p1 = xi[k++];
  x2 = pt->xcoord[k];
  p2 = xi[k];

  h0 = x2 - x1;
  h1 = x2 - x0;
  h2 = x1 - x0;

  chi0 = (r - x0) / h0;
  chi1 = (r - x1) / h1;
  chi2 = (r - x2) / h2;

  /* return the potential value */
  *grad = (chi2 / h1 + chi1 / h2) * p0 - (chi0 / h2 + chi2 / h0) * p1 + (chi0 / h1 + chi1 / h0) * p2;

  return chi1 * chi2 * p0 - chi0 * chi2 * p1 + chi0 * chi1 * p2;
}

#endif /* PARABOLA */
