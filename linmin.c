/****************************************************************
 *
 * linmin.c: Finds the minimum of a multivariable function along
 *	a certain direction
 *
 ****************************************************************
 *
 * Copyright 2002-2008 Peter Brommer, Franz G"ahler
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

/*** rewritten for double precision and zero-offset vectors and matrices 	***/
/*** adapted to Powell requrirements (return vector instead of value)... 	***/
/*** adapted to real variables (ITAP standard) by PB, ITAP, 2002-10-24 		***/
/*** by Peter Brommer, ITAP, 2002-10-10 					***/


#include "potfit.h"
#include "utils.h"
#include "bracket.h"
#include "powell_lsq.h"
#define TOL 1.0e-1

real *xicom, *delcom;

real linmin(real xi[], real del[], real fxi1, real *x1, real *x2,
  real *fret1, real *fret2)
/* takes vector del (direction of search), xi (originating point),
   n,m (dimensions),
   x1, x2 (two best locations),
   fret1, fret2 (return vectors) as arguments */
{
  int   j;
  static real *vecu = NULL;	/* Vector of location u */
  real  xx, fx, fb, bx, ax;
  real  fa = fxi1;
  real  xmin;
  real  xmin2;

  xicom = xi;
  delcom = del;
  ax = 0.0;			/*do not change without correcting fa, */
  /*saves 1 fcalc... */
  bx = .1;

  if (vecu == NULL)
    vecu = vect_real(ndimtot);
  for (j = 0; j < ndimtot; j++)
    vecu[j] = xicom[j] + bx * delcom[j];	/*set vecu */
  fb = (*calc_forces) (vecu, fret2, 0);

  bracket(&ax, &xx, &bx, &fa, &fx, &fb, fret1, fret2);

  fx = brent(ax, xx, bx, fx, TOL, &xmin, &xmin2, fret1, fret2);
  for (j = 0; j < ndimtot; j++) {
    del[j] *= xmin;
    xi[j] += del[j];
  }
  *x1 = xmin;
  *x2 = xmin2;
  return fx;
}

#undef TOL
