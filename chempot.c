/****************************************************************
 *
 * chempot.c: calculation of the chemical potential
 *
 ****************************************************************
 *
 * Copyright 2009-2011
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

#if defined APOT && defined PAIR

#include "potfit.h"

#include "functions.h"

int swap_chem_pot(int i, int j)
{
  double temp;

  if (i != j) {
    SWAP(apot_table.values[apot_table.number][i],
      apot_table.values[apot_table.number][j], temp);
    SWAP(compnodelist[i - ntypes], compnodelist[j - ntypes], temp);
    SWAP(apot_table.pmin[apot_table.number][i],
      apot_table.pmin[apot_table.number][j], temp);
    SWAP(apot_table.pmax[apot_table.number][i],
      apot_table.pmax[apot_table.number][j], temp);
    return 0;
  } else
    return -1;
}

int sort_chem_pot_2d()
{
  /* bubble sort */
  int   i, swapped;

  if (compnodes > 0)
    do {
      swapped = 0;
      for (i = 0; i < (compnodes - 1); i++) {
	if (compnodelist[i] > compnodelist[i + 1]) {
	  swap_chem_pot(ntypes + i, ntypes + i + 1);
	  swapped = 1;
	}
      }
    } while (swapped);

  return 0;
}

double chemical_potential_1d(int *n, double *mu)
{
  return n[0] * mu[0];
}

double chemical_potential_2d(int *n, double *mu)
{
  int   i = 0, ntot;
  double nfrac;

  ntot = n[0] + n[1];
  nfrac = (double)n[1] / ntot;

  if (nfrac == 0 || nfrac == 1 || compnodes == 0) {
    return n[0] * mu[0] + n[1] * mu[1];
  }

  while (nfrac > compnodelist[i] && i < compnodes) {
    i++;
  }

  double xl, xr, yl, yr, temp;

  if (i == 0) {
    xl = 0;
    xr = compnodelist[0];
    yl = mu[0];
    yr = mu[2];
  } else if (i == compnodes) {
    xr = 1;
    xl = compnodelist[compnodes - 1];
    yl = mu[ntypes + compnodes - 1];
    yr = mu[1];
  } else {
    xl = compnodelist[i - 1];
    xr = compnodelist[i];
    yl = mu[ntypes + i - 2];
    yr = mu[ntypes + i - 1];
  }

  temp = (yr - yl) / (xr - xl);
  temp = (yl + (nfrac - xl) * temp);

  return temp * ntot;
}

double chemical_potential_3d(int *n, double *mu, int dim)
{
  int   i;
  double temp = 0;

  for (i = 0; i < dim; i++)
    temp += n[i] * mu[i];
  return temp;
}

void init_chemical_potential(int dim)
{
  if (dim == 2)
    sort_chem_pot_2d();
}

double chemical_potential(int dim, int *n, double *mu)
{
  if (dim == 1)
    return chemical_potential_1d(n, mu);
  if (dim == 2)
    return chemical_potential_2d(n, mu);
  if (dim >= 3)
    return chemical_potential_3d(n, mu, dim);
  return 0.;
}

#endif /* APOT && PAIR */
