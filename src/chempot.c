/****************************************************************
 *
 * chempot.c: calculation of the chemical potential
 *
 ****************************************************************
 *
 * Copyright 2002-2014
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

#if defined APOT && defined PAIR

#include "functions.h"

int swap_chem_pot(int i, int j)
{
  double temp;

  if (i != j) {
    SWAP(g_pot.apot_table.values[g_pot.apot_table.number][i],
         g_pot.apot_table.values[g_pot.apot_table.number][j], temp);
    SWAP(g_pot.compnodelist[i - g_param.ntypes],
         g_pot.compnodelist[j - g_param.ntypes], temp);
    SWAP(g_pot.apot_table.pmin[g_pot.apot_table.number][i],
         g_pot.apot_table.pmin[g_pot.apot_table.number][j], temp);
    SWAP(g_pot.apot_table.pmax[g_pot.apot_table.number][i],
         g_pot.apot_table.pmax[g_pot.apot_table.number][j], temp);
    return 0;
  } else
    return -1;
}

int sort_chem_pot_2d()
{
  /* bubble sort */
  int i, swapped;

  if (g_param.compnodes > 0) do {
      swapped = 0;
      for (i = 0; i < (g_param.compnodes - 1); i++) {
        if (g_pot.compnodelist[i] > g_pot.compnodelist[i + 1]) {
          swap_chem_pot(g_param.ntypes + i, g_param.ntypes + i + 1);
          swapped = 1;
        }
      }
    } while (swapped);

  return 0;
}

double chemical_potential_1d(int *n, double *mu) { return n[0] * mu[0]; }

double chemical_potential_2d(int *n, double *mu)
{
  int i = 0, ntot;
  double nfrac;

  ntot = n[0] + n[1];
  nfrac = (double)n[1] / ntot;

  if (nfrac == 0 || nfrac == 1 || g_param.compnodes == 0) {
    return n[0] * mu[0] + n[1] * mu[1];
  }

  while (nfrac > g_pot.compnodelist[i] && i < g_param.compnodes) {
    i++;
  }

  double xl, xr, yl, yr, temp;

  if (i == 0) {
    xl = 0;
    xr = g_pot.compnodelist[0];
    yl = mu[0];
    yr = mu[2];
  } else if (i == g_param.compnodes) {
    xr = 1;
    xl = g_pot.compnodelist[g_param.compnodes - 1];
    yl = mu[g_param.ntypes + g_param.compnodes - 1];
    yr = mu[1];
  } else {
    xl = g_pot.compnodelist[i - 1];
    xr = g_pot.compnodelist[i];
    yl = mu[g_param.ntypes + i - 2];
    yr = mu[g_param.ntypes + i - 1];
  }

  temp = (yr - yl) / (xr - xl);
  temp = (yl + (nfrac - xl) * temp);

  return temp * ntot;
}

double chemical_potential_3d(int *n, double *mu, int dim)
{
  int i;
  double temp = 0;

  for (i = 0; i < dim; i++) temp += n[i] * mu[i];

  return temp;
}

void init_chemical_potential(int dim)
{
  if (dim == 2) sort_chem_pot_2d();
}

double chemical_potential(int dim, int *n, double *mu)
{
  if (dim == 1) return chemical_potential_1d(n, mu);
  if (dim == 2) return chemical_potential_2d(n, mu);
  if (dim >= 3) return chemical_potential_3d(n, mu, dim);

  return 0.0;
}

#endif /* APOT && PAIR */
