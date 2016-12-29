/****************************************************************
 *
 * chempot.c: calculation of the chemical potential
 *
 ****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
 *
 * https://www.potfit.net/
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

#include "chempot.h"

/****************************************************************
  swap_chem_pot
****************************************************************/

int swap_chem_pot(int idx1, int idx2)
{
  double temp;

  if (idx1 != idx2) {
    SWAP(g_pot.apot_table.values[g_pot.apot_table.number][idx1],
         g_pot.apot_table.values[g_pot.apot_table.number][idx2], temp);
    SWAP(g_pot.compnodelist[idx1 - g_param.ntypes],
         g_pot.compnodelist[idx2 - g_param.ntypes], temp);
    SWAP(g_pot.apot_table.pmin[g_pot.apot_table.number][idx1],
         g_pot.apot_table.pmin[g_pot.apot_table.number][idx2], temp);
    SWAP(g_pot.apot_table.pmax[g_pot.apot_table.number][idx1],
         g_pot.apot_table.pmax[g_pot.apot_table.number][idx2], temp);
    return 0;
  } else
    return -1;
}

/****************************************************************
  sort_chem_pot_2d
****************************************************************/

int sort_chem_pot_2d(void)
{
  int swapped = 0;

  if (g_param.compnodes > 0) {
    do {
      swapped = 0;
      for (int i = 0; i < (g_param.compnodes - 1); i++) {
        if (g_pot.compnodelist[i] > g_pot.compnodelist[i + 1]) {
          swap_chem_pot(g_param.ntypes + i, g_param.ntypes + i + 1);
          swapped = 1;
        }
      }
    } while (swapped);
  }

  return 0;
}

/****************************************************************
  chemical_potential_1d
****************************************************************/

double chemical_potential_1d(int* positions, double* values)
{
  return positions[0] * values[0];
}

/****************************************************************
  chemical_potential_2d
****************************************************************/

double chemical_potential_2d(int* positions, double* values)
{
  const int ntot = positions[0] + positions[1];
  const double nfrac = (double)positions[1] / ntot;

  if (nfrac == 0.0 || nfrac == 1.0 || g_param.compnodes == 0)
    return positions[0] * values[0] + positions[1] * values[1];

  int i = 0;
  while (nfrac > g_pot.compnodelist[i] && i < g_param.compnodes)
    i++;

  double xl = 0.0;
  double xr = 0.0;
  double yl = 0.0;
  double yr = 0.0;

  if (i == 0) {
    xl = 0.0;
    xr = g_pot.compnodelist[0];
    yl = values[0];
    yr = values[2];
  } else if (i == g_param.compnodes) {
    xr = 1.0;
    xl = g_pot.compnodelist[g_param.compnodes - 1];
    yl = values[g_param.ntypes + g_param.compnodes - 1];
    yr = values[1];
  } else {
    xl = g_pot.compnodelist[i - 1];
    xr = g_pot.compnodelist[i];
    yl = values[g_param.ntypes + i - 2];
    yr = values[g_param.ntypes + i - 1];
  }

  return (yl + (nfrac - xl) * (yr - yl) / (xr - xl)) * ntot;
}

/****************************************************************
  chemical_potential_multidim
****************************************************************/

double chemical_potential_multidim(int* positions, double* values, int dim)
{
  double temp = 0.0;

  for (int i = 0; i < dim; i++)
    temp += positions[i] * values[i];

  return temp;
}

/****************************************************************
  init_chemical_potential
****************************************************************/

void init_chemical_potential(int dim)
{
  if (dim == 2)
    sort_chem_pot_2d();
}

/****************************************************************
  chemical_potential
****************************************************************/

double chemical_potential(int dimension, int* positions, double* values)
{
  if (dimension == 1)
    return chemical_potential_1d(positions, values);
  if (dimension == 2)
    return chemical_potential_2d(positions, values);
  if (dimension >= 3)
    return chemical_potential_multidim(positions, values, dimension);

  return 0.0;
}
