/****************************************************************
* 
*  chempot.c: calculation of the chemical potential
*
*****************************************************************/
/*
*   Copyright 2009 Daniel Schopf
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
* $Revision: 1.1 $
* $Date: 2009/03/06 08:52:30 $
*****************************************************************/

#ifdef APOT

#include "potfit.h"

#define SWAP_REAL(x,y,t) t=x;x=y;y=t;

int swap_chem_pot(int i, int j)
{
  real  temp;

  if (i != j) {
    SWAP_REAL(apot_table.values[apot_table.number][i],
	      apot_table.values[apot_table.number][j], temp);
    SWAP_REAL(compnodelist[i - ntypes], compnodelist[j - ntypes], temp);
    SWAP_REAL(apot_table.pmin[apot_table.number][i],
	      apot_table.pmin[apot_table.number][j], temp);
    SWAP_REAL(apot_table.pmax[apot_table.number][i],
	      apot_table.pmax[apot_table.number][j], temp);
    return 0;
  } else
    return -1;
}

int sort_chem_pot_2d()
{
  /* bubble sort */
  int   i, swapped;

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

real chemical_potential_2d(int *n)
{
  real  nfrac, temp;
  real  xl, xr, yl, yr;
  int   i, ntot;

  ntot = n[0] + n[1];
  nfrac = (real)n[1] / ntot;

  if (nfrac == 0 || nfrac == 1 || compnodes == 0)
    return n[0] * apot_table.chempot[0] + n[1] * apot_table.chempot[1];

  i = 0;
  while (nfrac > compnodelist[i] && i < compnodes) {
    i++;
  }

  if (i == 0) {
    xl = 0;
    xr = compnodelist[0];
    yl = apot_table.chempot[0];
    yr = apot_table.chempot[2];
  } else if (i == compnodes) {
    xr = 1;
    xl = compnodelist[compnodes - 1];
    yl = apot_table.chempot[ntypes + compnodes - 1];
    yr = apot_table.chempot[1];
  } else {
    xl = compnodelist[i - 1];
    xr = compnodelist[i];
    yl = apot_table.chempot[ntypes + i - 2];
    yr = apot_table.chempot[ntypes + i - 1];
  }

  temp = (yr - yl) / (xr - xl);
  temp = (yl + (nfrac - xl) * temp);

  return temp * ntot;
}

real chemical_potential_3d(int *n)
{
  return n[0] * apot_table.chempot[0] + n[1] * apot_table.chempot[1] +
    n[2] * apot_table.chempot[2];
}

void init_chemical_potential(int dim)
{
  if (dim == 2)
    sort_chem_pot_2d();
  if (dim == 3)
    printf
      ("Chemical potentials for n>=3 is not implemented.\nFalling back to N_i*mu_i\n");
}

real chemical_potential(int dim, int *n)
{
  if (dim == 1)
    return n[0];
  if (dim == 2)
    return chemical_potential_2d(n);
  if (dim >= 3)
    return chemical_potential_3d(n);
  return 0.;
}

#endif
