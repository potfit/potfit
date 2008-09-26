/****************************************************************
*
* functions.c: Routines and function calls used for analytic potentials
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

#ifdef APOT

#include "potfit.h"

/*****************************************************************************
*
* return the number of parameters for a specific analytic potential
*
******************************************************************************/

int apot_parameters(char *name)
{
  if (strcmp(name, "lj") == 0) {
    return 2;
  } else if (strcmp(name, "test4") == 0) {
    return 4;
  } else if (strcmp(name, "test6") == 0) {
    return 6;
  }
  return -1;
}

/*****************************************************************************
*
* assign function pointers to corresponding functions
*
******************************************************************************/

int apot_assign_functions(apot_table_t *apt)
{
  int   i;

  for (i = 0; i < apt->number; i++) {
    if (strcmp(apt->names[i], "lj") == 0) {
      apt->fvalue[i] = &lj_value;
    } else if (strcmp(apt->names[i], "test4") == 0) {
      apt->fvalue[i] = &test4_value;
    } else if (strcmp(apt->names[i], "test6") == 0) {
      apt->fvalue[i] = &test6_value;
    } else
      return -1;

  }
  return 0;
}

/*****************************************************************************
*
* check if the given analytic potential is valid
*
******************************************************************************/

void apot_validate_functions(apot_table_t *apt)
{
/* 	 TODO: check if given function is valid */
}

/*****************************************************************************
*
* actual functions representing the analytic potentials
*
******************************************************************************/

void lj_value(real r, real *p, real *f)
{
  real  sig_d_rad6, sig_d_rad12;

  sig_d_rad6 = (p[1] * p[1]) / (r * r);
  sig_d_rad6 = sig_d_rad6 * sig_d_rad6 * sig_d_rad6;
  sig_d_rad12 = sig_d_rad6 * sig_d_rad6;
  /* Lennard-Jones is 4*epsilon*((sigma/r)^12-(sigma/r)^6) */
  *f = 4 * p[0] * (sig_d_rad12 - sig_d_rad6);
}

void test4_value(real r, real *p, real *f)
{
  *f = r + p[0] + p[1] + p[2] + p[3];
}

void test6_value(real r, real *p, real *f)
{
  *f = r + p[0] + p[1] + p[2] + p[3] + p[4] + p[5];
}

void debug_calc_pot(real *pot)
{
  int   i, j;
  for (j = 0; j < apot_table.number; j++)
    for (i = 0; i < APOT_STEPS + 2; i++)
      fprintf(stderr, "%d %f\n", j * (APOT_STEPS + 2) + i,
	      pot[j * (APOT_STEPS + 2) + i]);
  exit(2);
}

void debug_apot(apot_table_t *apt)
{
  int   i, j;

  printf("-----------------------------\nDebug_Apot:\n");
  printf("Found %d analytic potentials in apot_table_t\n", apt->number);
  for (i = 0; i < apt->number; i++) {
    printf("Potential #%d is a %s potential (%f-%f) with %d parameters\n",
	   i, apt->names[i], apt->begin[i], apt->end[i],
	   apot_parameters(apt->names[i]));
    if (smoothen[i])
      printf("Smoothening potential #%d.\n", i + 1);
    if (invar_pot[i])
      printf("potential #%d is fixed and will not be adjusted.\n", i);
    for (j = 0; j < apt->n_par[i]; j++) {
      printf("Parameter #%d is set to %f\n", j + 1, apt->values[i][j]);
    }
  }
  exit(2);
}

#endif /* APOT */
