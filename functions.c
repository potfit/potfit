/****************************************************************
*
* functions.c: Routines and function calls used for analytic potentials
*
*****************************************************************/
/*
*   Copyright 2008-2009 Daniel Schopf
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
* $Revision: 1.14 $
* $Date: 2009/05/13 10:11:19 $
*****************************************************************/

#ifdef APOT

#include "potfit.h"
#include <mkl_vml.h>

/*****************************************************************************
*
* return the number of parameters for a specific analytic potential
*
******************************************************************************/

int apot_parameters(char *name)
{
  if (strcmp(name, "lj") == 0) {
    return 2;
  } else if (strcmp(name, "eopp") == 0) {
    return 6;
  } else if (strcmp(name, "morse") == 0) {
    return 3;
  } else if (strcmp(name, "softshell") == 0) {
    return 2;
  } else if (strcmp(name, "eopp_exp") == 0) {
    return 6;
  } else if (strcmp(name, "meopp") == 0) {
    return 7;
  } else if (strcmp(name, "power_decay") == 0) {
    return 2;
  } else if (strcmp(name, "pohlong") == 0) {
    return 2;
  }
  else if (strcmp(name, "parabola") == 0) {
    return 3;
  }

  /* template for new potential function called newpot */

  else if (strcmp(name, "newpot") == 0) {
    return 2;
  }

  /* end of template */

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
    } else if (strcmp(apt->names[i], "eopp") == 0) {
      apt->fvalue[i] = &eopp_value;
    } else if (strcmp(apt->names[i], "morse") == 0) {
      apt->fvalue[i] = &morse_value;
    } else if (strcmp(apt->names[i], "softshell") == 0) {
      apt->fvalue[i] = &softshell_value;
    } else if (strcmp(apt->names[i], "eopp_exp") == 0) {
      apt->fvalue[i] = &eopp_exp_value;
    } else if (strcmp(apt->names[i], "meopp") == 0) {
      apt->fvalue[i] = &meopp_value;
    } else if (strcmp(apt->names[i], "power_decay") == 0) {
      apt->fvalue[i] = &power_decay_value;
    } else if (strcmp(apt->names[i], "pohlong") == 0) {
      apt->fvalue[i] = &pohlong_value;
    }
    else if (strcmp(apt->names[i], "parabola") == 0) {
      apt->fvalue[i] = &parabola_value;
    }

/* template for new potential function called newpot */

    else if (strcmp(apt->names[i], "newpot") == 0) {
      apt->fvalue[i] = &newpot_value;
    }

/* end of template */

    else
      return -1;
  }
  return 0;
}

/*****************************************************************************
*
* actual functions representing the analytic potentials
*
******************************************************************************/

/******************************************************************************
*
* lennard-jones potential
*
******************************************************************************/

void lj_value(real r, real *p, real *f)
{
  real  sig_d_rad6, sig_d_rad12;

  sig_d_rad6 = (p[1] * p[1]) / (r * r);
  sig_d_rad6 = sig_d_rad6 * sig_d_rad6 * sig_d_rad6;
  sig_d_rad12 = sig_d_rad6 * sig_d_rad6;

  *f = 4 * p[0] * (sig_d_rad12 - sig_d_rad6);
}

/******************************************************************************
*
* empirical oscillating pair potential
*
******************************************************************************/

void eopp_value(real r, real *p, real *f)
{
  static real x[2], y[2], power[2];

  x[0] = r;
  x[1] = r;
  y[0] = p[1];
  y[1] = p[3];
  vdPow(2, x, y, power);

  *f = p[0] / power[0] + (p[2] / power[1]) * cos(p[4] * r + p[5]);
}

/******************************************************************************
*
* morse potential
*
******************************************************************************/

void morse_value(real r, real *p, real *f)
{
  *f = p[0] * (exp(-2 * p[1] * (r - p[2])) - 2 * exp(-p[1] * (r - p[2])));
}

/******************************************************************************
*
* softshell potential
*
******************************************************************************/

void softshell_value(real r, real *p, real *f)
{
  static real x, y;

  x = p[0] / r;
  y = p[1];
  vdPow(1, &x, &y, f);
}

/******************************************************************************
*
* eopp_exp potential
*
******************************************************************************/

void eopp_exp_value(real r, real *p, real *f)
{
  static real x, y, power;

  x = r;
  y = p[3];
  vdPow(2, &x, &y, &power);

  *f = p[0] * exp(-p[1] * r) + (p[2] / power) * cos(p[4] * r + p[5]);
}

/******************************************************************************
*
* meopp potential
*
******************************************************************************/

void meopp_value(real r, real *p, real *f)
{
  static real x[2], y[2], power[2];

  x[0] = r - p[6];
  x[1] = r;
  y[0] = p[1];
  y[1] = p[3];
  vdPow(2, x, y, power);

  *f = p[0] / power[0] + (p[2] / power[1]) * cos(p[4] * r + p[5]);
}

/******************************************************************************
*
* power_decay potential
*
******************************************************************************/

void power_decay_value(real r, real *p, real *f)
{
  static real x, y, power;

  x = 1. / r;
  y = p[1];
  vdPow(1, &x, &y, &power);

  *f = p[0] * power;
}

/******************************************************************************
*
* pohlong potential
*
******************************************************************************/

void pohlong_value(real r, real *p, real *f)
{
  real  power;

  vdPow(1, &r, &p[1], &power);

  *f = p[0] * (1 - p[1] * log(r)) * power;
}

/******************************************************************************
*
* parabola potential
*
******************************************************************************/

void parabola_value(real r, real *p, real *f)
{
  *f = r*r * p[0] + r*p[1]+p[2];
}

/******************************************************************************
* 
* template for new potential function called mypotential 
* for further information plase have a look at the online documentation
* 
* http://www.itap.physik.uni-stuttgart.de/~imd/potfit/potfit.html
*
******************************************************************************/

/* template for new function */

/******************************************************************************
*
* newpot potential
*
******************************************************************************/

void newpot_value(real r, real *p, real *f)
{
  *f = r * p[0] + p[1];
}

/* end of template */

/******************************************************************************
*
* end of analytic potentials
*
******************************************************************************/

/*****************************************************************************
*
* check if the given analytic potential is valid
*
******************************************************************************/

int apot_validate(int param_index, real new_val)
{
  int   pot_index = apot_table.idxpot[param_index];
  real  x;

  if (pot_index < apot_table.number) {

    /* check if potential vanishes at 3*cutoff */
    apot_table.fvalue[pot_index] (3 * apot_table.end[pot_index],
				  apot_table.values[pot_index], &x);

    if (fabs(x) > 10e-2) {
#ifdef DEBUG
      printf("Validate failed\n");
#endif
      return 0;
    }
  }
  return 1;
}

/*****************************************************************************
*
*  debug function to print the potentials
*
******************************************************************************/

#ifdef DEBUG

void debug_apot()
{
  int   i, j;

  fflush(stdout);
  fprintf(stderr, "\n\n##############################################\n");
  fprintf(stderr, "###########      DEBUG OUTPUT      ###########\n");
  fprintf(stderr, "##############################################\n");
  for (i = 0; i < apot_table.number; i++) {
    fprintf(stderr, "\npotential #%d (type=%s, smooth=%d)\n", i + 1,
	    apot_table.names[i], smooth_pot[i]);
    for (j = 0; j < apot_table.n_par[i]; j++) {
      fprintf(stderr, "parameter %d: name=%s value=%f min=%f max=%f\n", j + 1,
	      apot_table.param_name[i][j], apot_table.values[i][j],
	      apot_table.pmin[i][j], apot_table.pmax[i][j]);
    }
  }
  if (disable_cp) {
    fprintf(stderr, "\nchemical potentials are DISABLED!\n");
  } else {
    fprintf(stderr, "\nchemical potentials:\n");
    for (i = 0; i < ntypes; i++)
      fprintf(stderr, "cp_%d=%f min=%f max=%f\n", i, apot_table.chempot[i],
	      apot_table.pmin[apot_table.number][i],
	      apot_table.pmax[apot_table.number][i]);
    if (compnodes > 0) {
      if (ntypes == 2) {
	fprintf(stderr, "composition nodes:\n");
	for (j = 0; j < compnodes; j++)
	  fprintf(stderr, "composition=%f value=%f min=%f max=%f\n",
		  compnodelist[j], apot_table.chempot[ntypes + j],
		  apot_table.pmin[apot_table.number][ntypes + j],
		  apot_table.pmax[apot_table.number][ntypes + j]);
      }
    }
  }
  exit(2);
}

#endif
#endif
