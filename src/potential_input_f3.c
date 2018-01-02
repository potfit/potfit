/****************************************************************
 *
 * potential_input_f3.c: Routines for reading a potential table
 *
 ****************************************************************
 *
 * Copyright 2002-2018 - the potfit development team
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

#include "memory.h"
#include "potential_input.h"
#include "utils.h"

#if defined(APOT)

void read_pot_table3(char const* a, FILE* b, potential_state* c)
{
  error(1, "Unsupported potential format in %s", a);
}

#else

void init_calc_table3();

/****************************************************************
 *
 *  read potential in third format:
 *
 *  Sampling points are equidistant.
 *
 *  Header:  one line for each function with
 *           rbegin rstart npoints
 *
 *  Table: Function values at sampling points,
 *         functions separated by blank lines
 *
 *  parameters:
 * 	pot_table_t * ... pointer to the potential table
 *  	int ... number of potential functions
 *  	char * ... name of the potential file (for error messages)
 *  	FILE * ... open file handle of the potential file
 *
 ****************************************************************/

void read_pot_table3(char const* potential_filename, FILE* pfile,
                     potential_state* pstate)
{
  int nvals[pstate->num_pots];

  pot_table_t* pt = &g_pot.opt_pot;

  /* read the info block of the function table */
  for (int i = 0; i < pstate->num_pots; i++) {
    if (3 > fscanf(pfile, "%lf %lf %d", &pt->begin[i], &pt->end[i], &nvals[i]))
      error(1, "Premature end of potential file %s\n", pstate->filename);

    pt->step[i] = (pt->end[i] - pt->begin[i]) / (nvals[i] - 1);
    pt->invstep[i] = 1.0 / pt->step[i];

    /* in the two slots between last[i-1] and first[i] the gradients
       of the respective functions are stored */
    if (i == 0)
      pt->first[i] = 2;
    else
      pt->first[i] = pt->last[i - 1] + 3;

    pt->last[i] = pt->first[i] + nvals[i] - 1;
    pt->len = pt->first[i] + nvals[i];
  }

  /* allocate the function table */
  pt->table = (double*)Malloc(pt->len * sizeof(double));
  pt->xcoord = (double*)Malloc(pt->len * sizeof(double));
  pt->d2tab = (double*)Malloc(pt->len * sizeof(double));
  pt->idx = (int*)Malloc(pt->len * sizeof(int));

  /* input loop */
  double* val = pt->table;
  int k = 0;
  int l = 0;

  /* read pair potentials */
  for (int i = 0; i < g_calc.paircol; i++) {
    if (pstate->have_gradient) {
      if (2 > fscanf(pfile, "%lf %lf\n", val, val + 1))
        error(1, "Premature end of potential file %s\n", pstate->filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.0;
    }
    val += 2;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (int j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(pfile, "%lf\n", val))
        error(1, "Premature end of potential file %s\n", pstate->filename);
      else
        val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!g_pot.invar_pot[i]) && (j < nvals[i] - 1))
        pt->idx[k++] = l++;
      else
        l++;
    }
  }

#if defined EAM || defined ADP || defined MEAM
  /* read EAM transfer function rho(r) */
  for (int i = g_calc.paircol; i < g_calc.paircol + g_param.ntypes; i++) {
    if (pstate->have_gradient) {
      if (2 > fscanf(pfile, "%lf %lf\n", val, val + 1))
        error(1, "Premature end of potential file %s\n", pstate->filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.0;
    }
    val += 2;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (int j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(pfile, "%lf\n", val))
        error(1, "Premature end of potential file %s\n", pstate->filename);
      else
        val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!g_pot.invar_pot[i]) && (j < nvals[i] - 1))
        pt->idx[k++] = l++;
      else
        l++;
    }
  }

  /* read EAM embedding function F(n) */
  for (int i = g_calc.paircol + g_param.ntypes;
       i < g_calc.paircol + 2 * g_param.ntypes; i++) {
    if (pstate->have_gradient) {
      if (2 > fscanf(pfile, "%lf %lf\n", val, val + 1))
        error(1, "Premature end of potential file %s\n", pstate->filename);
    } else {
      *val = 1.e30;
      *(val + 1) = 1.e30;
    }
    val += 2;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (int j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(pfile, "%lf\n", val))
        error(1, "Premature end of potential file %s\n", pstate->filename);
      else
        val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if (!g_pot.invar_pot[i])
        pt->idx[k++] = l++;
      else
        l++;
    }
  }

#if defined(TBEAM)
  /* read TBEAM transfer function rho(r) for the s-band */
  for (int i = g_calc.paircol + 2 * g_param.ntypes;
       i < g_calc.paircol + 3 * g_param.ntypes; i++) {
    if (pstate->have_gradient) {
      if (2 > fscanf(pfile, "%lf %lf\n", val, val + 1))
        error(1, "Premature end of potential file %s\n", pstate->filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.0;
    }
    val += 2;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (int j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(pfile, "%lf\n", val))
        error(1, "Premature end of potential file %s\n", pstate->filename);
      else
        val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!g_pot.invar_pot[i]) && (j < nvals[i] - 1))
        pt->idx[k++] = l++;
      else
        l++;
    }
  }

  /* read TBEAM embedding function F(n) for the s-band */
  for (int i = g_calc.paircol + 3 * g_param.ntypes;
       i < g_calc.paircol + 4 * g_param.ntypes; i++) {
    if (pstate->have_gradient) {
      if (2 > fscanf(pfile, "%lf %lf\n", val, val + 1))
        error(1, "Premature end of potential file %s\n", pstate->filename);
    } else {
      *val = 1.e30;
      *(val + 1) = 1.e30;
    }
    val += 2;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (int j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(pfile, "%lf\n", val))
        error(1, "Premature end of potential file %s\n", pstate->filename);
      else
        val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if (!g_pot.invar_pot[i])
        pt->idx[k++] = l++;
      else
        l++;
    }
  }
#endif  // TBEAM
#endif  // EAM || ADP || MEAM

#if defined(ADP)
  /* read ADP dipole function u(r) */
  for (int i = g_calc.paircol + 2 * g_param.ntypes;
       i < 2 * (g_calc.paircol + g_param.ntypes); i++) {
    if (pstate->have_gradient) {
      if (2 > fscanf(pfile, "%lf %lf\n", val, val + 1))
        error(1, "Premature end of potential file %s\n", pstate->filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.0;
    }
    val += 2;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (int j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(pfile, "%lf\n", val))
        error(1, "Premature end of potential file %s\n", pstate->filename);
      else
        val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!g_pot.invar_pot[i]) && (j < nvals[i] - 1))
        pt->idx[k++] = l++;
      else
        l++;
    }
  }

  /* read adp quadrupole function w(r) */
  for (int i = 2 * (g_calc.paircol + g_param.ntypes);
       i < 3 * g_calc.paircol + 2 * g_param.ntypes; i++) {
    if (pstate->have_gradient) {
      if (2 > fscanf(pfile, "%lf %lf\n", val, val + 1))
        error(1, "Premature end of potential file %s\n", pstate->filename);
    } else {
      *val = 1.e30;
      *(val + 1) = 1.e30;
    }
    val += 2;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (int j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(pfile, "%lf\n", val))
        error(1, "Premature end of potential file %s\n", pstate->filename);
      else
        val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!g_pot.invar_pot[i]) && (j < nvals[i] - 1))
        pt->idx[k++] = l++;
      else
        l++;
    }
  }
#endif  // ADP

#if defined(MEAM)
  /* read in second pair pot    f */
  for (int i = g_calc.paircol + 2 * g_param.ntypes;
       i < 2 * g_calc.paircol + 2 * g_param.ntypes; i++) {
    /* read gradient */
    if (pstate->have_gradient) {
      if (2 > fscanf(pfile, "%lf %lf\n", val, val + 1))
        error(1, "Premature end of potential file %s\n", pstate->filename);
      ;
    } else {
      *val = 1e30;
      *(val + 1) = 0.0;
    }
    val += 2;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    for (int j = 0; j < nvals[i]; j++) { /* read values */
      if (1 > fscanf(pfile, "%lf\n", val)) {
        error(1, "Premature end of potential file %s\n", pstate->filename);
        ;
      } else
        val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
// Clamp first spline knot in first f_ij potential only
// to remove degeneracy of f*f*g where f' = f/b and g' = b^2*g
#if !defined(MEAMf)
      if ((!g_pot.invar_pot[i]) &&
          (j < nvals[i] - 1 &&
           (j != 0 || i != g_calc.paircol + 2 * g_param.ntypes)))
#else
      if (!invar_pot[i])
#endif  // MEAMf
        pt->idx[k++] = l++;
      else
        l++;
    }
  }

  /* read in angl part */
  for (int i = 2 * g_calc.paircol + 2 * g_param.ntypes;
       i < 2 * g_calc.paircol + 3 * g_param.ntypes; i++) {
    /* read gradient */
    if (pstate->have_gradient) {
      if (2 > fscanf(pfile, "%lf %lf\n", val, val + 1))
        error(1, "Premature end of potential file %s\n", pstate->filename);
      ;
    } else {
      *val = 0;
      *(val + 1) = 0;
    }
    val += 2;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!g_pot.invar_pot[i]) && (g_pot.gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    for (int j = 0; j < nvals[i]; j++) { /* read values */
      if (1 > fscanf(pfile, "%lf\n", val)) {
        error(1, "Premature end of potential file %s\n", pstate->filename);
        ;
      } else
        val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if (!g_pot.invar_pot[i])
        pt->idx[k++] = l++;
      else
        l++;
    }
  }
#endif  // MEAM

  pt->idxlen = k;
  init_calc_table3();

  return;
}

/****************************************************************
 *
 *  init_calc_table3
 *      bla bla
 *
 ****************************************************************/

void init_calc_table3()
{
  g_pot.calc_pot.len = g_pot.opt_pot.len;
  g_pot.calc_pot.idxlen = g_pot.opt_pot.idxlen;
  g_pot.calc_pot.ncols = g_pot.opt_pot.ncols;
  g_pot.calc_pot.begin = g_pot.opt_pot.begin;
  g_pot.calc_pot.end = g_pot.opt_pot.end;
  g_pot.calc_pot.step = g_pot.opt_pot.step;
  g_pot.calc_pot.invstep = g_pot.opt_pot.invstep;
  g_pot.calc_pot.first = g_pot.opt_pot.first;
  g_pot.calc_pot.last = g_pot.opt_pot.last;
  g_pot.calc_pot.xcoord = g_pot.opt_pot.xcoord;
  g_pot.calc_pot.table = g_pot.opt_pot.table;
  g_pot.calc_pot.d2tab = g_pot.opt_pot.d2tab;
  g_pot.calc_pot.idx = g_pot.opt_pot.idx;
}

#endif  // APOT
