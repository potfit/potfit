/****************************************************************
 *
 * utils.c: potfit utilities (vectors and memory management)
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
 *****************************************************************/

#include "potfit.h"

#include <stdint.h>

/* 32-bit */
#if UINTPTR_MAX == 0xffffffff
#ifndef ACML
#include <mkl_vml.h>
#endif /* ACML */
#define _32BIT

/* 64-bit */
#elif UINTPTR_MAX == 0xffffffffffffffff
#ifndef ACML
#include <mkl_vml.h>
#elif defined ACML4
#include <acml_mv.h>
#elif defined ACML5
#include <amdlibm.h>
#endif /* ACML */

/* wtf */
#else
#error Unknown integer size

#endif /* UINTPTR_MAX */

#include "utils.h"

int  *vect_int(long dim)
{
  int  *vect, i;
  vect = (int *)malloc((size_t) (dim * sizeof(int)));
  if (vect == NULL)
    error(1, "Error in integer vector allocation");
  for (i = 0; i < dim; i++)
    vect[i] = 0;

  return vect;
}

double *vect_double(long dim)
{
  double *vect;
  int   i;
  vect = (double *)malloc((size_t) (dim * sizeof(double)));
  if (vect == NULL)
    error(1, "Error in double vector allocation");
  for (i = 0; i < dim; i++)
    vect[i] = 0.0;

  return vect;
}

double **mat_double(long rowdim, long coldim)
{
  long  i;
  double **matrix;

  /* matrix: array of array of pointers */
  /* matrix: pointer to rows */
  matrix = (double **)malloc((size_t) rowdim * sizeof(double *));
  if (matrix == NULL)
    error(1, "Error in double matrix row allocation");

  /* matrix[0]: pointer to elements */
  matrix[0] = (double *)malloc((size_t) rowdim * coldim * sizeof(double));
  if (matrix[0] == NULL)
    error(1, "Error in double matrix element allocation");

  for (i = 1; i < rowdim; i++)
    matrix[i] = matrix[i - 1] + coldim;

  int   j, k;
  for (j = 0; j < rowdim; j++)
    for (k = 0; k < coldim; k++)
      matrix[j][k] = 0.0;

  return matrix;
}

void free_vect_double(double *vect)
{
  free(vect);
}

void free_vect_int(int *vect)
{
  free(vect);
}

void free_mat_double(double **matrix)
{
  free(matrix[0]);
  free(matrix);
}

void init_atom(atom_t *atom)
{
  atom->type = 0;
  atom->num_neigh = 0;
  atom->pos.x = 0.0;
  atom->pos.y = 0.0;
  atom->pos.z = 0.0;
  atom->force.x = 0.0;
  atom->force.y = 0.0;
  atom->force.z = 0.0;
  atom->absforce = 0.0;
  atom->conf = 0;

#ifdef CONTRIB
  atom->contrib = 0;
#endif /* CONTRIB */

#if defined EAM || defined ADP || defined MEAM
  atom->rho = 0.0;
  atom->gradF = 0.0;
#ifdef TBEAM
  atom->rho_s = 0.0;
  atom->gradF_s = 0.0;
#endif /* TBEAM */
#endif /* EAM || ADP || MEAM */

#ifdef ADP
  atom->mu.x = 0.0;
  atom->mu.y = 0.0;
  atom->mu.z = 0.0;
  atom->lambda.xx = 0.0;
  atom->lambda.yy = 0.0;
  atom->lambda.zz = 0.0;
  atom->lambda.xy = 0.0;
  atom->lambda.yz = 0.0;
  atom->lambda.zx = 0.0;
  atom->nu = 0.0;
#endif /* ADP */

#ifdef DIPOLE
  atom->E_stat.x = 0.0;
  atom->E_stat.y = 0.0;
  atom->E_stat.z = 0.0;
  atom->p_sr.x = 0.0;
  atom->p_sr.y = 0.0;
  atom->p_sr.z = 0.0;
  atom->E_ind.x = 0.0;
  atom->E_ind.y = 0.0;
  atom->E_ind.z = 0.0;
  atom->p_ind.x = 0.0;
  atom->p_ind.y = 0.0;
  atom->p_ind.z = 0.0;
  atom->E_old.x = 0.0;
  atom->E_old.y = 0.0;
  atom->E_old.z = 0.0;
  atom->E_tot.x = 0.0;
  atom->E_tot.y = 0.0;
  atom->E_tot.z = 0.0;
#endif /* DIPOLE */

#ifdef THREEBODY
  atom->num_angles = 0;
#ifdef MEAM
  atom->rho_eam = 0.0;
#endif /* MEAM */
#endif /* MANYBODY */

  atom->neigh = NULL;
#ifdef THREEBODY
  atom->angle_part = NULL;
#endif /* THREEBODY */
}

void init_neigh(neigh_t *neigh)
{
  int   i = 0;

  neigh->type = 0;
  neigh->nr = 0;
  neigh->r = 0.0;
  neigh->r2 = 0.0;
  neigh->inv_r = 0.0;
  neigh->dist.x = 0.0;
  neigh->dist.y = 0.0;
  neigh->dist.z = 0.0;
  neigh->dist_r.x = 0.0;
  neigh->dist_r.y = 0.0;
  neigh->dist_r.z = 0.0;
  for (i = 0; i < SLOTS; i++) {
    neigh->slot[i] = 0;
    neigh->shift[i] = 0.0;
    neigh->step[i] = 0.0;
    neigh->col[i] = 0;
  }

#ifdef ADP
  neigh->sqrdist.xx = 0.0;
  neigh->sqrdist.yy = 0.0;
  neigh->sqrdist.zz = 0.0;
  neigh->sqrdist.xy = 0.0;
  neigh->sqrdist.yz = 0.0;
  neigh->sqrdist.zx = 0.0;
  neigh->u_val = 0.0;
  neigh->u_grad = 0.0;
  neigh->w_val = 0.0;
  neigh->w_grad = 0.0;
#endif /* APOT */

#ifdef COULOMB
  neigh->fnval_el = 0.0;
  neigh->grad_el = 0.0;
  neigh->ggrad_el = 0.0;
#endif /* COULOMB */

#ifdef THREEBODY
  neigh->f = 0.0;
  neigh->df = 0.0;
  neigh->ijk_start = 0;
#endif /* THREEBODY */

#ifdef MEAM
  neigh->drho = 0.0;
#endif /* MEAM */

#ifdef TERSOFF
  neigh->dzeta.x = 0.0;
  neigh->dzeta.y = 0.0;
  neigh->dzeta.z = 0.0;
#endif /* TERSOFF */
}

#ifdef THREEBODY
void init_angle(angle_t * angle)
{
  angle->cos = 0.0;

#ifdef MEAM
  angle->slot = 0;
  angle->shift = 0.0;
  angle->step = 0.0;
  angle->g = 0.0;
  angle->dg = 0.0;
#endif /* MEAM */
}
#endif /* THREEBODY */

void reg_for_free(void *p, char *name, ...)
{
  va_list ap;

  pointer_names = (char **)realloc(pointer_names, (num_pointers + 1) * sizeof(char *));
  pointer_names[num_pointers] = (char *)malloc((strlen(name) + 10) * sizeof(char));
  va_start(ap, name);
  vsprintf(pointer_names[num_pointers], name, ap);
  va_end(ap);
  all_pointers = (void **)realloc(all_pointers, (num_pointers + 1) * sizeof(void *));
  all_pointers[num_pointers] = p;
  num_pointers++;
}

void free_all_pointers()
{
  int   i;

  for (i = (num_pointers - 1); i >= 0; i--) {
    free(all_pointers[i]);
    free(pointer_names[i]);
  }
  free(all_pointers);
  free(pointer_names);
}

/* vector product */
vector vec_prod(vector u, vector v)
{
  vector w;

  w.x = u.y * v.z - u.z * v.y;
  w.y = u.z * v.x - u.x * v.z;
  w.z = u.x * v.y - u.y * v.x;

  return w;
}

/****************************************************************
 *
 *  double eqdist(): Returns an equally distributed random number in [0,1[
 * 	Uses dsfmt PRNG to generate a random number.
 *
 ****************************************************************/

inline double eqdist()
{
  return dsfmt_genrand_close_open(&dsfmt);
}

/****************************************************************
 *
 *  double normdist(): Returns a normally distributed random variable
 *
 ****************************************************************/

double normdist()
{
  static int have = 0;
  static double nd2;
  double x1, x2, sqr, cnst;

  if (!(have)) {
    do {
      x1 = 2.0 * eqdist() - 1.0;
      x2 = 2.0 * eqdist() - 1.0;
      sqr = x1 * x1 + x2 * x2;
    } while (!(sqr <= 1.0 && sqr > 0));
    /* Box Muller Transformation */
    cnst = sqrt(-2.0 * log(sqr) / sqr);
    nd2 = x2 * cnst;
    have = 1;
    return x1 * cnst;
  } else {
    have = 0;
    return nd2;
  }
}

/****************************************************************
 *
 *  square functions for integer and double values
 *
 ****************************************************************/

inline int isquare(int i)
{
  return i * i;
}

/* modified */
/* inline double dsquare(double d)
*/
double dsquare(double d)
/* modified ends */
{
  return d * d;
}

/****************************************************************
 *
 *  higher powers in one and more dimensions
 *
 ****************************************************************/

void power_1(double *result, double *x, double *y)
{
#ifdef _32BIT
  *result = pow(*x, *y);
#else
#ifndef ACML
  vdPow(1, x, y, result);
#elif defined ACML4
  *result = fastpow(*x, *y);
#elif defined ACML5
  *result = pow(*x, *y);
#endif /* ACML */
#endif /* _32BIT */
}

void power_m(int dim, double *result, double *x, double *y)
{
#ifdef _32BIT
  int   i = 0;
  for (i = 0; i < dim; i++)
    result[i] = pow(x[i], y[i]);
#else
#ifndef ACML
  vdPow(dim, x, y, result);
#elif defined ACML4
  int   i;
  for (i = 0; i < dim; i++)
    *(result + i) = fastpow(*(x + i), *(y + i));
#elif defined ACML5
  int   i;
  for (i = 0; i < dim; i++)
    *(result + i) = pow(*(x + i), *(y + i));
#endif /* ACML */
#endif /* _32BIT */
}

#if defined APOT && defined EVO

/****************************************************************
 *
 *  quicksort algorithm for opposition-based diff_evo
 *
 ****************************************************************/

void quicksort(double *x, int low, int high, double **p)
{
  int   newIndex;
  if (low < high) {
    int   index = (low + high) / 2;
    newIndex = partition(x, low, high, index, p);
    quicksort(x, low, newIndex - 1, p);
    quicksort(x, newIndex + 1, high, p);
  }
}

int partition(double *x, int low, int high, int index, double **p)
{
  int   i, store;
  double ind_val = x[index], temp;

  SWAP(x[index], x[high], temp);
  swap_population(p[index], p[high]);

  store = low;

  for (i = low; i < high; i++)
    if (x[i] <= ind_val) {
      SWAP(x[i], x[store], temp);
      swap_population(p[i], p[store]);
      store++;
    }
  SWAP(x[store], x[high], temp);
  swap_population(p[store], p[high]);

  return store;
}

void swap_population(double *a, double *b)
{
  int   i;
  double temp;
  for (i = 0; i < ndimtot + 2; i++) {
    SWAP(a[i], b[i], temp);
  }
}

#endif /* APOT && EVO */

#ifdef _32BIT
#undef _32BIT
#endif /* _32BIT */
