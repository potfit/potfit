/****************************************************************
 *
 * utils.c: potfit utilities (vectors and memory management)
 *
 ****************************************************************
 *
 * Copyright 2002-2011
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
 *****************************************************************/

#ifndef ACML
#include <mkl_vml.h>
#else
#include <acml_mv.h>
#endif /* ACML */

#include "potfit.h"
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
    vect[i] = 0.;
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
      matrix[j][k] = 0.;

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

void reg_for_free(void *p, char *name, ...)
{
  va_list ap;

  pointer_names =
    (char **)realloc(pointer_names, (num_pointers + 1) * sizeof(char *));
  pointer_names[num_pointers] =
    (char *)malloc((strlen(name) + 10) * sizeof(char));
  va_start(ap, name);
  vsprintf(pointer_names[num_pointers], name, ap);
  va_end(ap);
  all_pointers =
    (void **)realloc(all_pointers, (num_pointers + 1) * sizeof(void *));
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
 * 	Uses dsfmt PRNG to generate a random number.
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

inline double dsquare(double d)
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
#ifndef ACML
  vdPow(1, x, y, result);
#else
  *result = fastpow(*x, *y);
#endif /* ACML */
}

void power_m(int dim, double *result, double *x, double *y)
{
#ifndef ACML
  vdPow(dim, x, y, result);
#else
  int   i;
  for (i = 0; i < dim; i++) {
    *(result + i) = fastpow(*(x + i), *(y + i));
  }
#endif /* ACML */
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
