/****************************************************************
 *
 * utils.c: potfit utilities (vectors and memory management)
 *
 ****************************************************************
 *
 * Copyright 2002-2010 Peter Brommer, Daniel Schopf
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
 *****************************************************************/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "potfit.h"
#include "utils.h"

inline real sqrreal(real x)
{
  return x * x;
}

int  *vect_int(long dim)
{
  int  *vect, i;
  vect = (int *)malloc((size_t) (dim * sizeof(int)));
  if (vect == NULL)
    error("Error in integer vector allocation");
  for (i = 0; i < dim; i++)
    vect[i] = 0;
  return vect;
}

real *vect_real(long dim)
{
  real *vect;
  int   i;
  vect = (real *)malloc((size_t) (dim * sizeof(real)));
  if (vect == NULL)
    error("Error in real vector allocation");
  for (i = 0; i < dim; i++)
    vect[i] = 0.;
  return vect;
}

real **mat_real(long rowdim, long coldim)
{
  long  i;
  real **matrix;

  /* matrix: array of array of pointers */
  /* matrix: pointer to rows */
  matrix = (real **)malloc((size_t) rowdim * sizeof(real *));
  if (matrix == NULL)
    error("Error in real matrix row allocation");

  /* matrix[0]: pointer to elements */
  matrix[0] = (real *)malloc((size_t) rowdim * coldim * sizeof(real));
  if (matrix[0] == NULL)
    error("Error in real matrix element allocation");

  for (i = 1; i < rowdim; i++)
    matrix[i] = matrix[i - 1] + coldim;

  return matrix;
}

void free_vect_real(real *vect)
{
  free(vect);
}

void free_vect_int(int *vect)
{
  free(vect);
}

void free_mat_real(real **matrix)
{
  free(matrix[0]);
  free(matrix);
}

void reg_for_free(void *p, char *name)
{
  pointer_names =
    (char **)realloc(pointer_names, (num_pointers + 1) * sizeof(char *));
  pointer_names[num_pointers] =
    (char *)malloc((strlen(name) + 1) * sizeof(char));
  strcpy(pointer_names[num_pointers], name);
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
 *  real normdist(): Returns a normally distributed random variable
 * 	Uses dsfmt PRNG to generate a random number.
 *
 *****************************************************************/

real normdist()
{
  static int have = 0;
  static real nd2;
  real  x1, x2, sqr, cnst;

  if (!(have)) {
    do {
      x1 = 2.0 * dsfmt_genrand_close_open(&dsfmt) - 1.0;
      x2 = 2.0 * dsfmt_genrand_close_open(&dsfmt) - 1.0;
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

#ifdef APOT

/****************************************************************
 *
 *  quicksort algorithm for opposition-based diff_evo
 *
 *****************************************************************/

void quicksort(real *x, int low, int high, real **p)
{
  int   newIndex;
  if (low < high) {
    int   index = (low + high) / 2;
    newIndex = partition(x, low, high, index, p);
    quicksort(x, low, newIndex - 1, p);
    quicksort(x, newIndex + 1, high, p);
  }
}

int partition(real *x, int low, int high, int index, real **p)
{
  int   i, store;
  real  ind_val = x[index], temp;

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

void swap_population(real *a, real *b)
{
  int   i;
  real  temp;
  for (i = 0; i < ndimtot + 2; i++) {
    SWAP(a[i], b[i], temp);
  }
}

#endif /* APOT */
