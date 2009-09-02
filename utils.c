/****************************************************************
*
*  utils.c: potfit utilities (vectors etc.)
*
*****************************************************************/
/*
*   Copyright 2002-2005 Peter Brommer
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
* $Revision: 1.3 $
* $Date: 2009/09/02 14:16:19 $
*****************************************************************/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "potfit.h"
#include "utils.h"

int  *vect_int(long dim)
{
  int  *vect;
  vect = (int *)malloc((size_t) (dim * sizeof(int)));
  if (vect == NULL)
    error("Error in integer vector allocation");
  return vect;
}

real *vect_real(long dim)
{
  real *vect;
  vect = (real *)malloc((size_t) (dim * sizeof(real)));
  if (vect == NULL)
    error("Error in real vector allocation");
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
