/****************************************************************
 *
 * utils.c: potfit utilities (vectors and memory management)
 *
 ****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
 *
 * http://potfit.sourceforge.net/
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
 *****************************************************************/

#include "potfit.h"

#include <stdint.h>

/* 32-bit */
#if UINTPTR_MAX == 0xffffffff
#if !defined(ACML)
#include <mkl_vml.h>
#endif  // ACML
#define _32BIT

/* 64-bit */
#elif UINTPTR_MAX == 0xffffffffffffffff
#if !defined(ACML)
#include <mkl_vml.h>
#elif defined(ACML4)
#include <acml_mv.h>
#elif defined(ACML5)
#include <amdlibm.h>
#endif  // ACML

/* wtf */
#else
#error Unknown integer size

#endif  // UINTPTR_MAX

#include "utils.h"

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
 *  higher powers in one and more dimensions
 *
 ****************************************************************/

void power_1(double* result, const double* x, const double* y)
{
#if defined(_32BIT)
  *result = pow(*x, *y);
#else
#if !defined(ACML)
  vdPow(1, x, y, result);
#elif defined(ACML4)
  *result = fastpow(*x, *y);
#elif defined(ACML5)
  *result = pow(*x, *y);
#endif  // ACML
#endif  // _32BIT
}

void power_m(int dim, double* result, const double* x, const double* y)
{
#if defined(_32BIT)
  int i = 0;
  for (i = 0; i < dim; i++)
    result[i] = pow(x[i], y[i]);
#else
#if !defined(ACML)
  vdPow(dim, x, y, result);
#elif defined(ACML4)
  int i;
  for (i = 0; i < dim; i++)
    *(result + i) = fastpow(*(x + i), *(y + i));
#elif defined(ACML5)
  int i;
  for (i = 0; i < dim; i++)
    *(result + i) = pow(*(x + i), *(y + i));
#endif  // ACML
#endif  // _32BIT
}

#if defined(_32BIT)
#undef _32BIT
#endif  // _32BIT
