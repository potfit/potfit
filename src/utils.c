/****************************************************************
 *
 * utils.c: potfit utilities (vectors and memory management)
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
 *****************************************************************/

#include <stdint.h>

#include "potfit.h"

#include "memory.h"
#include "utils.h"

#if UINTPTR_MAX == 0xffffffff

// 32-bit
#if defined(MKL)
#include <mkl_vml.h>
#endif  // MKL
#define _32BIT

#elif UINTPTR_MAX == 0xffffffffffffffff

// 64-bit
#if defined(MKL)
#include <mkl_vml.h>
#elif defined(__ACCELERATE__)
#include <Accelerate/Accelerate.h>
#else
#error No math library defined!
#endif

#else

// wtf
#error Unknown integer size

#endif  // UINTPTR_MAX

// vector product
vector vec_prod(vector u, vector v)
{
  vector w;

  w.x = u.y * v.z - u.z * v.y;
  w.y = u.z * v.x - u.x * v.z;
  w.z = u.x * v.y - u.y * v.x;

  return w;
}

#if defined(__ACCELERATE__)
static const int g_dim = 1;
#endif // __ACCELERATE__

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
#if defined(MKL)
  vdPow(1, x, y, result);
#elif defined(__ACCELERATE__)
  vvpow(result, y, x, &g_dim);
#endif
#endif  // _32BIT
}

void power_m(int dim, double* result, const double* x, const double* y)
{
#if defined(_32BIT)
  int i = 0;
  for (i = 0; i < dim; i++)
    result[i] = pow(x[i], y[i]);
#else
#if defined(MKL)
  vdPow(dim, x, y, result);
#elif defined(__ACCELERATE__)
  vvpow(result, y, x, &dim);
#endif
#endif  // _32BIT
}

#if defined(_32BIT)
#undef _32BIT
#endif  // _32BIT

char* fgets_potfit(char* buffer, int len, FILE* f)
{
  char* p = fgets(buffer, len, f);
  if (!p)
    return p;

  // iterate the buffer to find \r\n
  int found = 0;
  char* c = p;
  while (*c) {
    if (*c == '\r' && *(c+1) == '\n') {
      found = 1;
      break;
    }
    ++c;
  }

  if (!found)
    return p;

  warning("Your input files contain CRLF line endings!\n");

  char* new_p = malloc(len);

  int pos = 0;
  for (int i = 0; i < strlen(p); ++i) {
    if (p[i] != '\r') {
      new_p[pos++] = p[i];
    }
  }

  if (pos < len)
    new_p[pos] = '\0';

  memcpy(buffer, new_p, pos);

  free(new_p);

  return buffer;
}
