/****************************************************************
 *
 * bracket.c: Brackets a minimum of a function.
 *
 ****************************************************************
 *
 * Copyright 1996, 1997, 1998, 1999, 2000 - Fabrice Rossi (gsl/min/bracketing.c)
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

#include "bracket.h"
#include "force.h"
#include "memory.h"
#include "utils.h"

void bracket(double* x_lower, double* x_minimum, double* x_upper,
             double* f_lower, double* f_minimum, double* f_upper,
             double* f_vec1, double* f_vec2)
{
  double f_left = *f_lower;
  double f_right = *f_upper;
  double f_center = 0.0;

  double x_left = *x_lower;
  double x_right = *x_upper;
  double x_center = 0.0;

  double* p_left = NULL;
  double* p_right = NULL;
  double* p_center = NULL;
  double* p_temp = NULL;

  int last = 0; /* indicates whether upwards is left or right */
  long nb_eval = 0;

  static double* vecu = NULL;
  static double* f_vec3 = NULL;

  if (vecu == NULL)
    vecu = (double*)Malloc(g_calc.ndimtot * sizeof(double));

  if (f_vec3 == NULL)
    f_vec3 = (double*)Malloc(g_calc.mdim * sizeof(double));

  p_left = f_vec1;
  p_right = f_vec2;
  p_center = f_vec3;

  if (f_right >= f_left) {
    x_center = x_left;
    f_center = f_left;
    SWAP(p_center, p_left, p_temp);
    x_left = -(x_right - x_center) / CGOLD + x_right;
    nb_eval++;
    for (int i = 0; i < g_calc.ndimtot; i++)
      vecu[i] = xicom[i] + x_left * delcom[i]; /*set vecu */
    f_left = calc_forces(vecu, p_left, 0);
  } else {
    x_center = x_right;
    f_center = f_right;
    SWAP(p_center, p_right, p_temp);
    x_right = (x_center - x_left) / CGOLD + x_left;
    nb_eval++;
    for (int i = 0; i < g_calc.ndimtot; i++)
      vecu[i] = xicom[i] + x_right * delcom[i]; /*set vecu */
    f_right = calc_forces(vecu, p_right, 0);
  }

  do {
    if (f_center < f_left) {
      if (f_center < f_right) {
        /* SUCCESS: Minimum in bracket! */
        *x_lower = x_left;
        *x_upper = x_right;
        *x_minimum = x_center;
        *f_lower = f_left;
        *f_upper = f_right;
        *f_minimum = f_center;
        for (int i = 0; i < g_calc.mdim; i++)
          f_vec1[i] = p_center[i];
        return;
      } else if (f_center > f_right) {
        /* OK, go right! */
        x_left = x_center;
        f_left = f_center;
        SWAP(p_left, p_center, p_temp);
        x_center = x_right;
        f_center = f_right;
        SWAP(p_center, p_right, p_temp);
        x_right = (x_center - x_left) / CGOLD + x_left;
        nb_eval++;
        for (int i = 0; i < g_calc.ndimtot; i++)
          vecu[i] = xicom[i] + x_right * delcom[i];
        f_right = calc_forces(vecu, p_right, 0);
      } else { /* f_center == f_right */

/* Pathological: Search between center and right */
/* This means a change from original algorithm */
#if defined(DEBUG)
        warning("Pathological  @%li %f %f %f! center-right!\n", nb_eval, x_left,
                x_center, x_right);
#endif  // DEBUG
        x_right = (x_right - x_left) * CGOLD + x_right;
        nb_eval++;
        for (int i = 0; i < g_calc.ndimtot; i++)
          vecu[i] = xicom[i] + x_right * delcom[i];
        f_right = calc_forces(vecu, p_right, 0);
        last = 1;
      }
    } else if (f_center > f_left)
    /* Search to the left */
    {
      x_right = x_center;
      f_right = f_center;
      SWAP(p_right, p_center, p_temp);
      x_center = x_left;
      f_center = f_left;
      SWAP(p_center, p_left, p_temp);
      x_left = -(x_right - x_center) / CGOLD + x_right;
      nb_eval++;
      for (int i = 0; i < g_calc.ndimtot; i++)
        vecu[i] = xicom[i] + x_left * delcom[i];
      f_left = calc_forces(vecu, p_left, 0);
    } else { /* f_center == f_left */

      if (f_center < f_right) {
/* between center and left */
#if defined(DEBUG)
        warning("Pathological  @%li %f %f %f! center-left!\n", nb_eval, x_left,
                x_center, x_right);
#endif  // DEBUG
        x_left = -(x_right - x_left) * CGOLD + x_left;
        nb_eval++;
        for (int i = 0; i < g_calc.ndimtot; i++)
          vecu[i] = xicom[i] + x_left * delcom[i];
        f_left = calc_forces(vecu, p_left, 0);
        last = 2;
      } else if (f_center > f_right) {
        /* Search to the right */
        x_left = x_center;
        f_left = f_center;
        SWAP(p_left, p_center, p_temp);
        x_center = x_right;
        f_center = f_right;
        SWAP(p_center, p_right, p_temp);
        x_right = (x_center - x_left) / CGOLD + x_left;
        nb_eval++;
        for (int i = 0; i < g_calc.ndimtot; i++)
          vecu[i] = xicom[i] + x_right * delcom[i];
        f_right = calc_forces(vecu, p_right, 0);
      } else { /* f_center==f_left==f_right */

        /* Kind of pathological case: go left/right in turns */
        if (last == 2) {
/* go further to left, it goes up towards the right */
#if defined(DEBUG)
          warning("Pathological  @%li %f %f %f! Go left!\n", nb_eval, x_left,
                  x_center, x_right);
#endif  // DEBUG
          x_left = -(x_right - x_left) / CGOLD + x_left;
          nb_eval++;
          for (int i = 0; i < g_calc.ndimtot; i++)
            vecu[i] = xicom[i] + x_left * delcom[i];
          f_left = calc_forces(vecu, p_left, 0);
          last = 1;
        } else { /* go further to the right, to left it went up */

#if defined(DEBUG)
          warning("Pathological @%li %f %f %f! Go right!\n", nb_eval, x_left,
                  x_center, x_right);
#endif  // DEBUG
          x_right = (x_right - x_left) / CGOLD + x_right;
          nb_eval++;
          for (int i = 0; i < g_calc.ndimtot; i++)
            vecu[i] = xicom[i] + x_right * delcom[i];
          f_right = calc_forces(vecu, p_right, 0);
          last = 2;
        }
      }
    }
  } while (nb_eval < MAX_IT);
#if defined(DEBUG)
  error(0, "Problems with bracketing minimum in %li tries:\n", nb_eval);
  error(1, "F(%.16g)=%.16g, F(%.16g)=%.16g, F(%.16g)=%.16g.\n", x_left, f_left,
        x_center, f_center, x_right, f_right);
#else
  error(1, "Problems with bracketing of minimum, aborting\n");
#endif  // DEBUG
}
