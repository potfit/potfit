/****************************************************************
 *
 * rescale.c: Routines used to automatically rescale
 *	EAM potential.
 *
 *****************************************************************
 *
 * Copyright 2002-2014
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.sourceforge.net/
 *
 *****************************************************************
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

#if defined RESCALE && !defined APOT

#include "forces.h"
#include "splines.h"

/* Doesn't make much sense without EAM or ADP  */

#if defined EAM || defined ADP

/****************************************************************
 *
 * rescale: Routine used to automatically rescale
 *     EAM potential. Flag indicates whether to force update...
 *     upper is upper limit of electron density.
 *
 ****************************************************************/

double rescale(pot_table_t *pt, double upper, int flag)
{
  int mincol, maxcol, col, col2, first, vals, h, i, j, typ1, typ2, sign,
      dimneuxi;
  double *xi, *neuxi, *neuord, *neustep, *maxrho, *minrho, *left, *right;
  atom_t *atom;
  neigh_t *neigh;
  double fnval, pos, grad, a;
  double min = 1e100, max = -1e100;

  xi = pt->table;
  dimneuxi = pt->last[g_calc.paircol + 2 * g_param.ntypes - 1] -
             pt->last[g_calc.paircol + g_param.ntypes - 1];
  neuxi = (double *)malloc(dimneuxi * sizeof(double));
  neuord = (double *)malloc(dimneuxi * sizeof(double));
  neustep = (double *)malloc(g_param.ntypes * sizeof(double));
  maxrho = (double *)malloc(g_param.ntypes * sizeof(double));
  minrho = (double *)malloc(g_param.ntypes * sizeof(double));
  left = (double *)malloc(g_param.ntypes * sizeof(double));
  right = (double *)malloc(g_param.ntypes * sizeof(double));
  for (i = 0; i < g_param.ntypes; i++) {
    maxrho[i] = -1e100;
    minrho[i] = 1e100;
  }
  /* find Max/Min rho  */
  /* init splines - better safe than sorry */
  /* init second derivatives for splines */
  for (col = 0; col < g_calc.paircol; col++) { /* just pair potentials */
    first = pt->first[col];
    if (g_pot.format == 3 || g_pot.format == 0)
      spline_ed(pt->step[col], pt->table + first, pt->last[col] - first + 1,
                *(pt->table + first - 2), 0.0, pt->d2tab + first);
    else /* format == 4 ! */
      spline_ne(pt->xcoord + first, pt->table + first,
                pt->last[col] - first + 1, *(pt->table + first - 2), 0.0,
                pt->d2tab + first);
  }
  for (col = g_calc.paircol; col < g_calc.paircol + g_param.ntypes;
       col++) { /* rho */
    first = pt->first[col];
    if (g_pot.format == 3)
      spline_ed(pt->step[col], xi + first, pt->last[col] - first + 1,
                *(xi + first - 2), 0.0, pt->d2tab + first);
    else /* format == 4 ! */
      spline_ne(pt->xcoord + first, xi + first, pt->last[col] - first + 1,
                *(xi + first - 2), 0.0, pt->d2tab + first);
  }
  for (col = g_calc.paircol + g_param.ntypes;
       col < g_calc.paircol + 2 * g_param.ntypes; col++) { /* F */
    first = pt->first[col];
    /* gradient 0 at r_cut */
    if (g_pot.format == 3)
      spline_ed(pt->step[col], xi + first, pt->last[col] - first + 1,
                *(xi + first - 2), *(xi + first - 1), pt->d2tab + first);
    else /* format == 4 */
      spline_ne(pt->xcoord + first, xi + first, pt->last[col] - first + 1,
                *(xi + first - 2), *(xi + first - 1), pt->d2tab + first);
  }

  /* re-calculate atom_rho (might be a waste...) */
  for (h = 0; h < g_config.nconf; h++) {
    for (i = 0; i < g_config.inconf[h]; i++)
      g_config.atoms[g_config.cnfstart[h] + i].rho = 0.0;
    for (i = 0; i < g_config.inconf[h]; i++) {
      atom = g_config.atoms + i + g_config.cnfstart[h];
      typ1 = atom->type;
      for (j = 0; j < atom->num_neigh; j++) {
        neigh = atom->neigh + j;
        if (neigh->nr > i + g_config.cnfstart[h]) {
          typ2 = neigh->type;
          col2 = g_calc.paircol + typ2;
          if (typ2 == typ1) {
            if (neigh->r < pt->end[col2]) {
              fnval = splint_dir(pt, xi, neigh->slot[1], neigh->shift[1],
                                 neigh->step[1]);
              atom->rho += fnval;
              g_config.atoms[neigh->nr].rho += fnval;
            }
          } else {
            col = g_calc.paircol + typ1;
            if (neigh->r < pt->end[col2]) {
              atom->rho += splint_dir(pt, xi, neigh->slot[1], neigh->shift[1],
                                      neigh->step[1]);
            }
            if (neigh->r < pt->end[col])
              g_config.atoms[neigh->nr].rho +=
                  (*g_splint)(pt, xi, col, neigh->r);
          }
        }
      }
      maxrho[typ1] = MAX(maxrho[typ1], atom->rho);
      minrho[typ1] = MIN(minrho[typ1], atom->rho);
    }
  }
  for (i = 0; i < g_param.ntypes; i++) {
    if (maxrho[i] > max) {
      max = maxrho[i];
      maxcol = i;
    }
    if (minrho[i] < min) {
      min = minrho[i];
      mincol = i;
    }
  }
  /* determine dominant side */
  sign = (max >= -min) ? 1 : -1;

  /* determine new left and right boundary, add 40 percent... */
  for (i = 0; i < g_param.ntypes; i++) {
    j = g_calc.paircol + g_param.ntypes + i;
    left[i] = minrho[i] - 0.3 * pt->step[j];
    right[i] = maxrho[i] + 0.3 * pt->step[j];
    /* is expansion necessary? */
    if (flag || minrho[i] - pt->begin[j] < 0.0 ||
        minrho[i] - pt->begin[j] > 0.95 * pt->step[j] ||
        maxrho[i] - pt->end[j] > 0 ||
        maxrho[i] - pt->end[j] < -0.95 * pt->step[j])
      flag = 1;
  }

  /* determine scaling factor  */
  a = (sign == 1) ? upper / right[maxcol] : upper / left[mincol];

  if (flag || fabs(a) > 1.05 || fabs(a) < 0.95) flag = 1;

  /* update needed? */
  if (!flag) return 0.0; /* no */

  /* Let's update... */
  /* expand potential  */
  h = 0;
  for (i = 0; i < g_param.ntypes; i++) {
    col = g_calc.paircol + g_param.ntypes + i; /* 1. embedding function */
    vals = pt->last[col] - pt->first[col];
    neustep[i] = (right[i] - left[i]) / (double)vals;
    pos = left[i];
    for (j = 0; j <= vals; j++) {
      neuxi[h] = splint_ne(pt, xi, col, pos); /* inter- or extrapolation */

      neuord[h] = pos;
      h++;
      pos += neustep[i];
    }
    /* correct gradient */
    if (*(xi + pt->first[col] - 2) < 1.e30)
      *(xi + pt->first[col] - 2) = splint_grad_ne(pt, xi, col, left[i]);
    if (*(xi + pt->first[col] - 1) < 1.e30)
      *(xi + pt->first[col] - 1) = splint_grad_ne(pt, xi, col, right[i]);
  }

  /* write back values */
  col = 0; /* first value to be changed */
  for (j = g_calc.paircol + g_param.ntypes;
       j < g_calc.paircol + 2 * g_param.ntypes; j++)
    for (i = pt->first[j]; i <= pt->last[j]; i++) {
      xi[i] = neuxi[col];
      pt->xcoord[i] = neuord[col];
      col++;
    }
  printf("Scaling factor %f\n", a);
  /* scale */
  for (i = g_calc.paircol; i < g_calc.paircol + g_param.ntypes; i++) {
    for (j = pt->first[i]; j <= pt->last[i]; j++) {
      pt->table[j] *= a;
    }
    if (*(xi + pt->first[i] - 2) < 1.e30) *(xi + pt->first[i] - 2) *= a;
  }

  /* rescale all embed. by a */
  if (sign == 1) {
    j = 0;
    for (i = g_calc.paircol + g_param.ntypes;
         i < g_calc.paircol + 2 * g_param.ntypes; i++) {
      pt->begin[i] = a * left[j];
      pt->end[i] = a * right[j];
      pt->step[i] = a * neustep[j];
      pt->invstep[i] = 1.0 / pt->step[i];
      /* gradient correction */
      if (xi[pt->first[i] - 2] < 1.e30) xi[pt->first[i] - 2] /= a;
      if (xi[pt->first[i] - 1] < 1.e30) xi[pt->first[i] - 1] /= a;
      pos = pt->begin[i];
      for (h = pt->first[i]; h <= pt->last[i]; h++) {
        pt->xcoord[h] = pos;
        pos += pt->step[i];
      }
      j++;
    }
  } else { /* reverse - a negativ */
    j = 0;
    for (i = g_calc.paircol + g_param.ntypes;
         i < g_calc.paircol + 2 * g_param.ntypes; i++) {
      pt->begin[i] = a * right[j];
      pt->end[i] = a * left[j];
      pt->step[i] = -a * neustep[j];
      pt->invstep[i] = 1.0 / pt->step[i];
      /* gradient correction and exchange */
      if (xi[pt->first[i] - 2] < 1.e30)
        grad = -xi[pt->first[i] - 2] / a;
      else
        grad = 1.e30;
      if (xi[pt->first[i] - 1] < 1.e30)
        xi[pt->first[i] - 2] = -xi[pt->first[i] - 1] / a;
      else
        xi[pt->first[i] - 2] = 1.e30;
      xi[pt->first[i] - 1] = grad;
      pos = pt->begin[i];
      for (h = pt->first[i]; h <= pt->last[i]; h++) {
        pt->xcoord[h] = pos;
        pos += pt->step[i];
      }
      j++;
    }
    h = 0;
    for (i = 0; i < g_param.ntypes; i++) { /* values in reverse order */
      col = g_calc.paircol + g_param.ntypes + i;
      for (j = pt->last[col]; j >= pt->first[col]; j--) {
        neuxi[h] = xi[j];
        h++;
      }
    }
    col = 0; /* and write back */
    for (j = g_calc.paircol + g_param.ntypes;
         j < g_calc.paircol + 2 * g_param.ntypes; j++)
      for (i = pt->first[j]; i <= pt->last[j]; i++) {
        xi[i] = neuxi[col];
        col++;
      }
  }
  /* re-initialise splines */
  for (col = g_calc.paircol; col < g_calc.paircol + g_param.ntypes;
       col++) { /* rho */
    first = pt->first[col];
    if (g_pot.format == 3)
      spline_ed(pt->step[col], xi + first, pt->last[col] - first + 1,
                *(xi + first - 2), 0.0, pt->d2tab + first);
    else /* format == 4 ! */
      spline_ne(pt->xcoord + first, xi + first, pt->last[col] - first + 1,
                *(xi + first - 2), 0.0, pt->d2tab + first);
  }

  for (col = g_calc.paircol + g_param.ntypes;
       col < g_calc.paircol + 2 * g_param.ntypes; col++) { /* F */
    first = pt->first[col];
    /* gradient 0 at r_cut */
    if (g_pot.format == 3)
      spline_ed(pt->step[col], xi + first, pt->last[col] - first + 1,
                *(xi + first - 2), *(xi + first - 1), pt->d2tab + first);
    else /* format == 4 */
      spline_ne(pt->xcoord + first, xi + first, pt->last[col] - first + 1,
                *(xi + first - 2), *(xi + first - 1), pt->d2tab + first);
  }

  /* correct gauge: U'(n_mean)=0 */
  for (i = 0; i < g_param.ntypes; i++) {
    g_todo.lambda[i] = (*g_splint_grad)(
        &g_pot.opt_pot, pt->table, g_calc.paircol + g_param.ntypes + i,
        0.5 * (pt->begin[g_calc.paircol + g_param.ntypes + i] +
               pt->end[g_calc.paircol + g_param.ntypes + i]));
  }
  for (i = 0; i < g_param.ntypes; i++)
    printf("lambda[%d] = %f\n", i, g_todo.lambda[i]);
  i = 0;

  for (col = 0; col < g_param.ntypes; col++)
    for (col2 = col; col2 < g_param.ntypes; col2++) {
      for (j = pt->first[i]; j <= pt->last[i]; j++)
        pt->table[j] +=
            (pt->xcoord[j] < pt->end[g_calc.paircol + col2]
                 ? g_todo.lambda[col] * splint_ne(pt, pt->table,
                                                  g_calc.paircol + col2,
                                                  pt->xcoord[j])
                 : 0.0) +
            (pt->xcoord[j] < pt->end[g_calc.paircol + col]
                 ? g_todo.lambda[col2] * splint_ne(pt, pt->table,
                                                   g_calc.paircol + col,
                                                   pt->xcoord[j])
                 : 0.0);
      /* Gradient */
      if (pt->table[pt->first[i] - 2] < 1e29) /* natural spline */
        pt->table[pt->first[i] - 2] +=
            (pt->begin[i] < pt->end[g_calc.paircol + col2]
                 ? g_todo.lambda[col] * (*g_splint_grad)(pt, pt->table,
                                                         g_calc.paircol + col2,
                                                         pt->begin[i])
                 : 0.0) +
            (pt->begin[i] < pt->end[g_calc.paircol + col]
                 ? g_todo.lambda[col2] * (*g_splint_grad)(pt, pt->table,
                                                          g_calc.paircol + col,
                                                          pt->begin[i])
                 : 0.0);
      if (pt->table[pt->first[i] - 1] < 1e29) /* natural spline */
        pt->table[pt->first[i] - 1] +=
            (pt->end[i] < pt->end[g_calc.paircol + col2]
                 ? g_todo.lambda[col] * (*g_splint_grad)(pt, pt->table,
                                                         g_calc.paircol + col2,
                                                         pt->end[i])
                 : 0.0) +
            (pt->end[i] < pt->end[g_calc.paircol + col]
                 ? g_todo.lambda[col2] * (*g_splint_grad)(pt, pt->table,
                                                          g_calc.paircol + col,
                                                          pt->end[i])
                 : 0.0);
      i++;
    }
  for (i = 0; i < g_param.ntypes; i++) {
    for (j = pt->first[g_calc.paircol + g_param.ntypes + i];
         j <= pt->last[g_calc.paircol + g_param.ntypes + i]; j++)
      pt->table[j] -= pt->xcoord[j] * g_todo.lambda[i];
    /* Gradients */
    if (pt->table[pt->first[g_calc.paircol + g_param.ntypes + i] - 2] <
        1e29) /* natural spline */
      pt->table[pt->first[g_calc.paircol + g_param.ntypes + i] - 2] -=
          g_todo.lambda[i];
    if (pt->table[pt->first[g_calc.paircol + g_param.ntypes + i] - 1] <
        1e29) /* natural spline */
      pt->table[pt->first[g_calc.paircol + g_param.ntypes + i] - 1] -=
          g_todo.lambda[i];
    g_todo.lambda[i] = 0.0;
  }

  /* init second derivatives for splines */
  for (col = 0; col < g_calc.paircol; col++) { /* just pair potentials */
    first = pt->first[col];
    if (g_pot.format == 3)
      spline_ed(pt->step[col], pt->table + first, pt->last[col] - first + 1,
                *(pt->table + first - 2), 0.0, pt->d2tab + first);
    else /* format == 4 ! */
      spline_ne(pt->xcoord + first, pt->table + first,
                pt->last[col] - first + 1, *(pt->table + first - 2), 0.0,
                pt->d2tab + first);
  }

  free(neuxi);
  free(neustep);
  free(maxrho);
  free(minrho);
  free(left);
  free(right);

  /* return factor */
  return a;
}

/****************************************************************
 *
 * embed_shift: Shift embedding function to U(0)=0
 *
 ****************************************************************/

void embed_shift(pot_table_t *pt)
{
  double shift;
  double *xi;
  int i, j, first;
  xi = pt->table;
  for (i = g_calc.paircol + g_param.ntypes;
       i < g_calc.paircol + 2 * g_param.ntypes; i++) {
    first = pt->first[i];
    /* init splines - better safe than sorry */
    /* gradient 0 at r_cut */
    /********* DANGER ****************/
    /** NOT FOOLPROOF ***************/
    if (pt->begin[i] <= 0) { /* 0 in domain of U(n) */
      if (g_pot.format == 3)
        spline_ed(pt->step[i], xi + first, pt->last[i] - first + 1,
                  *(xi + first - 2), *(xi + first - 1), pt->d2tab + first);
      else /* format == 4 ! */
        spline_ne(pt->xcoord + first, xi + first, pt->last[i] - first + 1,
                  *(xi + first - 2), *(xi + first - 1), pt->d2tab + first);
      shift = (*g_splint)(pt, xi, i, 0.0);
#ifdef DEBUG
      printf("shifting by %f\n", shift);
#endif /* DEBUG */
    } else
      shift = xi[first];
    for (j = first; j <= pt->last[i]; j++) xi[j] -= shift;
  }
}
#endif /* EAM || ADP */
#endif /* RESCALE && !APOT */
