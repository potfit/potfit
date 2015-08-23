/****************************************************************
 *
 * rescale_meam.c: Routines used to automatically rescale MEAM potentials
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

#if defined(RESCALE) && !defined(APOT)

#include "memory.h"
#include "splines.h"

/* Doesn't make much sense without MEAM */

#if defined(MEAM)

/****************************************************************
 *
 * rescale: Routine used to automatically rescale
 *     MEAM potential. Flag indicates whether to force update...
 *     upper is upper limit of electron density.
 *
 ****************************************************************/

double rescale(pot_table_t* pt, double upper, int flag)
{
  int i, j, h, col, col2, type1, type2, mincol, maxcol, first, vals, sign;
  double min = 1e100, max = -1e100;
  double pos, grad, a;
  atom_t* atom;
  neigh_t* neigh;

  int jj, kk, ijk;
  angle_t* angle;
  neigh_t* neigh_j, *neigh_k;

  /* Set potential array in xi */
  double* const xi = pt->table;

  /* Last index of F(type n) - Last index of rho(type n)
     The total number of splines and gradients for F and all of its
     types for an alloy */
  const int dimnewxi = pt->last[g_calc.paircol + 2 * g_param.ntypes - 1] -
                       pt->last[g_calc.paircol + g_param.ntypes - 1];

  /* Allocate memory for dynamic arrays that will store the new F */
  double* const newxi =
      (double*)Malloc_Local(dimnewxi * sizeof(double));  // Size of F array
  double* const neword =
      (double*)Malloc_Local(dimnewxi * sizeof(double));  // Size of F array
  double* const newstep =
      (double*)Malloc_Local(g_param.ntypes * sizeof(double));  // # of cols of F
  double* const maxrho =
      (double*)Malloc_Local(g_param.ntypes * sizeof(double));  // # of cols of F
  double* const minrho =
      (double*)Malloc_Local(g_param.ntypes * sizeof(double));  // # of cols of F
  double* const left =
      (double*)Malloc_Local(g_param.ntypes * sizeof(double));  // # of cols of F
  double* const right =
      (double*)Malloc_Local(g_param.ntypes * sizeof(double));  // # of cols of F

  // Initialize max and min rho's for each column in F
  for (i = 0; i < g_param.ntypes; i++)
  {
    maxrho[i] = -1e100;
    minrho[i] = 1e100;
  }

  // We need to find the max and min rho's for each column in F
  // Unfortunately we need to resolve for the rho's of each atom
  // since in MPI this info isn't passed to the root node
  //////////////////////////////////////////////////////////////////

  // Initialize the 2nd derivs for splines, so that we can interpolate
  // in the future
  for (col = 0; col < 2 * g_calc.paircol + 3 * g_param.ntypes; ++col)
  {
    // Pointer to first entry
    first = pt->first[col];

    // Get 2nd derivatives
    // step = width of spline knots (known as h)
    // xi+first = array with spline values
    // pt->last[col1] - first + 1 = num of spline pts
    // *(xi + first - 2) = value of endpoint gradient (default: 1e30)
    // *(xi + first - 1) = value of other endpoint gradient
    //          (default: phi=0.0, rho=0.0, F=1e30)
    // g_pot.calc_pot.d2tab + first = array to hold 2nd deriv
    spline_ed(pt->step[col], xi + first, pt->last[col] - first + 1, *(xi + first - 2),
              *(xi + first - 1), pt->d2tab + first);
  }

  // The dreaded recalculation of atom->rho for each atom
  // LOOP OVER EACH CONFIG
  for (h = 0; h < g_config.nconf; ++h)
  {
    // Reset the rho for each atom in config
    for (i = 0; i < g_config.inconf[h]; ++i)
      g_config.atoms[g_config.cnfstart[h] + i].rho = 0.0;

    // BEGIN LOOP OVER EACH ATOM i CALCULTING THE RHO
    for (i = 0; i < g_config.inconf[h]; ++i)
    {
      // Set temporary atom structure
      atom = g_config.atoms + i + g_config.cnfstart[h];

      // Get type of atom i
      type1 = atom->type;

      // LOOP OVER EACH NEIGHBOR OF ATOM i
      for (j = 0; j < atom->num_neigh; ++j)
      {
        // Store neighbor to temp variable
        neigh = atom->neigh + j;

        // Store type of neighbor
        type2 = neigh->type;

        // BEGIN COMPUTING rho values contributed from rho_ij
        // Store rho_ij column that is used for neighbor
        col2 = g_calc.paircol + type2;

        // Check that neighbor lies in range of rho_ij potential
        if (neigh->r < pt->end[col2])
        {
          // Compute rho value and store it for atom i
          atom->rho +=
              splint_dir(pt, xi, neigh->slot[1], neigh->shift[1], neigh->step[1]);
        }
        // BEGIN COMPUTING rho values for f_ij potential
        // Get column for atom j (it behaves similarly to pair potential, phi)
        // just offset by "g_calc.paircol + 2*ntypes"
        col2 = (type1 <= type2)
                   ? g_calc.paircol + 2 * g_param.ntypes + type1 * g_param.ntypes +
                         type2 - ((type1 * (type1 + 1)) / 2)
                   : g_calc.paircol + 2 * g_param.ntypes + type2 * g_param.ntypes +
                         type1 - ((type2 * (type2 + 1)) / 2);

        // Check that neighbor lies in range of f_ij potential
        if (neigh->r < pt->end[col2])
        {
          // Store the f(r_ij) value
          neigh->f = splint_dir(pt, xi, neigh->slot[2], neigh->shift[2], neigh->step[2]);
        }
        else
        {
          // Set f(r_ij) to 0 so it doesn't multiply in later
          neigh->f = 0;
        }

      }  // END OF LOOP OVER NEIGHBORS

      // Loop over every angle formed by neighbors
      // N(N-1)/2 possible combinations
      // Used in computing angular part g_ijk
      ijk = 0;  // count number of angles
      for (jj = 0; jj < atom->num_neigh - 1; ++jj)
      {
        // Get pointer to neighbor jj
        neigh_j = atom->neigh + jj;
        for (kk = jj + 1; kk < atom->num_neigh; ++kk)
        {
          // Store pointer to angular part (g)
          angle = atom->angle_part + ijk;

          // Get pointer to neighbor kk
          neigh_k = atom->neigh + kk;

          // The cos(theta) should always lie inside -1 ... 1
          // So store the g and g' without checking bounds
          angle->g = splint_dir(pt, xi, angle->slot, angle->shift, angle->step);

          // Sum up rho piece for atom i caused by j and k
          // f_ij * f_ik * m_ijk
          atom->rho += neigh_j->f * neigh_k->f * angle->g;
          ++ijk;
        }  // END OF INNER LOOP OVER TRIPLETS
      }    // END OF OUTER LOOP OVER TRIPLETS
    }      // END OF LOOP OVER ATOM i

    // BEGIN LOOP OVER EACH ATOM i FINDING MAX/MIN RHO
    for (i = 0; i < g_config.inconf[h]; ++i)
    {
      // Set temporary atom structure
      atom = g_config.atoms + i + g_config.cnfstart[h];

      // Get type of atom i
      type1 = atom->type;

      // Store max/min rho for this column of F seen by atom i
      maxrho[type1] = MAX(maxrho[type1], atom->rho);
      minrho[type1] = MIN(minrho[type1], atom->rho);
    }
  }  // END OF LOOP OVER CONFIGURATIONS

  // Loop through each alloy's F looking for max/min column
  for (i = 0; i < g_param.ntypes; ++i)
  {
    if (maxrho[i] > max)
    {
      max = maxrho[i];
      maxcol = i;
    }
    if (minrho[i] < min)
    {
      min = minrho[i];
      mincol = i;
    }
  }

  // Determine the dominant side
  sign = (max >= -min) ? 1 : -1;

  // Determine new left and right boundary, add 40 per cent...
  // Loop through each column of F
  for (i = 0; i < g_param.ntypes; i++)
  {
    // Column of F being looped through
    j = g_calc.paircol + g_param.ntypes + i;

    // For each column find the bounds of the x-axis of F
    // that just encompass the min and maxrho by 30% of a
    // step, h
    left[i] = minrho[i] - 0.3 * pt->step[j];
    right[i] = maxrho[i] + 0.3 * pt->step[j];

    // Check if expansion is necessary
    // Is minrho outside of range of F's x-axis
    // OR is minrho too far inside the range of F
    // OR same for maxrho
    if (flag || minrho[i] - pt->begin[j] < 0.0 ||
        minrho[i] - pt->begin[j] > 0.95 * pt->step[j] || maxrho[i] - pt->end[j] > 0 ||
        maxrho[i] - pt->end[j] < -0.95 * pt->step[j])
      flag = 1;  // Continue with scaling
  }

  // Determine the scaling factor
  a = (sign == 1) ? upper / right[maxcol] : upper / left[mincol];

  // See if update is still needed
  // Scaling factor will increase it by at least 5% in
  // either direction
  if (flag || fabs(a) > 1.05 || fabs(a) < 0.95)
    flag = 1;  // then continue with scaling

  // If no update needed then die
  if (!flag)
    return 0;

  // BEGIN UPDATING F and other potentials
  ///////////////////////////////////////////////

  // Expand the potential
  h = 0;

  // Loop over columns of F
  for (i = 0; i < g_param.ntypes; ++i)
  {
    // Column of F being looped through
    col = g_calc.paircol + g_param.ntypes + i;

    // Number of splines
    vals = pt->last[col] - pt->first[col];

    // Steps between spline points for new F
    // with bounds right and left
    newstep[i] = (right[i] - left[i]) / (double)vals;

    // Temp variable storing location of leftmost spline knot location
    pos = left[i];

    // Loop through each spline knot
    for (j = 0; j <= vals; ++j)
    {
      // Interpolate or extrapolate value at "pos"
      // Non-equidistant spline points
      newxi[h] = splint_ne(pt, xi, col, pos);

      // Set x-coord for this point in the new potential
      neword[h] = pos;

      // Increment to next potential spot in the array
      ++h;

      // Set the new x-coord position to the next new step point
      pos += newstep[i];
    }  // END OF LOOP OVER SPLINE KNOTS

    // Correct the gradients as well
    // If NOT natural spline (2nd deriv == 0)
    // then reset gradient to extrapolated or interpolated
    // gradient at left and right new bounds
    if (*(xi + pt->first[col] - 2) < 1.e30)
      *(xi + pt->first[col] - 2) = splint_grad_ne(pt, xi, col, left[i]);
    if (*(xi + pt->first[col] - 1) < 1.e30)
      *(xi + pt->first[col] - 1) = splint_grad_ne(pt, xi, col, right[i]);
  }  // END OF LOOP OVER F cols

  // Now that we have the new spline values for F
  // we need to write back the new F to the old potential table
  col = 0;  // Loop over F columns in the potential table
  for (j = g_calc.paircol + g_param.ntypes; j < g_calc.paircol + 2 * g_param.ntypes; ++j)
  {
    // Loop over each spline point in the old table
    for (i = pt->first[j]; i <= pt->last[j]; ++i)
    {
      // Set the value in the old table to the new value
      xi[i] = newxi[col];

      // Set the xcoord to the new value
      pt->xcoord[i] = neword[col];

      // Go to next point in the new F array
      ++col;
    }
  }

  // Scale the actual values of rho_ij to keep physics the same
  // Loop through each rho_ij column for an alloy
  for (i = g_calc.paircol; i < g_calc.paircol + g_param.ntypes; ++i)
  {
    // Loop through each spline point
    for (j = pt->first[i]; j <= pt->last[i]; ++j)
    {
      pt->table[j] *= a;  // Set new value
    }

    // Scale gradient if not natural spline
    if (*(xi + pt->first[i] - 2) < 1.e30)
      *(xi + pt->first[i] - 2) *= a;
    if (*(xi + pt->first[i] - 1) < 1.e30)
      *(xi + pt->first[i] - 1) *= a;
  }

  // In MEAM you have a*f*f*g, where 'a' scale factor is
  // carried over from scaling the total density
  // We multiply this 'a' to the g potential

  // Loop through each g column for an alloy
  for (i = 2 * g_calc.paircol + 2 * g_param.ntypes;
       i < 2 * g_calc.paircol + 3 * g_param.ntypes; ++i)
  {
    // Loop through each spline point
    for (j = pt->first[i]; j <= pt->last[i]; ++j)
    {
      pt->table[j] *= a;  // Set new value
    }

    // Scale gradient if not natural spline
    if (*(xi + pt->first[i] - 2) < 1.e30)
      *(xi + pt->first[i] - 2) *= a;
    if (*(xi + pt->first[i] - 1) < 1.e30)
      *(xi + pt->first[i] - 1) *= a;
  }

  // Output some details
  printf("Scaling factor %f\n", a);

  // Rescele all F columns by 'a' now
  // Be careful, if sign is negative then reverse potential
  if (sign == 1)
  {
    j = 0;

    // Loop through each F column
    for (i = g_calc.paircol + g_param.ntypes; i < g_calc.paircol + 2 * g_param.ntypes;
         ++i)
    {
      // Reset first/last x-coord value in the table to scaled value
      pt->begin[i] = a * left[j];
      pt->end[i] = a * right[j];

      // Rescale the step size and inverse step size
      pt->step[i] = a * newstep[j];
      pt->invstep[i] = 1.0 / pt->step[i];

      // Rescale the gradients as well unless natural spline
      if (xi[pt->first[i] - 2] < 1.e30)
        xi[pt->first[i] - 2] /= a;
      if (xi[pt->first[i] - 1] < 1.e30)
        xi[pt->first[i] - 1] /= a;

      // Set pos to new beginning x-coord value of F
      pos = pt->begin[i];

      // Loop through each spline knot of F
      for (h = pt->first[i]; h <= pt->last[i]; ++h)
      {
        // Set x-coord to new position
        pt->xcoord[h] = pos;

        // Increment position by new step size
        pos += pt->step[i];
      }

      // Go to next pseudo-column in large array
      ++j;
    }  // END OF LOOP THROUGH EACH COL IN F
  }
  else
  {  // Reverse the F potential
    j = 0;

    // Loop through each F column
    for (i = g_calc.paircol + g_param.ntypes; i < g_calc.paircol + 2 * g_param.ntypes;
         ++i)
    {
      // Reset first/last x-coord value in the table to scaled value
      pt->begin[i] = a * right[j];
      pt->end[i] = a * left[j];

      // Rescale the step size and inverse step size
      pt->step[i] = -a * newstep[j];
      pt->invstep[i] = 1.0 / pt->step[i];

      // Rescale the gradients as well unless natural spline
      // Need to store grad of left side since we swap the grads
      // on either side after scaling
      if (xi[pt->first[i] - 2] < 1.e30)
        grad = -xi[pt->first[i] - 2] / a;

      else
        grad = 1.e30;
      if (xi[pt->first[i] - 1] < 1.e30)
        // WARNING: POSSIBLE MISTAKE IN OLDER CODE!!!!
        xi[pt->first[i] - 2] = -xi[pt->first[i] - 1] / a;

      else
        xi[pt->first[i] - 2] = 1.e30;

      // Finish the swap
      xi[pt->first[i] - 1] = grad;

      // Set pos to new beginning x-coord value of F
      pos = pt->begin[i];

      // Loop through each spline knot of F
      for (h = pt->first[i]; h <= pt->last[i]; ++h)
      {
        // Set x-coord to new position
        pt->xcoord[h] = pos;

        // Increment position by new step size
        pos += pt->step[i];
      }

      // Go to next pseudo-column in large array
      ++j;
    }  // END OF LOOP THROUGH EACH COL IN F
    h = 0;

    // Loop through each column in F
    // We will be setting the values in reverse order
    for (i = 0; i < g_param.ntypes; ++i)
    {
      col = g_calc.paircol + g_param.ntypes + i;

      // Loop through each spline point in REVERSE ORDER for this column
      for (j = pt->last[col]; j >= pt->first[col]; --j)
      {
        // Store these values for swapping later
        newxi[h] = xi[j];
        ++h;
      }
    }
    col = 0;

    // Loop through each column in F
    // We will be setting back the values we stored earlier in reverse order
    for (j = g_calc.paircol + g_param.ntypes; j < g_calc.paircol + 2 * g_param.ntypes;
         ++j)
    {
      // Loop through each spline knot in CORRECT ORDER for this column
      for (i = pt->first[j]; i <= pt->last[j]; ++i)
      {
        // Swap back value now in reverse order
        xi[i] = newxi[col];

        // Store the x-coord for this point
        // WARNING: We set this earlier, is this necessary???
        // pt->xcoord[i] = neword[col];
        ++col;
      }
    }
  }  // END OF REVERSAL OF F POTENTIAL

// Only worry about EAM below
#ifdef EAM

  // RE-Initialize the 2nd derivs for splines, so that we can interpolate
  // in the future: only initialize rho, F, and f since they changed
  for (col = g_calc.paircol; col < g_calc.paircol + 2 * g_param.ntypes; ++col)
  {
    // Pointer to first entry
    first = pt->first[col];

    // Get 2nd derivatives
    // step = width of spline knots (known as h)
    // xi+first = array with spline values
    // pt->last[col1] - first + 1 = num of spline pts
    // *(xi + first - 2) = value of endpoint gradient (default: 1e30)
    // *(xi + first - 1) = value of other endpoint gradient
    //          (default: phi=0.0, rho=0.0, F=1e30)
    // g_pot.calc_pot.d2tab + first = array to hold 2nd deriv
    spline_ed(pt->step[col], xi + first, pt->last[col] - first + 1, *(xi + first - 2),
              *(xi + first - 1), pt->d2tab + first);
  }

  // Now we worry about the gauge conditions and set F'(rho-mean) = 0 FOR EAM
  // we set F(rho-mean) = 0 for MEAM
  ///////////////////////////////////////////////////////////////////////////////

  // Loop through each column in F
  for (i = 0; i < g_param.ntypes; ++i)
  {
    // Store gradient for F' at the point halfway between the domain of F
    // at each column in lambda
    lambda[i] = splint_grad(&opt_pot, pt->table, g_calc.paircol + g_param.ntypes + i,
                            0.5 * (pt->begin[g_calc.paircol + g_param.ntypes + i] +
                                   pt->end[g_calc.paircol + g_param.ntypes + i]));
  }

  // Now output the gradient in each column of F at this point
  for (i = 0; i < g_param.ntypes; i++)
    printf("lambda[%d] = %f\n", i, lambda[i]);
  i = 0;

  // Loop through all columns of pair potential, phi
  for (col = 0; col < g_param.ntypes; ++col)
  {
    for (col2 = col; col2 < g_param.ntypes; ++col2)
    {
      // Loop through each spline point for phi
      for (j = pt->first[i]; j <= pt->last[i]; ++j)
      {
        // Gauge transformation to keep F' at halfway point have 0 gradient
        // Phi' = Phi + lambda_i * rho_j + lambda_j * rho_i
        pt->table[j] +=
            (pt->xcoord[j] < pt->end[g_calc.paircol + col2]
                 ? lambda[col] *
                       splint_ne(pt, pt->table, g_calc.paircol + col2, pt->xcoord[j])
                 : 0.0) +
            (pt->xcoord[j] < pt->end[g_calc.paircol + col]
                 ? lambda[col2] *
                       splint_ne(pt, pt->table, g_calc.paircol + col, pt->xcoord[j])
                 : 0.0);
      }  // End of Loop over spline points

      // Now we fix Phi's gradients
      if (pt->table[pt->first[i] - 2] < 1.e30)  // natural spline
        pt->table[pt->first[i] - 2] +=
            (pt->begin[i] < pt->end[g_calc.paircol + col2]
                 ? lambda[col] *
                       splint_grad(pt, pt->table, g_calc.paircol + col2, pt->begin[i])
                 : 0.0) +
            (pt->begin[i] < pt->end[g_calc.paircol + col]
                 ? lambda[col2] *
                       splint_grad(pt, pt->table, g_calc.paircol + col, pt->begin[i])
                 : 0.0);
      if (pt->table[pt->first[i] - 1] < 1.e30)
        pt->table[pt->first[i] - 1] +=
            (pt->end[i] < pt->end[g_calc.paircol + col2]
                 ? lambda[col] *
                       splint_grad(pt, pt->table, g_calc.paircol + col2, pt->end[i])
                 : 0.0) +
            (pt->end[i] < pt->end[g_calc.paircol + col]
                 ? lambda[col2] *
                       splint_grad(pt, pt->table, g_calc.paircol + col, pt->end[i])
                 : 0.0);
      ++i;
    }  // END OF INNER LOOP OVER PHI
  }    // END OF OUTER LOOP OVER PHI

  // Update F potential using gauge condition
  for (i = 0; i < g_param.ntypes; i++)
  {
    // Loop through spline points of F
    // Use gauge condition F'_i = F_i - lambda_i * rho_i
    for (j = pt->first[g_calc.paircol + g_param.ntypes + i];
         j <= pt->last[g_calc.paircol + g_param.ntypes + i]; ++j)
      pt->table[j] -= pt->xcoord[j] * lambda[i];

    // Fix the gradients as well
    if (pt->table[pt->first[g_calc.paircol + g_param.ntypes + i] - 2] <
        1.e30)  // natural spline
      pt->table[pt->first[g_calc.paircol + g_param.ntypes + i] - 2] -= lambda[i];
    if (pt->table[pt->first[g_calc.paircol + g_param.ntypes + i] - 1] < 1.e30)
      pt->table[pt->first[g_calc.paircol + g_param.ntypes + i] - 1] -= lambda[i];
    lambda[i] = 0.0;
  }

  // Initialize 2nd derivatives just for pair potential (not sure why though)
  for (col = 0; col < g_calc.paircol; ++col)
  {
    first = pt->first[col];

    spline_ed(pt->step[col], pt->table + first, pt->last[col] - first + 1,
              *(pt->table + first - 2), 0.0, pt->d2tab + first);
  }

#endif /* EAM */

  free_local_memory();

  return a;
}

#endif /* MEAM */
#endif /* RESCALE || !APOT */
