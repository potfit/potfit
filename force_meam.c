/****************************************************************
*
* force_meam.c: Routine used for calculating meam forces/energies
* 	in various interpolation schemes.
*
*****************************************************************/
/*
*   Copyright 2002-2010 Peter Brommer, Franz G"ahler, Daniel Schopf
*             Institute for Theoretical and Applied Physics
*             University of Stuttgart, D-70550 Stuttgart, Germany
*             http://www.itap.physik.uni-stuttgart.de/
*	      MEAM-potential: Jeremy Nicklas
*                  Ohio State University, Physics Department
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
* $Revision: 1.x $
* $Date: 2010/07/15 12:24:43 $
*****************************************************************/

#if defined MEAM
#include "potfit.h"

/****************************************************************
 *
 *  compute forces using eam potentials with spline interpolation
 *
 *  returns sum of squares of differences between calculated and reference
 *     values
 *
 *  arguments: *xi - pointer to potential
 *             *forces - pointer to forces calculated from potential
 *             flag - used for special tasks
 *
 * When using the mpi-parallelized version of potfit, all processes but the
 * root process jump into this function immediately after initialization and
 * stay in here for an infinite loop, to exit only when a certain flag value
 * is passed from process 0. When a set of forces needs to be calculated,
 * the root process enters the function with a flag value of 0, broadcasts
 * the current potential table xi and the flag value to the other processes,
 * thus initiating a force calculation. Whereas the root process returns with
 * the result, the other processes stay in the loop. If the root process is
 * called with flag value 1, all processes exit the function without
 * calculating the forces.
 * If anything changes about the potential beyond the values of the parameters,
 * e.g. the location of the sampling points, these changes have to be broadcast
 * from rank 0 process to the higher ranked processes. This is done when the
 * root process is called with flag value 2. Then a potsync function call is
 * initiated by all processes to get the new potential from root.
 *
 * xi_opt is the array storing the potential parameters (usually it is the
 *     opt_pot.table - part of the struct opt_pot, but it can also be
 *     modified from the current potential.
 *
 * forces is the array storing the deviations from the reference data, not
 *     only for forces, but also for energies, stresses or dummy constraints
 *     (if applicable).
 *
 * flag is an integer controlling the behaviour of calc_forces_eam.
 *    flag == 1 will cause all processes to exit calc_forces_eam after
 *             calculation of forces.
 *    flag == 2 will cause all processes to perform a potsync (i.e. broadcast
 *             any changed potential parameters from process 0 to the others)
 *             before calculation of forces
 *    all other values will cause a set of forces to be calculated. The root
 *             process will return with the sum of squares of the forces,
 *             while all other processes remain in the function, waiting for
 *             the next communication initiating another force calculation
 *             loop
 *
 ****************************************************************/

// DO NOT USE:
// APOT, PARABEL, MPI, OPENMP, WZERO

double calc_forces_meam(double *xi_opt, double *forces, int flag)
{
  // Set temporary pot table
  double *xi = xi_opt;

  // Some useful temp variables
  int   h, col, first;
  double rho_sum_loc, tmpsum;
  double rho_sum = 0., sum = 0.;

  // Only works for format=3 for now
  if (format != 3)
    error(1, "Needs format = 3.\n");

  while (1) {

    // Reset tmpsum and rho_sum_loc
    // tmpsum = Sum of all the forces, energies and constraints
    // rho_sum_loc = Sum of density, rho, for all atoms
    tmpsum = 0.;
    rho_sum_loc = 0.;

#ifdef MPI
    // Broadcast the potential and flag to all the procs
    // The flag can be used to kill the slave procs
    MPI_Bcast(xi, calc_pot.len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Potentials have been rescaled so broadcast them
    if (flag == 2)
      potsync();
    // Kill slave procs
    if (flag == 1)
      break;
#endif // MPI

    // First step is to initialize 2nd derivatives for splines

    // Pair potential (phi), density (rho), embedding funtion (F)
    // where paircol is number of pair potential columns
    // and ntypes is number of rho columns
    // and ntypes is number of F columns
    for (col = 0; col < 2 * paircol + 3 * ntypes; ++col) {
      // Pointer to first entry
      first = calc_pot.first[col];

      // Get 2nd derivatives
      // step = width of spline knots (known as h)
      // xi+first = array with spline values
      // calc_pot.last[col1] - first + 1 = num of spline pts
      // *(xi + first - 2) = value of endpoint gradient (default: 1e30)
      // *(xi + first - 1) = value of other endpoint gradient
      //          (default: phi=0.0, rho=0.0, F=1e30)
      // calc_pot.d2tab + first = array to hold 2nd deriv
      spline_ed(calc_pot.step[col], xi + first, calc_pot.last[col] - first + 1,
	*(xi + first - 2), *(xi + first - 1), calc_pot.d2tab + first);
    }

#ifndef MPI
    myconf = nconf;
#endif // MPI

    {
      // Temp variables
      atom_t *atom;		// atom type
      int   i, j, k, l;

#ifdef STRESS
      int   stresses;
#endif // STRESS

      // Some useful temp struct variable types
      neigh_t *neigh;		// neighbor type

      // Pair variables
      double phi_val, phi_grad;
      vector tmp_force;		// x,y,z vector

      // Eam variables
      int   col_F;
      double eam_force;
#ifdef NORESCALE
      double rho_val;
#endif

      // Meam variables
      int   m, jj, kk, ijk;
      double dV3j, dV3k, V3, vlj, vlk, vv3j, vv3k, dfj[3], dfk[3];
      neigh_t *neigh_j, *neigh_k;
      angl *n_angl;		// angle type

      // Loop over configurations
      for (h = firstconf; h < firstconf + myconf; h++) {

	// Reset energies
	forces[energy_p + h] = 0.;

#ifdef STRESS
	// Reset stresses
	for (i = 0; i < 6; ++i)
	  forces[stress_p + 6 * h + i] = 0.;
#endif // STRESS

#ifndef NORESCALE
	// Set limiting constraints
	forces[limit_p + h] = -force_0[limit_p + h];
#endif // NORESCALE

	// FIRST LOOP: Reset forces and densities for each atom
	for (i = 0; i < inconf[h]; ++i) {

	  // Skip every 3 spots in force array starting from
	  // position of first atom
	  k = 3 * (cnfstart[h] + i);

	  // Set initial forces to negative of user given forces
	  // so we can take difference
	  forces[k] = -force_0[k];
	  forces[k + 1] = -force_0[k + 1];
	  forces[k + 2] = -force_0[k + 2];

	  // Reset the density for each atom
	  conf_atoms[cnfstart[h] - firstatom + i].rho = 0.0;

	}			// END OF FIRST LOOP


	// SECOND LOOP: Calculate pair forces and energies, atomic densities
	for (i = 0; i < inconf[h]; ++i) {

	  // Set pointer to temp atom pointer
	  atom = conf_atoms + (cnfstart[h] - firstatom + i);

	  // Skip every 3 spots for force array
	  k = 3 * (cnfstart[h] + i);

	  // Loop over neighbors
	  for (j = 0; j < atom->n_neigh; ++j) {

	    // Set pointer to temp neighbor pointer
	    neigh = atom->neigh + j;

	    // Find the correct column in the potential table for pair potential: phi_ij
	    // For Binary Alloy: 0 = phi_AA, 1 = phi_AB, phi_BA, 2 = phi_BB
	    // where typ = A = 0 and typ = B = 1

	    // We need to check that neighbor atom exists inside pair
	    // potential's radius
	    if (neigh->r < calc_pot.end[neigh->col[0]]) {
	      // Compute phi and phi' value given radial distance
	      // NOTE: slot = spline point index right below radial distance
	      // shift = % distance from 'slot' spline pt
	      // step = width of spline points (given as 'h' in books)
	      // 0 means the pair potential spline info???
	      phi_val =
		0.5 * splint_comb_dir(&calc_pot, xi, neigh->slot[0],
		neigh->shift[0], neigh->step[0], &phi_grad);

	      // Only half of the gradient contributes to the force
	      // as well as half of the energy since we are double counting
	      phi_grad *= 0.5;

	      // Add in piece contributed by neighbor to energy
	      forces[energy_p + h] += phi_val;

	      // Compute tmp force values
	      tmp_force.x = neigh->dist.x * phi_grad;
	      tmp_force.y = neigh->dist.y * phi_grad;
	      tmp_force.z = neigh->dist.z * phi_grad;

	      // Add in force on atom i from atom j
	      forces[k] += tmp_force.x;
	      forces[k + 1] += tmp_force.y;
	      forces[k + 2] += tmp_force.z;

	      // Subtract off force on atom j from atom i
	      // Newton's law: action = -reaction
	      l = 3 * neigh->nr;
	      forces[l] -= tmp_force.x;
	      forces[l + 1] -= tmp_force.y;
	      forces[l + 2] -= tmp_force.z;

#ifdef STRESS
	      tmp_force.x *= neigh->r;
	      tmp_force.y *= neigh->r;
	      tmp_force.z *= neigh->r;
	      stresses = stress_p + 6 * h;
	      forces[stresses] -= neigh->dist.x * tmp_force.x;
	      forces[stresses + 1] -= neigh->dist.y * tmp_force.y;
	      forces[stresses + 2] -= neigh->dist.z * tmp_force.z;
	      forces[stresses + 3] -= neigh->dist.x * tmp_force.y;
	      forces[stresses + 4] -= neigh->dist.y * tmp_force.z;
	      forces[stresses + 5] -= neigh->dist.z * tmp_force.x;
#endif // STRESS

	    }			// END IF STMNT: NEIGH LIES INSIDE CUTOFF FOR PAIR POTENTIAL

	    // Find the correct column in the potential table for atomic density, rho_ij
	    // paircol = number of pair potential columns
	    // Binary Alloy: paircol = 2... 2 = rho_A, 3 = rho_B
	    // where A, B are atom type for the neighbor

	    // Compute rho rho value and sum them up
	    // Need to play tricks so that rho values are put in the correct
	    // columns if alloy. If atom j is A or B, fn value needs to be
	    // in correct rho_A or rho_B respectively, it doesn't depend
	    // on atom i.

	    // Check that atom j lies inside rho_typ2
	    if (neigh->r < calc_pot.end[neigh->col[1]]) {

	      // Store gradient in the neighbor for the pair r_ij
	      // to be used in the future when computing forces
	      // and sum up rho for atom i
	      atom->rho +=
		splint_comb_dir(&calc_pot, xi, neigh->slot[1], neigh->shift[1], neigh->step[1], &neigh->drho);
	    } else {

	      // If the pair distance does not lie inside rho_typ2
	      // We set the grad to 0 so it doesn't sum into the net
	      // force
	      neigh->drho = 0.;
	    }

	    // Compute the f_ij values and store the fn and grad
	    // in each neighbor struct for easy access later
	    //////////////////////////////////////////////////////////

	    // Find the correct column in the potential table for "f": f_ij
	    // For Binary Alloy: 0 = f_AA, 1 = f_AB, f_BA, 2 = f_BB
	    // where typ = A = 0 and typ = B = 1
	    // Note: it is "paircol+2*ntypes" spots away in the array

	    // Check that atom j lies inside f_col2
	    if (neigh->r < calc_pot.end[neigh->col[2]]) {

	      // Store the f(r_ij) value and the gradient for future use
	      neigh->f =
		splint_comb_dir(&calc_pot, xi, neigh->slot[2], neigh->shift[2], neigh->step[2], &neigh->df);
	    } else {

	      // Store f and f' = 0 if doesn't lie in boundary
	      // to be used later when calculating forces
	      neigh->f = 0.;
	      neigh->df = 0.;
	    }

	  }			// END LOOP OVER NEIGHBORS

	  // Find the correct column in the potential table for angle part: g_ijk
	  // Binary Alloy: 0 = g_A, 1 = g_B
	  // where A, B are atom type for the main atom i
	  // Note: it is now "2*paircol+2*ntypes" from beginning column
	  // to account for phi(paircol)+rho(nytpes)+F(ntypes)+f(paircol)
	  //col2 = 2 * paircol + 2 * ntypes + typ1;


	  // Loop over every angle formed by neighbors
	  // N(N-1)/2 possible combinations
	  // Used in computing angular part g_ijk
	  ijk = 0;		// count number of angles
	  for (jj = 0; jj < atom->n_neigh - 1; ++jj) {

	    // Get pointer to neighbor jj
	    neigh_j = atom->neigh + jj;

	    for (kk = jj + 1; kk < atom->n_neigh; ++kk) {

	      // Store pointer to angular part (g)
	      n_angl = atom->angl_part + ijk;

	      // Get pointer to neighbor kk
	      neigh_k = atom->neigh + kk;

	      // The cos(theta) should always lie inside -1 ... 1
	      // So store the g and g' without checking bounds
	      n_angl->g =
		splint_comb_dir(&calc_pot, xi, n_angl->slot, n_angl->shift, n_angl->step, &n_angl->dg);

	      // Sum up rho piece for atom i caused by j and k
	      // f_ij * f_ik * m_ijk
	      atom->rho += neigh_j->f * neigh_k->f * n_angl->g;

	      ++ijk;
	    }
	  }

	  // Column for embedding function, F
	  col_F = paircol + ntypes + atom->typ;

#ifndef NORESCALE
	  // Compute energy, gradient for embedding function F
	  // Check if rho lies short of inner cutoff of F(rho)
	  if (atom->rho < calc_pot.begin[col_F]) {

	    // Punish this potential for having rho lie outside of F
	    forces[limit_p + h] += DUMMY_WEIGHT * 10. * dsquare(calc_pot.begin[col_F] - atom->rho);

	    // Set the atomic density to the first rho in the spline F
	    atom->rho = calc_pot.begin[col_F];

	  } else if (atom->rho > calc_pot.end[col_F]) {	// rho is to the right of the spline

	    // Punish this potential for having rho lie outside of F
	    forces[limit_p + h] += DUMMY_WEIGHT * 10. * dsquare(atom->rho - calc_pot.end[col_F]);

	    // Set the atomic density to the last rho in the spline F
	    atom->rho = calc_pot.end[col_F];
	  }
	  // Compute energy piece from F, and store the gradient for later use
	  forces[energy_p + h] += splint_comb(&calc_pot, xi, col_F, atom->rho, &atom->gradF);

#elif defined NORESCALE		// !NORESCALE
	  // Compute energy, gradient for embedding function F
	  // Check if rho lies short of inner cutoff of F(rho)
	  if (atom->rho < calc_pot.begin[col_F]) {

	    // Linear extrapolate values to left to get F_i(rho)
	    // This gets value and grad of initial spline point
	    rho_val = splint_comb(&calc_pot, xi, col_F, calc_pot.begin[col_F], &atom->gradF);

	    // Sum this to the total energy for this configuration
	    // Linear extrapolate this energy
	    forces[energy_p + h] += rho_val + (atom->rho - calc_pot.begin[col_F]) * atom->gradF;

	  } else if (atom->rho > calc_pot.end[col_F]) {	// rho is to the right of the spline

	    // Get value and grad at 1/2 the width from the final spline point
	    rho_val =
	      splint_comb(&calc_pot, xi, col_F,
	      calc_pot.end[col_F] - .5 * calc_pot.step[col_F], &atom->gradF);

	    // Linear extrapolate to the right to get energy
	    forces[energy_p + h] += rho_val + (atom->rho - calc_pot.end[col_F]) * atom->gradF;
	  } else {

	    // Get energy value from within spline and store the grad
	    forces[energy_p + h] += splint_comb(&calc_pot, xi, col_F, atom->rho, &atom->gradF);
	  }
#endif // !NORESCALE

	  // Sum up rho for future MPI use
	  rho_sum_loc += atom->rho;



	  ///////////////////////////////////////////////////////////
	  // Calculate remaining forces from embedding function
	  ///////////////////////////////////////////////////////////

	  // Loop over neighbors
	  for (j = 0; j < atom->n_neigh; ++j) {

	    // Set pointer to temp neighbor pointer and record type
	    neigh = atom->neigh + j;

	    // Check that radial distance between pair is within
	    // cutoff distance of either possible rho_A or rho_B
	    // for alloys, where A or B stands for atom i
	    // WARNING: Double check this!!! May not need this
	    // since drho will be 0 otherwise
	    if (neigh->r < calc_pot.end[neigh->col[1]]) {

	      // Calculate eam force
	      eam_force = neigh->drho * atom->gradF;

	      // Multiply the eamforce with x/r to get real force
	      tmp_force.x = neigh->dist.x * eam_force;
	      tmp_force.y = neigh->dist.y * eam_force;
	      tmp_force.z = neigh->dist.z * eam_force;

	      // Sum up forces acting on atom i from atom j
	      forces[k] += tmp_force.x;
	      forces[k + 1] += tmp_force.y;
	      forces[k + 2] += tmp_force.z;

	      // Subtract off forces acting on atom j from atom i
	      l = 3 * neigh->nr;
	      forces[l] -= tmp_force.x;
	      forces[l + 1] -= tmp_force.y;
	      forces[l + 2] -= tmp_force.z;

#ifdef STRESS
	      tmp_force.x *= neigh->r;
	      tmp_force.y *= neigh->r;
	      tmp_force.z *= neigh->r;
	      stresses = stress_p + 6 * h;
	      forces[stresses] -= neigh->dist.x * tmp_force.x;
	      forces[stresses + 1] -= neigh->dist.y * tmp_force.y;
	      forces[stresses + 2] -= neigh->dist.z * tmp_force.z;
	      forces[stresses + 3] -= neigh->dist.x * tmp_force.y;
	      forces[stresses + 4] -= neigh->dist.y * tmp_force.z;
	      forces[stresses + 5] -= neigh->dist.z * tmp_force.x;
#endif // STRESS

	    }			// END IF STMT: Inside reach of rho cutoff

	  }			// END LOOP OVER NEIGHBORS


	  // Compute MEAM Forces
	  //////////////////////////////////

	  // Loop over every angle formed by neighbors
	  // N(N-1)/2 possible combinations
	  // Used in computing angular part g_ijk
	  ijk = 0;		// count number of angles
	  for (jj = 0; jj < atom->n_neigh - 1; ++jj) {

	    // Get pointer to neighbor j
	    neigh_j = atom->neigh + jj;

	    // Force location for atom j
	    l = 3 * neigh_j->nr;

	    for (kk = jj + 1; kk < atom->n_neigh; ++kk) {

	      // Store pointer to angular part (g)
	      n_angl = atom->angl_part + ijk;

	      // Get pointer to neighbor k
	      neigh_k = atom->neigh + kk;

	      // Force location for atom k
	      m = 3 * neigh_k->nr;

	      // Some tmp variables to clean up force fn below
	      dV3j = n_angl->g * neigh_j->df * neigh_k->f;
	      dV3k = n_angl->g * neigh_j->f * neigh_k->df;
	      V3 = neigh_j->f * neigh_k->f * n_angl->dg;

	      vlj = V3 * neigh_j->recip;
	      vlk = V3 * neigh_k->recip;
	      vv3j = dV3j - vlj * n_angl->cos;
	      vv3k = dV3k - vlk * n_angl->cos;

	      dfj[0] = vv3j * neigh_j->dist.x + vlj * neigh_k->dist.x;
	      dfj[1] = vv3j * neigh_j->dist.y + vlj * neigh_k->dist.y;
	      dfj[2] = vv3j * neigh_j->dist.z + vlj * neigh_k->dist.z;

	      dfk[0] = vv3k * neigh_k->dist.x + vlk * neigh_j->dist.x;
	      dfk[1] = vv3k * neigh_k->dist.y + vlk * neigh_j->dist.y;
	      dfk[2] = vv3k * neigh_k->dist.z + vlk * neigh_j->dist.z;

	      // Force on atom i from j and k
	      forces[k] += atom->gradF * (dfj[0] + dfk[0]);
	      forces[k + 1] += atom->gradF * (dfj[1] + dfk[1]);
	      forces[k + 2] += atom->gradF * (dfj[2] + dfk[2]);

	      // Reaction force on atom j from i and k
	      forces[l] -= atom->gradF * dfj[0];
	      forces[l + 1] -= atom->gradF * dfj[1];
	      forces[l + 2] -= atom->gradF * dfj[2];

	      // Reaction force on atom k from i and j
	      forces[m] -= atom->gradF * dfk[0];
	      forces[m + 1] -= atom->gradF * dfk[1];
	      forces[m + 2] -= atom->gradF * dfk[2];

#ifdef STRESS
	      // Force from j on atom i (guessing here)
	      tmp_force.x = atom->gradF * dfj[0] * neigh_j->r;
	      tmp_force.y = atom->gradF * dfj[1] * neigh_j->r;
	      tmp_force.z = atom->gradF * dfj[2] * neigh_j->r;
	      stresses = stress_p + 6 * h;
	      forces[stresses] -= neigh_j->dist.x * tmp_force.x;
	      forces[stresses + 1] -= neigh_j->dist.y * tmp_force.y;
	      forces[stresses + 2] -= neigh_j->dist.z * tmp_force.z;
	      forces[stresses + 3] -= neigh_j->dist.x * tmp_force.y;
	      forces[stresses + 4] -= neigh_j->dist.y * tmp_force.z;
	      forces[stresses + 5] -= neigh_j->dist.z * tmp_force.x;

	      // Force from k on atom i (also guessing here)
	      tmp_force.x = atom->gradF * dfk[0] * neigh_k->r;
	      tmp_force.y = atom->gradF * dfk[1] * neigh_k->r;
	      tmp_force.z = atom->gradF * dfk[2] * neigh_k->r;
	      stresses = stress_p + 6 * h;
	      forces[stresses] -= neigh_k->dist.x * tmp_force.x;
	      forces[stresses + 1] -= neigh_k->dist.y * tmp_force.y;
	      forces[stresses + 2] -= neigh_k->dist.z * tmp_force.z;
	      forces[stresses + 3] -= neigh_k->dist.x * tmp_force.y;
	      forces[stresses + 4] -= neigh_k->dist.y * tmp_force.z;
	      forces[stresses + 5] -= neigh_k->dist.z * tmp_force.x;
#endif // STRESS

	      ++ijk;
	    }			// End inner loop over angles (neighbor atom kk)
	  }			// End outer loop over angles (neighbor atom jj)
	}			// END OF SECOND LOOP OVER ATOM i

	// 3RD LOOP OVER ATOM i
	// Sum up the square of the forces for each atom
	// then multiply it by the weight for this config
	for (i = 0; i < inconf[h]; ++i) {

	  // Skip every 3 spots for force array
	  k = 3 * (cnfstart[h] + i);
	  forces[k] *= conf_weight[h];
	  forces[k + 1] *= conf_weight[h];
	  forces[k + 2] *= conf_weight[h];
	  tmpsum += (dsquare(forces[k]) + dsquare(forces[k + 1]) + dsquare(forces[k + 2]));

	}			// END OF THIRD LOOP OVER ATOM i

	// Add in the energy per atom and its weight to the sum
	// First divide by num atoms
	forces[energy_p + h] /= (double)inconf[h];

	// Then subtract off the cohesive energy given to use by user
	forces[energy_p + h] -= force_0[energy_p + h];

	// Sum up square of this new energy term for each config
	// multiplied by its respective weight
	forces[energy_p + h] *= conf_weight[h] * eweight;
	tmpsum += dsquare(forces[energy_p + h]);

#ifdef STRESS
	// LOOP OVER STRESSES
	for (i = 0; i < 6; ++i) {

	  // Multiply weight to stresses and divide by volume
	  forces[stress_p + 6 * h + i] /= conf_vol[h - firstconf];

	  // Subtract off user supplied stresses
	  forces[stress_p + 6 * h + i] -= force_0[stress_p + 6 * h + i];

	  // Sum in the square of each stress component with config weight
	  forces[stress_p + 6 * h + i] *= conf_weight[h] * sweight;
	  tmpsum += dsquare(forces[stress_p + 6 * h + i]);

	}			// END LOOP OVER 6 STRESSES
#endif // STRESS

#ifndef NORESCALE
	// Add in the square of the limiting constraints for each config
	// This is punishment from going out of bounds for F(rho)
	// if NORESCALE is not defined
	forces[limit_p + h] *= conf_weight[h];
	tmpsum += dsquare(forces[limit_p + h]);
#endif // not NORESCALE

      }				// END MAIN LOOP OVER CONFIGURATIONS

    }

#ifdef MPI
    // Reduce the rho_sum into root node
    MPI_Reduce(&rho_sum_loc, &rho_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    rho_sum = rho_sum_loc;
#endif // MPI

    // Root node
    if (myid == 0) {

#ifdef NORESCALE
      // Calculate the average rho_sum per atom
      // NOTE: This gauge constraint exists for both EAM and MEAM
      rho_sum /= (double)natoms;

      // Another constraint for the gauge conditions
      // this sets the avg rho per atom to 1
      // Please read the other constraint on gauge conditions
      // above.
      forces[dummy_p + ntypes] = DUMMY_WEIGHT * (rho_sum - 1.);
      tmpsum += dsquare(forces[dummy_p + ntypes]);
#endif // NORESCALE
    }				// END ROOT NODE PROCESSING OF CONSTRAINTS ON F and F'

    // Set tmpsum to sum - only matters when not running MPI
    sum = tmpsum;

#ifdef MPI
    // Reduce the global sum from all the tmpsum's
    sum = 0.;
    MPI_Reduce(&tmpsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // gather forces, energies, stresses
    // forces
    MPI_Gatherv(forces + firstatom * 3, myatoms, MPI_VECTOR, forces, atom_len,
      atom_dist, MPI_VECTOR, 0, MPI_COMM_WORLD);
    // energies
    MPI_Gatherv(forces + natoms * 3 + firstconf, myconf, MPI_DOUBLE,
      forces + natoms * 3, conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // stresses
    MPI_Gatherv(forces + natoms * 3 + nconf + 6 * firstconf, myconf, MPI_STENS,
      forces + natoms * 3 + nconf, conf_len, conf_dist, MPI_STENS, 0, MPI_COMM_WORLD);
    // punishment constraints for going out of bounds for F(rho)
    // only really used when rescaling
    MPI_Gatherv(forces + natoms * 3 + 7 * nconf + firstconf, myconf, MPI_DOUBLE,
      forces + natoms * 3 + 7 * nconf, conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // no need to pick up dummy constraints - are already @ root
    // they are constraints on F and F'
#endif // MPI

    // Root process only
    if (myid == 0) {

      // Increment function calls
      ++fcalls;

      // If total sum is NAN return large number instead
      if (isnan(sum))
	return 10e10;
      else
	return sum;
    }

  }				// END OF INFINITE LOOP

  // Kill off other procs
  return -1;
}

#endif /* MEAM */
