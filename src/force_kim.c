/****************************************************************
 *
 * force_kim.c: Routines used for calculating forces/energies via KIM
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
 ****************************************************************/

#if !defined(KIM)
#error force_kim.c compiled without KIM support
#endif

#include "potfit.h"

#include "force.h"
#include "kim.h"
#include "utils.h"

/****************************************************************
  init_force
****************************************************************/

void init_force(int is_worker)
{
  // nothing to do here for KIM potentials
}

/****************************************************************
  update_kim_model
****************************************************************/

void update_kim_model(double *xi_opt)
{
  int res = 0;

  // set new parameters

  int pos = 0;

  for (int i = 0; i < g_kim.nparams; ++i) {
    for (int j = 0; j < g_kim.params[i].extent; ++j) {
      if (xi_opt[pos] != g_kim.params[i].values[j].d) {
        res = KIM_Model_SetParameterDouble(g_kim.model, i, j, xi_opt[pos]);
        if (res)
          error(1, "Error updating KIM model parameter (%d. %d): %d\n", i, j, res);
      }
      ++pos;
    }
  }

  // refresh the model - should recalculate internal data structures if necessary

  res = KIM_Model_ClearThenRefresh(g_kim.model);
  if (res)
    error(1, "Error calling model refresh: %d\n", res);

  // check if cutoff has changed

  int num_lists = 0;
  const double* cutoffs = NULL;
  const int* data = NULL; // modelWillNotRequestNeighborsOfNoncontributingParticles
  double max_cutoff = 0.0;

  KIM_Model_GetInfluenceDistance(g_kim.model, &max_cutoff);

  KIM_Model_GetNeighborListPointers(g_kim.model, &num_lists, &cutoffs, &data);

  if (g_kim.cutoffs[num_lists] != max_cutoff) {
    if (max_cutoff < g_kim.cutoffs[num_lists]) {
      warning("KIM influence distance has decreased from %.10f to %.10f!\n", g_kim.cutoffs[num_lists], max_cutoff);
      g_kim.cutoffs[num_lists] = max_cutoff;
    } else
      error(1, "KIM influence distance has increased from %.10f to %.10f\n", g_kim.cutoffs[num_lists], max_cutoff);
  }

  for (int i = 0; i < num_lists; ++i) {
    if (!data[0])
      error(1, "KIM model does request neighbors of ghost atoms - NYI in potfit\n");
    if (cutoffs[i] > max_cutoff)
      max_cutoff = cutoffs[i];
    if (g_kim.cutoffs[i] != cutoffs[i]) {
      if (cutoffs[i] < g_kim.cutoffs[i]) {
        warning("KIM cutoff distance %d has decreased from %.10f to %.10f!\n", i, g_kim.cutoffs[i], cutoffs[i]);
        g_kim.cutoffs[i] = cutoffs[i];
      } else
        error(1, "KIM cutoff distance %d has increased from %.10f to %.10f!\n", i, g_kim.cutoffs[i], cutoffs[i]);
    }
  }
}

/****************************************************************
 *
 *  compute forces via KIM
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
 *     g_pot.opt_pot.table - part of the struct g_opt.opt_pot, but it can also be
 *     modified from the current potential.
 *
 * forces is the array storing the deviations from the reference data, not
 *     only for forces, but also for energies, stresses or dummy constraints
 *     (if applicable).
 *
 * flag is an integer controlling the behaviour of calc_forces_pair.
 *    flag == 1 will cause all processes to exit calc_forces_pair after
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

double calc_forces(double *xi_opt, double *forces, int flag)
{
  // This is the start of an infinite loop
  while (1) {
    double error_sum = 0.0; // sum of squares of local process

#if !defined(MPI)
    g_mpi.myconf = g_config.nconf;
#endif // MPI

    update_kim_model(xi_opt);

    // loop over configurations
    for (int h = g_mpi.firstconf; h < g_mpi.firstconf + g_mpi.myconf; h++) {
      double kim_energy = 0.0;
      int ret = KIM_ComputeArguments_SetArgumentPointerDouble(g_kim.arguments[h], KIM_COMPUTE_ARGUMENT_NAME_partialEnergy, &kim_energy);
      if (ret)
        error(1, "Error setting energy pointer\n");

      double* kim_forces = g_kim.helpers[h].forces;
      int uf = g_config.conf_uf[h - g_mpi.firstconf];
      // forces pointer has already been set when setting up the compute argument

#if defined(STRESS)
      double virial[6] = {0};
      int us = g_config.conf_us[h - g_mpi.firstconf];
      if (us) {
        int ret = KIM_ComputeArguments_SetArgumentPointerDouble(g_kim.arguments[h], KIM_COMPUTE_ARGUMENT_NAME_partialVirial, virial);
        if (ret)
          error(1, "Error setting virial pointer\n");
      }
#endif // STRESS

      // Calculate forces from KIM (general forces, including forces, virial and energy)

      ret = KIM_Model_Compute(g_kim.model, g_kim.arguments[h]);
      if (ret)
        error(1, "Error calling KIM_Model_Compute\n");

      // forces contributation
      double weight = sqrt(g_config.conf_weight[h]);

      // gather all KIM mirror image contributions

      for (int i = 0; i < g_config.number_of_particles[h]; ++i) {
        int idx = 3 * g_config.source_atom[h][i];
        forces[idx + 0] += kim_forces[3 * i + 0];
        forces[idx + 1] += kim_forces[3 * i + 1];
        forces[idx + 2] += kim_forces[3 * i + 2];
      }

      for (int i = 0; i < g_config.inconf[h]; i++) {

        int n_i = 3 * (g_config.cnfstart[h] + i);
        if (uf) {
          forces[n_i + 0] = weight * (forces[n_i + 0] - g_config.force_0[n_i + 0]);
          forces[n_i + 1] = weight * (forces[n_i + 1] - g_config.force_0[n_i + 1]);
          forces[n_i + 2] = weight * (forces[n_i + 2] - g_config.force_0[n_i + 2]);
        } else {
          forces[n_i + 0] = 0.0;
          forces[n_i + 1] = 0.0;
          forces[n_i + 2] = 0.0;
        }

#if defined(FWEIGHT) || defined(CONTRIB)
        atom_t* atom = g_config.conf_atoms + i + g_config.cnfstart[h] - g_mpi.firstatom;
#endif // FWEIGHT || CONTRIB
#if defined(FWEIGHT)
        // weigh by absolute value of force
        forces[n_i + 0] /= FORCE_EPS + atom->absforce;
        forces[n_i + 1] /= FORCE_EPS + atom->absforce;
        forces[n_i + 2] /= FORCE_EPS + atom->absforce;
#endif // FWEIGHT
#if defined(CONTRIB)
        if (atom->contrib)
#endif // CONTRIB
          error_sum += (dsquare(forces[n_i + 0]) + dsquare(forces[n_i + 1]) + dsquare(forces[n_i + 2]));
      }

      // reset energies and stresses
      forces[g_calc.energy_p + h] = 0.0;
#if defined(STRESS)
      int stresses_offset = g_calc.stress_p + 6 * h;
      for (int i = 0; i < 6; i++)
        forces[stresses_offset + i] = 0.0;
#endif /* STRESS */

      weight = sqrt(g_config.conf_weight[h] * g_param.eweight);
      forces[g_calc.energy_p + h] = weight * kim_energy;
      forces[g_calc.energy_p + h] -=  weight * g_config.force_0[g_calc.energy_p + h];
      error_sum += dsquare(forces[g_calc.energy_p + h]);

#if defined(STRESS)
      // stress contributions
      if (uf && us) {
        weight = sqrt(g_config.conf_weight[h] * g_param.sweight);
        for (int i = 0; i < 6; i++) {
          forces[stresses_offset + i]  = weight * virial[i];
          forces[stresses_offset + i] /= g_config.conf_vol[h - g_mpi.firstconf];
          forces[stresses_offset + i] -= weight*g_config.force_0[stresses_offset + i];
          error_sum += dsquare(forces[stresses_offset + i]);
        }
      }
#endif // STRESS

    } // loop over configurations

    // add punishment for out of bounds (mostly for powell_lsq)
    if (g_mpi.myid == 0) {
       // loop over individual parameters
      for (int i = 0; i < g_calc.ndim; i++) {
        double min = g_pot.apot_table.pmin[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]];
        double max = g_pot.apot_table.pmax[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]];
        // punishment for out of bounds
        if (xi_opt[g_pot.opt_pot.idx[i]] < min) {
          double x = xi_opt[g_pot.opt_pot.idx[i]] - min;
          error_sum += APOT_PUNISH * x * x;
          forces[g_calc.punish_par_p + i] = APOT_PUNISH * x * x;
        } else if (xi_opt[g_pot.opt_pot.idx[i]] > max) {
          double x = xi_opt[g_pot.opt_pot.idx[i]] - max;
          error_sum += APOT_PUNISH * x * x;
          forces[g_calc.punish_par_p + i] = APOT_PUNISH * x * x;
        }
      }
    }

    gather_forces(&error_sum, forces);

    // root process exits this function now
    if (g_mpi.myid == 0) {
      g_calc.fcalls++;     // Increase function call counter
      if (isnan(error_sum)) {
#if defined(DEBUG)
        printf("\n--> Force is nan! <--\n\n");
#endif // DEBUG
        return 10e10;
      } else
        return error_sum;
    }
  } // end of infinite loop

  // once a non-root process arrives here, all is done
  return -1.0;
}

int get_neigh(const void* const puser,
              const int numberOfNeighborLists,
              const double* const cutoffs,
              const int neighborListIndex,
              const int particleNumber,
              int* const numberOfNeighbors,
              const int** const neighborsOfParticle)
{
  potfit_compute_helper_t* helper = (potfit_compute_helper_t*)puser;

  if (particleNumber > g_config.inconf[helper->config])
    error(1, "KIM get_neigh function requested a neighbor list for padding atoms which is not supported!\n");

  const int atom_idx = g_config.cnfstart[helper->config] + particleNumber;

  *numberOfNeighbors = g_config.atoms[atom_idx].num_neigh;
  *neighborsOfParticle = g_config.atoms[atom_idx].kim_neighbors;

  return 0;
}

int process_DEDr(const void* const dataObject, const double de, const double r,
                 const double* const dx, const int i, const int j)
{
  // don't do anything here
  // this function is just implemented for the case a model would require it
  return 0;
}

int process_D2EDr2(const void * const dataObject, const double de, const double* const r,
                    const double* const dx, const int* const i, const int* const j)
{
  // don't do anything here
  // this function is just implemented for the case a model would require it
  return 0;
}
