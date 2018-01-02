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
    called after all parameters and potentials are read
    additional assignments and initializations can be performed here
****************************************************************/

void init_force(int is_worker)
{
  // nothing to do here for KIM potentials
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
  int status;

  // publish KIM parameters
  for (int i = 0; i < g_config.nconf; ++i) {
    int idx = 0;
    for (int j = 0; j < g_pot.apot_table.total_par; ++j) {
      double* data = ((double**)g_kim.param_value[i])[g_pot.apot_table.idxparam[j]];
      for (int k = 0; k < g_kim.freeparams.size[g_kim.idx_opt_param[j]]; ++k) {
        data[k] = xi_opt[g_pot.opt_pot.idx[idx++]];
      }
    }
    status = KIM_API_model_reinit(g_kim.pkim[i]);
    if (KIM_STATUS_OK > status)
      error(1, "KIM_API_model_reinit failed!\n");
  }

  // This is the start of an infinite loop
  while (1) {
    double error_sum = 0.0; // sum of squares of local process

#if !defined(MPI)
    g_mpi.myconf = g_config.nconf;
#endif // MPI

    // loop over configurations
    for (int h = g_mpi.firstconf; h < g_mpi.firstconf + g_mpi.myconf; h++) {
      double* kim_energy = NULL;
      double* kim_forces = NULL;
      int uf = g_config.conf_uf[h - g_mpi.firstconf];
#if defined(STRESS)
      double* virial = NULL;
      int us = g_config.conf_us[h - g_mpi.firstconf];
#endif // STRESS

      // Calculate forces from KIM (general forces, including forces, virial and energy)

      // get KIM data
#if defined(STRESS)
      KIM_API_getm_data(g_kim.pkim[h], &status, 3*3,
                        "energy", &kim_energy, 1,
                        "forces", &kim_forces,  uf,
                        "virial", &virial, us);
      if (KIM_STATUS_OK > status)
        error(1, "KIM_API_getm_data failed!\n");
#else
      KIM_API_getm_data(g_kim.pkim[h], &status, 2*3,
                        "energy", &kim_energy, 1,
                        "forces", &kim_forces,  uf);
      if (KIM_STATUS_OK > status)
        error(1, "KIM_API_getm_data failed!\n");
#endif

      // Call model compute
      status = KIM_API_model_compute(g_kim.pkim[h]);
      if (KIM_STATUS_OK > status)
        error(1, "KIM_API_model_compute failed!\n");

      // forces contributation
      double weight = sqrt(g_config.conf_weight[h]);
      for (int i = 0; i < g_config.inconf[h]; i++) {
#if defined(FWEIGHT) || defined(CONTRIB)
        atom_t* atom = g_config.conf_atoms + i + g_config.cnfstart[h] - g_mpi.firstatom;
#endif // FWEIGHT || CONTRIB
        int n_i = DIM * (g_config.cnfstart[h] + i);
        if (uf) {
          forces[n_i + 0] = weight * (kim_forces[DIM * i + 0] - g_config.force_0[n_i + 0]);
          forces[n_i + 1] = weight * (kim_forces[DIM * i + 1] - g_config.force_0[n_i + 1]);
          forces[n_i + 2] = weight * (kim_forces[DIM * i + 2] - g_config.force_0[n_i + 2]);
        } else {
          forces[n_i + 0] = 0.0;
          forces[n_i + 1] = 0.0;
          forces[n_i + 2] = 0.0;
        }

#if defined(FWEIGHT)
        // weigh by absolute value of force
        forces[n_i + 0] /= FORCE_EPS + atom->absforce;
        forces[n_i + 1] /= FORCE_EPS + atom->absforce;
        forces[n_i + 2] /= FORCE_EPS + atom->absforce;
#endif // FWEIGHT

#if defined(CONTRIB)
        if (atom->contrib)
#endif // CONTRIB
          error_sum += (dsquare(forces[n_i + 0])
                  + dsquare(forces[n_i + 1])
                  + dsquare(forces[n_i + 2]) );
      }

      // reset energies and stresses
      forces[g_calc.energy_p + h] = 0.0;
#if defined(STRESS)
      int stresses_offset = g_calc.stress_p + 6 * h;
      for (int i = 0; i < 6; i++)
        forces[stresses_offset + i] = 0.0;
#endif /* STRESS */

      weight = sqrt(g_config.conf_weight[h] * g_param.eweight);
      forces[g_calc.energy_p + h] = weight * (*kim_energy) / (double)g_config.inconf[h];
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
