/****************************************************************
 *
 * force_tersoff.c: Routines used for calculating Tersoff forces/energies
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
 ****************************************************************/

#if !defined(TERSOFF)
#error force_tersoff.c compiled without TERSOFF support
#endif

#include "potfit.h"

#include "force.h"
#include "functions.h"
#include "memory.h"
#include "potential_input.h"
#include "potential_output.h"
#include "splines.h"
#include "utils.h"

/****************************************************************
 *
 *  init_forces
 *      called after all parameters and potentials are read
 *      additional assignments and initializations can be made here
 *
 ****************************************************************/

void init_force(int is_worker) {}

/****************************************************************
 *
 *  compute forces using pair potentials with spline interpolation
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
 *     g_pot.opt_pot.table - part of the struct g_pot.opt_pot, but it can also
 *be
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

#if !defined(TERSOFFMOD)

double calc_forces(double* xi_opt, double* forces, int flag)
{
  tersoff_t const* ters = &g_pot.apot_table.tersoff;

#if !defined(MPI)
  g_mpi.myconf = g_config.nconf;
#endif  // !MPI

  // This is the start of an infinite loop
  while (1) {
    // sum of squares of local process
    double error_sum = 0.0;

#if defined(MPI)
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (flag == 1)
      break; // Exception: flag 1 means clean up

    if (g_mpi.myid == 0)
      apot_check_params(xi_opt);
    MPI_Bcast(xi_opt, g_calc.ndimtot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    apot_check_params(xi_opt);
#endif  // MPI

    update_tersoff_pointers(xi_opt);

    // loop over configurations
    for (int config_idx = g_mpi.firstconf; config_idx < g_mpi.firstconf + g_mpi.myconf; config_idx++) {
      int uf = g_config.conf_uf[config_idx - g_mpi.firstconf];

      // reset energies and stresses
      forces[g_calc.energy_p + config_idx] = 0.0;

#if defined(STRESS)
      int us = g_config.conf_us[config_idx - g_mpi.firstconf];
      int stress_idx = g_calc.stress_p + 6 * config_idx;
      memset(forces + stress_idx, 0, 6 * sizeof(double));
#endif  // STRESS

      // first loop over all atoms: reset forces
      for (int atom_idx = 0; atom_idx < g_config.inconf[config_idx]; atom_idx++) {
        int n_i = 3 * (g_config.cnfstart[config_idx] + atom_idx);
        if (uf) {
          forces[n_i + 0] = -g_config.force_0[n_i + 0];
          forces[n_i + 1] = -g_config.force_0[n_i + 1];
          forces[n_i + 2] = -g_config.force_0[n_i + 2];
        } else {
          memset(forces + n_i, 0, 3 * sizeof(double));
        }
      }

      // second loop: calculate cutoff function f_c for all neighbor
      for (int atom_idx = 0; atom_idx < g_config.inconf[config_idx]; atom_idx++) {
        atom_t* atom = g_config.conf_atoms + atom_idx + g_config.cnfstart[config_idx] - g_mpi.firstatom;
        int n_i = 3 * (g_config.cnfstart[config_idx] + atom_idx);

        // loop over neighbors
        for (int neigh_idx = 0; neigh_idx < atom->num_neigh; neigh_idx++) {
          neigh_t* neigh_j = atom->neigh + neigh_idx;
          int col_j = neigh_j->col[0];
          // check if we are within the cutoff range
          if (neigh_j->r < ters->S[col_j][0]) {
            int self = (neigh_j->nr == atom_idx + g_config.cnfstart[config_idx]) ? 1 : 0;

            // calculate cutoff function f_c and store it for every neighbor
            double cut_tmp = M_PI / (ters->S[col_j][0] - ters->R[col_j][0]);
            double cut_tmp_j = cut_tmp * (neigh_j->r - ters->R[col_j][0]);
            if (neigh_j->r < ters->R[col_j][0]) {
              neigh_j->f = 1.0;
              neigh_j->df = 0.0;
            } else {
              neigh_j->f = 0.5 * (1.0 + cos(cut_tmp_j));
              neigh_j->df = -0.5 * cut_tmp * sin(cut_tmp_j);
            }

            // calculate pair part f_c*A*exp(-lambda*r) and the derivative
            double tmp = exp(-ters->lambda[col_j][0] * neigh_j->r);
            double phi_val = neigh_j->f * ters->A[col_j][0] * tmp;
            double phi_grad = neigh_j->df - ters->lambda[col_j][0] * neigh_j->f;
            phi_grad *= ters->A[col_j][0] * tmp;

            // avoid double counting if atom is interacting with itself
            if (self) {
              phi_val *= 0.5;
              phi_grad *= 0.5;
            }

            // only half cohesive energy because we have a full neighbor list
            forces[g_calc.energy_p + config_idx] += 0.5 * phi_val;

            if (uf) {
              // calculate pair forces
              vector tmp_force;
              tmp_force.x = neigh_j->dist_r.x * phi_grad;
              tmp_force.y = neigh_j->dist_r.y * phi_grad;
              tmp_force.z = neigh_j->dist_r.z * phi_grad;
              forces[n_i + 0] += tmp_force.x;
              forces[n_i + 1] += tmp_force.y;
              forces[n_i + 2] += tmp_force.z;
#if defined(STRESS)
              // also calculate pair stresses
              if (us) {
                forces[stress_idx + 0] -= 0.5 * neigh_j->dist.x * tmp_force.x;
                forces[stress_idx + 1] -= 0.5 * neigh_j->dist.y * tmp_force.y;
                forces[stress_idx + 2] -= 0.5 * neigh_j->dist.z * tmp_force.z;
                forces[stress_idx + 3] -= 0.5 * neigh_j->dist.x * tmp_force.y;
                forces[stress_idx + 4] -= 0.5 * neigh_j->dist.y * tmp_force.z;
                forces[stress_idx + 5] -= 0.5 * neigh_j->dist.z * tmp_force.x;
              }
#endif  // STRESS
            }
          } else {
            neigh_j->f = 0.0;
            neigh_j->df = 0.0;
          }
        } // loop over neighbors

        // calculate threebody part
        for (int neigh_j_idx = 0; neigh_j_idx < atom->num_neigh; neigh_j_idx++) {
          neigh_t* neigh_j = atom->neigh + neigh_j_idx;
          int col_j = neigh_j->col[0];
          // check if we are within the cutoff range
          if (neigh_j->r < ters->S[col_j][0]) {
            int ijk = neigh_j->ijk_start;
            int n_j = 3 * neigh_j->nr;

            // skip neighbor if coefficient is zero
            if (ters->B[col_j][0] == 0.0)
              continue;

            // reset variables for each neighbor
            double zeta = 0.0;
            vector dzeta_i = {0.0, 0.0, 0.0};
            vector dzeta_j = {0.0, 0.0, 0.0};

            // inner loop over neighbors
            for (int neigh_k_idx = 0; neigh_k_idx < atom->num_neigh; neigh_k_idx++) {
              if (neigh_k_idx == neigh_j_idx)
                continue;
              neigh_t* neigh_k = atom->neigh + neigh_k_idx;
              int col_k = neigh_k->col[0];
              angle_t* angle = atom->angle_part + ijk++;
              if (neigh_k->r < ters->S[col_k][0]) {
                double tmp_jk = 1.0 / (neigh_j->r * neigh_k->r);

                double tmp_1 = ters->h[col_j][0] - angle->cos;
                double tmp_2 = 1.0 / (ters->d2[col_j] + tmp_1 * tmp_1);
                double g_theta = 1.0 + ters->c2[col_j] / ters->d2[col_j] -
                          ters->c2[col_j] * tmp_2;

                // zeta
                zeta += neigh_k->f * ters->omega[col_k][0] * g_theta;

                double tmp_j2 = angle->cos / (neigh_j->r * neigh_j->r);
                double tmp_k2 = angle->cos / (neigh_k->r * neigh_k->r);

                vector dcos_j;
                dcos_j.x = tmp_jk * neigh_k->dist.x - tmp_j2 * neigh_j->dist.x;
                dcos_j.y = tmp_jk * neigh_k->dist.y - tmp_j2 * neigh_j->dist.y;
                dcos_j.z = tmp_jk * neigh_k->dist.z - tmp_j2 * neigh_j->dist.z;

                vector dcos_k;
                dcos_k.x = tmp_jk * neigh_j->dist.x - tmp_k2 * neigh_k->dist.x;
                dcos_k.y = tmp_jk * neigh_j->dist.y - tmp_k2 * neigh_k->dist.y;
                dcos_k.z = tmp_jk * neigh_j->dist.z - tmp_k2 * neigh_k->dist.z;

                double tmp_3 = 2.0 * ters->c2[col_j] * tmp_1 * tmp_2 * tmp_2 * neigh_k->f * ters->omega[col_k][0];

                double tmp_grad = neigh_k->df / neigh_k->r * g_theta * ters->omega[col_k][0];

                neigh_k->dzeta.x = tmp_grad * neigh_k->dist.x - tmp_3 * dcos_k.x;
                neigh_k->dzeta.y = tmp_grad * neigh_k->dist.y - tmp_3 * dcos_k.y;
                neigh_k->dzeta.z = tmp_grad * neigh_k->dist.z - tmp_3 * dcos_k.z;

                dzeta_i.x -= neigh_k->dzeta.x;
                dzeta_i.y -= neigh_k->dzeta.y;
                dzeta_i.z -= neigh_k->dzeta.z;

                dzeta_j.x -= tmp_3 * dcos_j.x;
                dzeta_j.y -= tmp_3 * dcos_j.y;
                dzeta_j.z -= tmp_3 * dcos_j.z;
              }
            } // neigh_k_idx

            double phi_a = 0.5 * ters->B[col_j][0] * exp(-ters->mu[col_j][0] * neigh_j->r);

            double tmp_pow_1 = ters->gamma[col_j][0] * zeta;
            double tmp_4 = 0.0;
            power_1(&tmp_4, &tmp_pow_1, ters->n[col_j]);

            tmp_pow_1 = 1.0 + tmp_4;
            double tmp_pow_2 = -1.0 / (2.0 * ters->n[col_j][0]);
            double b_ij = 0.0;
            power_1(&b_ij, &tmp_pow_1, &tmp_pow_2);

            double phi_val = -b_ij * phi_a;

            forces[g_calc.energy_p + config_idx] += neigh_j->f * phi_val;

            double tmp_5 = 0.0;
            if (zeta != 0.0)
              tmp_5 = -b_ij * neigh_j->f * phi_a * tmp_4 / (2.0 * zeta * (1.0 + tmp_4));
            double tmp_6 = (neigh_j->f * phi_a * ters->mu[col_j][0] * b_ij + neigh_j->df * phi_val) / neigh_j->r;

            vector force_j;
            force_j.x = -tmp_6 * neigh_j->dist.x + tmp_5 * dzeta_j.x;
            force_j.y = -tmp_6 * neigh_j->dist.y + tmp_5 * dzeta_j.y;
            force_j.z = -tmp_6 * neigh_j->dist.z + tmp_5 * dzeta_j.z;

            for (int neigh_k_idx = 0; neigh_k_idx < atom->num_neigh; neigh_k_idx++) {
              if (neigh_k_idx != neigh_j_idx) {
                neigh_t* neigh_k = atom->neigh + neigh_k_idx;
                int col_k = neigh_k->col[0];
                if (neigh_k->r < ters->S[col_k][0]) {
                  int n_k = 3 * neigh_k->nr;
                  // update force on particle k
                  forces[n_k + 0] += tmp_5 * neigh_k->dzeta.x;
                  forces[n_k + 1] += tmp_5 * neigh_k->dzeta.y;
                  forces[n_k + 2] += tmp_5 * neigh_k->dzeta.z;

#if defined(STRESS)
                  if (us) {
                    // Distribute stress among atoms
                    forces[stress_idx + 0] += neigh_k->dist.x * tmp_5 * neigh_k->dzeta.x;
                    forces[stress_idx + 1] += neigh_k->dist.y * tmp_5 * neigh_k->dzeta.y;
                    forces[stress_idx + 2] += neigh_k->dist.z * tmp_5 * neigh_k->dzeta.z;
                    forces[stress_idx + 3] += 0.5 * tmp_5 * (neigh_k->dist.x * neigh_k->dzeta.y + neigh_k->dist.y * neigh_k->dzeta.x);
                    forces[stress_idx + 4] += 0.5 * tmp_5 * (neigh_k->dist.y * neigh_k->dzeta.z + neigh_k->dist.z * neigh_k->dzeta.y);
                    forces[stress_idx + 5] += 0.5 * tmp_5 * (neigh_k->dist.z * neigh_k->dzeta.x + neigh_k->dist.x * neigh_k->dzeta.z);
                  }
#endif  // STRESS
                }
              } // neigh_k_idx != neigh_j_idx
            }   // neigh_k_idx loop

            // update force on particle j
            forces[n_j + 0] += force_j.x;
            forces[n_j + 1] += force_j.y;
            forces[n_j + 2] += force_j.z;

            // update force on particle i
            forces[n_i + 0] += tmp_5 * dzeta_i.x - force_j.x;
            forces[n_i + 1] += tmp_5 * dzeta_i.y - force_j.y;
            forces[n_i + 2] += tmp_5 * dzeta_i.z - force_j.z;

#if defined(STRESS)
            if (us) {
              // Distribute stress among atoms
              forces[stress_idx + 0] += neigh_j->dist.x * force_j.x;
              forces[stress_idx + 1] += neigh_j->dist.y * force_j.y;
              forces[stress_idx + 2] += neigh_j->dist.z * force_j.z;
              forces[stress_idx + 3] += 0.5 * (neigh_j->dist.x * force_j.y + neigh_j->dist.y * force_j.x);
              forces[stress_idx + 4] += 0.5 * (neigh_j->dist.y * force_j.z + neigh_j->dist.z * force_j.y);
              forces[stress_idx + 5] += 0.5 * (neigh_j->dist.z * force_j.x + neigh_j->dist.x * force_j.z);
            }
#endif  // STRESS
          } // neigh_j_idx
        }
      } // end second loop over all atoms

      // third loop over all atoms, sum up forces
      if (uf) {
        for (int atom_idx = 0; atom_idx < g_config.inconf[config_idx]; atom_idx++) {
          atom_t* atom = g_config.conf_atoms + atom_idx + g_config.cnfstart[config_idx] - g_mpi.firstatom;
          int n_i = 3 * (g_config.cnfstart[config_idx] + atom_idx);
#if defined(FWEIGHT)
          // Weigh by absolute value of force
          forces[n_i + 0] /= FORCE_EPS + atom->absforce;
          forces[n_i + 1] /= FORCE_EPS + atom->absforce;
          forces[n_i + 2] /= FORCE_EPS + atom->absforce;
#endif  // FWEIGHT

// sum up forces
#if defined(CONTRIB)
          if (atom->contrib)
#endif  // CONTRIB
            error_sum += g_config.conf_weight[config_idx] * (dsquare(forces[n_i + 0]) + dsquare(forces[n_i + 1]) + dsquare(forces[n_i + 2]));
        }
      } // end third loop over all atoms

      // energy contributions
      forces[g_calc.energy_p + config_idx] /= (double)g_config.inconf[config_idx];
      forces[g_calc.energy_p + config_idx] -= g_config.force_0[g_calc.energy_p + config_idx];
      error_sum += g_config.conf_weight[config_idx] * g_param.eweight * dsquare(forces[g_calc.energy_p + config_idx]);

#if defined(STRESS)
      // stress contributions
      if (uf && us) {
        for (int i = 0; i < 6; i++) {
          forces[stress_idx + i] /= g_config.conf_vol[config_idx - g_mpi.firstconf];
          forces[stress_idx + i] -= g_config.force_0[stress_idx + i];
          error_sum += g_config.conf_weight[config_idx] * g_param.sweight * dsquare(forces[stress_idx + i]);
        }
      }
#endif  // STRESS
    } // loop over configurations

    // add punishment for out of bounds (mostly for powell_lsq)
    if (g_mpi.myid == 0)
      error_sum += apot_punish(xi_opt, forces);

    gather_forces(&error_sum, forces);

    // root process exits this function now
    if (g_mpi.myid == 0) {
      // Increase function call counter
      g_calc.fcalls++;
      if (isnan(error_sum)) {
#if defined(DEBUG)
        printf("\n--> Force is nan! <--\n\n");
#endif  // DEBUG
        return 10e10;
      } else
        return error_sum;
    }

  } // end of infinite loop

  // once a non-root process arrives here, all is done
  return -1.0;
}

#else  // !TERSOFFMOD

/****************************************************************
  calc_forces
****************************************************************/

double calc_forces(double* xi_opt, double* forces, int flag)
{
  const tersoff_t* ters = &g_pot.apot_table.tersoff;

#if !defined(MPI)
  g_mpi.myconf = g_config.nconf;
#endif  // !MPI

  // This is the start of an infinite loop
  while (1) {
    // sum of squares of local process
    double error_sum = 0.0;

#if defined(MPI)
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (flag == 1)
      break; // Exception: flag 1 means clean up

    if (g_mpi.myid == 0)
      apot_check_params(xi_opt);
    MPI_Bcast(xi_opt, g_calc.ndimtot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    apot_check_params(xi_opt);
#endif  // MPI

    update_tersoff_pointers(xi_opt);

    // loop over configurations
    for (int config_idx = g_mpi.firstconf; config_idx < g_mpi.firstconf + g_mpi.myconf; config_idx++) {
      int uf = g_config.conf_uf[config_idx - g_mpi.firstconf];

      // reset energies and stresses
      forces[g_calc.energy_p + config_idx] = 0.0;

#if defined(STRESS)
      int us = g_config.conf_us[config_idx - g_mpi.firstconf];
      int stress_idx = g_calc.stress_p + 6 * config_idx;
      memset(forces + stress_idx, 0, 6 * sizeof(double));
#endif  // STRESS

      // first loop over all atoms: reset forces
      for (int atom_idx = 0; atom_idx < g_config.inconf[config_idx]; atom_idx++) {
        int n_i = 3 * (g_config.cnfstart[config_idx] + atom_idx);
        if (uf) {
          forces[n_i + 0] = -g_config.force_0[n_i + 0];
          forces[n_i + 1] = -g_config.force_0[n_i + 1];
          forces[n_i + 2] = -g_config.force_0[n_i + 2];
        } else {
          memset(forces + n_i, 0, 3 * sizeof(double));
        }
      }

      // second loop: calculate cutoff function f_c for all neighbors
      for (int atom_idx = 0; atom_idx < g_config.inconf[config_idx]; atom_idx++) {
        atom_t* atom = g_config.conf_atoms + atom_idx + g_config.cnfstart[config_idx] - g_mpi.firstatom;
        int n_i = 3 * (g_config.cnfstart[config_idx] + atom_idx);

        // loop over neighbors
        // calculate pair potential part: f*A*exp(-lambda*r)
        for (int neigh_idx = 0; neigh_idx < atom->num_neigh; neigh_idx++) {
          neigh_t* neigh_j = atom->neigh + neigh_idx;
          int col_j = neigh_j->col[0];
          // check if we are within the cutoff range
          if (neigh_j->r < ters->R2[col_j][0]) {
            int self = (neigh_j->nr == atom_idx + g_config.cnfstart[config_idx]) ? 1 : 0;

            // calculate cutoff function f_c and store it for every neighbor
            double tmp_1 = M_PI / (ters->R2[col_j][0] - ters->R1[col_j][0]);
            double tmp_2 = tmp_1 * (neigh_j->r - ters->R1[col_j][0]);
            if (neigh_j->r < ters->R1[col_j][0]) {
              neigh_j->f = 1.0;
              neigh_j->df = 0.0;
            } else {
              neigh_j->f = 0.5 * (1.0 + 1.125 * cos(tmp_2) - 0.125 * cos(3.0 * tmp_2));
              neigh_j->df = -0.5 * tmp_1 * (1.125 * sin(tmp_2) - 0.375 * sin(3.0 * tmp_2));
            }

            // calculate pair part f_c*A*exp(-lambda*r) and the derivative
            tmp_1 = exp(-ters->lambda[col_j][0] * neigh_j->r);
            double phi_val = neigh_j->f * ters->A[col_j][0] * tmp_1;
            double phi_grad = neigh_j->df - ters->lambda[col_j][0] * neigh_j->f;
            phi_grad *= ters->A[col_j][0] * tmp_1;

            // avoid double counting if atom is interacting with itself
            if (self) {
              phi_val *= 0.5;
              phi_grad *= 0.5;
            }

            // only half cohesive energy because we have a full neighbor list
            forces[g_calc.energy_p + config_idx] += 0.5 * phi_val;

            if (uf) {
              // calculate pair forces
              vector tmp_force;
              tmp_force.x = neigh_j->dist_r.x * phi_grad;
              tmp_force.y = neigh_j->dist_r.y * phi_grad;
              tmp_force.z = neigh_j->dist_r.z * phi_grad;
              forces[n_i + 0] += tmp_force.x;
              forces[n_i + 1] += tmp_force.y;
              forces[n_i + 2] += tmp_force.z;
#if defined(STRESS)
              // also calculate pair stresses
              if (us) {
                forces[stress_idx + 0] -= 0.5 * neigh_j->dist.x * tmp_force.x;
                forces[stress_idx + 1] -= 0.5 * neigh_j->dist.y * tmp_force.y;
                forces[stress_idx + 2] -= 0.5 * neigh_j->dist.z * tmp_force.z;
                forces[stress_idx + 3] -= 0.5 * neigh_j->dist.x * tmp_force.y;
                forces[stress_idx + 4] -= 0.5 * neigh_j->dist.y * tmp_force.z;
                forces[stress_idx + 5] -= 0.5 * neigh_j->dist.z * tmp_force.x;
              }
#endif  // STRESS
            }
          } else {
            neigh_j->f = 0.0;
            neigh_j->df = 0.0;
          }
        } // loop over neighbors

        // loop over neighbors
        // calculate threebody part
        for (int neigh_j_idx = 0; neigh_j_idx < atom->num_neigh; neigh_j_idx++) {
          neigh_t* neigh_j = atom->neigh + neigh_j_idx;
          int col_j = neigh_j->col[0];
          // check if we are within the cutoff range
          if (neigh_j->r < ters->R2[col_j][0]) {
            int ijk = neigh_j->ijk_start;
            int n_j = 3 * neigh_j->nr;

            // skip neighbor if coefficient is zero
            if (ters->B[col_j][0] == 0.0)
              continue;

            // reset variables for each neighbor
            double zeta = 0.0;
            double dzeta_ij = 0.0;
            vector dzeta_i;
            dzeta_i.x = 0.0;
            dzeta_i.y = 0.0;
            dzeta_i.z = 0.0;
            vector dzeta_j;
            dzeta_j.x = 0.0;
            dzeta_j.y = 0.0;
            dzeta_j.z = 0.0;

            // inner loop over neighbors
            for (int neigh_k_idx = 0; neigh_k_idx < atom->num_neigh; neigh_k_idx++) {
              if (neigh_j_idx == neigh_k_idx)
                continue;
              neigh_t* neigh_k = atom->neigh + neigh_k_idx;
              int col_k = neigh_k->col[0];
              angle_t* angle = atom->angle_part + ijk++;

              if (neigh_k->r < ters->R2[col_k][0]) {
                vector dcos_j;
                dcos_j.x = (neigh_k->dist_r.x - neigh_j->dist_r.x * angle->cos) / neigh_j->r;
                dcos_j.y = (neigh_k->dist_r.y - neigh_j->dist_r.y * angle->cos) / neigh_j->r;
                dcos_j.z = (neigh_k->dist_r.z - neigh_j->dist_r.z * angle->cos) / neigh_j->r;
                vector dcos_k;
                dcos_k.x = (neigh_j->dist_r.x - neigh_k->dist_r.x * angle->cos) / neigh_k->r;
                dcos_k.y = (neigh_j->dist_r.y - neigh_k->dist_r.y * angle->cos) / neigh_k->r;
                dcos_k.z = (neigh_j->dist_r.z - neigh_k->dist_r.z * angle->cos) / neigh_k->r;

                // g(theta)
                double tmp_1 = ters->h[col_j][0] - angle->cos;
                double tmp_2 = 1.0 / (ters->c3[col_j][0] + tmp_1 * tmp_1);
                double tmp_3 = ters->c4[col_j][0] * exp(-ters->c5[col_j][0] * tmp_1 * tmp_1);

                double g_theta = ters->c1[col_j][0] + ters->c2[col_j][0] * tmp_1 * tmp_1 * tmp_2 * (1.0 + tmp_3);
                double dg_theta = 2.0 * ters->c2[col_j][0] * tmp_1 * tmp_2 * (ters->c5[col_j][0] * tmp_1 * tmp_1 * tmp_3 - ters->c3[col_j][0] * tmp_2 * (1.0 + tmp_3));

                tmp_1 = neigh_j->r - neigh_k->r;
                tmp_2 = ters->alpha[col_j][0] * ters->beta[col_j][0] * pow(tmp_1, ters->beta[col_j][0] - 1.0);
                tmp_3 = exp(ters->alpha[col_j][0] * pow(tmp_1, ters->beta[col_j][0]));

                double dzeta_ik = (neigh_k->df - neigh_k->f * tmp_2) * g_theta * tmp_3;

                // zeta
                double tmp_4 = neigh_k->f * g_theta * tmp_3;

                zeta += tmp_4;
                dzeta_ij += tmp_4 * tmp_2;
                double dzeta_cos = neigh_k->f * dg_theta * tmp_3;

                neigh_k->dzeta.x = dzeta_cos * dcos_k.x + dzeta_ik * neigh_k->dist_r.x;
                neigh_k->dzeta.y = dzeta_cos * dcos_k.y + dzeta_ik * neigh_k->dist_r.y;
                neigh_k->dzeta.z = dzeta_cos * dcos_k.z + dzeta_ik * neigh_k->dist_r.z;

                dzeta_i.x -= neigh_k->dzeta.x;
                dzeta_i.y -= neigh_k->dzeta.y;
                dzeta_i.z -= neigh_k->dzeta.z;

                dzeta_j.x += dzeta_cos * dcos_j.x;
                dzeta_j.y += dzeta_cos * dcos_j.y;
                dzeta_j.z += dzeta_cos * dcos_j.z;
              }
            } // neigh_k_idx

            dzeta_j.x += dzeta_ij * neigh_j->dist_r.x;
            dzeta_j.y += dzeta_ij * neigh_j->dist_r.y;
            dzeta_j.z += dzeta_ij * neigh_j->dist_r.z;

            double tmp_1 = pow(zeta, ters->eta[col_j][0]);
            double b = pow(1.0 + tmp_1, -ters->delta[col_j][0]);

            double tmp_2 = 0.5 * b * ters->B[col_j][0] * exp(-ters->mu[col_j][0] * neigh_j->r);

            double tmp_3 = 0.0;
            if (zeta != 0.0)
              tmp_3 = tmp_2 * neigh_j->f * ters->eta[col_j][0] * ters->delta[col_j][0] * tmp_1 / ((1.0 + tmp_1) * zeta);

            double phi_val = -tmp_2;
            double tmp_4 = -tmp_2 * (neigh_j->df - *(ters->mu[col_j]) * neigh_j->f);

            forces[g_calc.energy_p + config_idx] += neigh_j->f * phi_val;

            vector force_j;
            force_j.x = -tmp_4 * neigh_j->dist_r.x - tmp_3 * dzeta_j.x;
            force_j.y = -tmp_4 * neigh_j->dist_r.y - tmp_3 * dzeta_j.y;
            force_j.z = -tmp_4 * neigh_j->dist_r.z - tmp_3 * dzeta_j.z;

            // update force on particle j
            forces[n_j + 0] += force_j.x;
            forces[n_j + 1] += force_j.y;
            forces[n_j + 2] += force_j.z;

            // update force on particle i
            forces[n_i + 0] -= tmp_3 * dzeta_i.x + force_j.x;
            forces[n_i + 1] -= tmp_3 * dzeta_i.y + force_j.y;
            forces[n_i + 2] -= tmp_3 * dzeta_i.z + force_j.z;

#if defined(STRESS)
            if (us) {
              // Distribute stress among atoms
              double tmp = neigh_j->dist.x * force_j.x;
              forces[stress_idx + 0] += tmp;
              tmp = neigh_j->dist.y * force_j.y;
              forces[stress_idx + 1] += tmp;
              tmp = neigh_j->dist.z * force_j.z;
              forces[stress_idx + 2] += tmp;
              tmp = 0.5 * (neigh_j->dist.x * force_j.y + neigh_j->dist.y * force_j.x);
              forces[stress_idx + 3] += tmp;
              tmp = 0.5 * (neigh_j->dist.y * force_j.z + neigh_j->dist.z * force_j.y);
              forces[stress_idx + 4] += tmp;
              tmp = 0.5 * (neigh_j->dist.z * force_j.x + neigh_j->dist.x * force_j.z);
              forces[stress_idx + 5] += tmp;
            }
#endif  // STRESS

            for (int neigh_k_idx = 0; neigh_k_idx < atom->num_neigh; neigh_k_idx++) {
              if (neigh_k_idx != neigh_j_idx) {
                neigh_t* neigh_k = atom->neigh + neigh_k_idx;
                int col_k = neigh_k->col[0];
                if (neigh_k->r < ters->R2[col_k][0]) {
                  int n_k = 3 * neigh_k->nr;
                  // update force on particle k
                  forces[n_k + 0] -= tmp_3 * neigh_k->dzeta.x;
                  forces[n_k + 1] -= tmp_3 * neigh_k->dzeta.y;
                  forces[n_k + 2] -= tmp_3 * neigh_k->dzeta.z;

#if defined(STRESS)
                  if (us) {
                    // Distribute stress among atoms
                    double tmp = neigh_k->dist.x * tmp_3 * neigh_k->dzeta.x;
                    forces[stress_idx + 0] -= tmp;
                    tmp = neigh_k->dist.y * tmp_3 * neigh_k->dzeta.y;
                    forces[stress_idx + 1] -= tmp;
                    tmp = neigh_k->dist.z * tmp_3 * neigh_k->dzeta.z;
                    forces[stress_idx + 2] -= tmp;
                    tmp = 0.5 * tmp_3 * (neigh_k->dist.x * neigh_k->dzeta.y + neigh_k->dist.y * neigh_k->dzeta.x);
                    forces[stress_idx + 3] -= tmp;
                    tmp = 0.5 * tmp_3 * (neigh_k->dist.y * neigh_k->dzeta.z + neigh_k->dist.z * neigh_k->dzeta.y);
                    forces[stress_idx + 4] -= tmp;
                    tmp = 0.5 * tmp_3 * (neigh_k->dist.z * neigh_k->dzeta.x + neigh_k->dist.x * neigh_k->dzeta.z);
                    forces[stress_idx + 5] -= tmp;
                  }
#endif  // STRESS
                }
              } // neigh_k_idx != neigh_j_idx
            }   // neigh_k_idx loop
          }     // neigh_f_idx
        }
      } // end second loop over all atoms

      // third loop over all atoms, sum up forces
      if (uf) {
        for (int atom_idx = 0; atom_idx < g_config.inconf[config_idx]; atom_idx++) {
          atom_t* atom = g_config.conf_atoms + atom_idx + g_config.cnfstart[config_idx] - g_mpi.firstatom;
          int n_i = 3 * (g_config.cnfstart[config_idx] + atom_idx);
#if defined(FWEIGHT)
          // Weigh by absolute value of force
          forces[n_i + 0] /= FORCE_EPS + atom->absforce;
          forces[n_i + 1] /= FORCE_EPS + atom->absforce;
          forces[n_i + 2] /= FORCE_EPS + atom->absforce;
#endif  // FWEIGHT

// sum up forces
#if defined(CONTRIB)
          if (atom->contrib)
#endif  // CONTRIB
            error_sum += g_config.conf_weight[config_idx] * (dsquare(forces[n_i + 0]) + dsquare(forces[n_i + 1]) + dsquare(forces[n_i + 2]));
        }
      } // end third loop over all atoms

      // energy contributions
      forces[g_calc.energy_p + config_idx] /= (double)g_config.inconf[config_idx];
      forces[g_calc.energy_p + config_idx] -= g_config.force_0[g_calc.energy_p + config_idx];
      error_sum += g_config.conf_weight[config_idx] * g_param.eweight *
                dsquare(forces[g_calc.energy_p + config_idx]);

#if defined(STRESS)
      // stress contributions
      if (uf && us) {
        for (int i = 0; i < 6; i++) {
          forces[stress_idx + i] /= g_config.conf_vol[config_idx - g_mpi.firstconf];
          forces[stress_idx + i] -= g_config.force_0[stress_idx + i];
          error_sum += g_config.conf_weight[config_idx] * g_param.sweight * dsquare(forces[stress_idx + i]);
        }
      }
#endif  // STRESS
    } // loop over configurations

    // add punishment for out of bounds (mostly for powell_lsq)
    if (g_mpi.myid == 0)
      error_sum += apot_punish(xi_opt, forces);

    gather_forces(&error_sum, forces);

    // root process exits this function now
    if (g_mpi.myid == 0) {
      // Increase function call counter
      g_calc.fcalls++;
      if (isnan(error_sum)) {
#if defined(DEBUG)
        printf("\n--> Force is nan! <--\n\n");
#endif  // DEBUG
        return 10e10;
      } else
        return error_sum;
    }
  } // end of infinite loop

  // once a non-root process arrives here, all is done
  return -1.0;
}

#endif  // !TERSOFFMOD

/****************************************************************
  update_tersoff_pointers
****************************************************************/

void update_tersoff_pointers(double* xi)
{
  double* index = xi + 2;
  tersoff_t* tersoff = &g_pot.apot_table.tersoff;

  // allocate if this has not been done
  if (tersoff->init == 0) {
    tersoff->A = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->B = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->lambda = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->mu = (double**)Malloc(g_calc.paircol * sizeof(double*));
#if defined(TERSOFFMOD)
    tersoff->eta = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->delta = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->alpha = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->beta = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->c1 = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->c2 = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->c3 = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->c4 = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->c5 = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->h = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->R1 = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->R2 = (double**)Malloc(g_calc.paircol * sizeof(double*));
#else
    tersoff->gamma = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->n = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->c = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->d = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->h = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->S = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->R = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->chi = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->omega = (double**)Malloc(g_calc.paircol * sizeof(double*));
    tersoff->c2 = (double*)Malloc(g_calc.paircol * sizeof(double));
    tersoff->d2 = (double*)Malloc(g_calc.paircol * sizeof(double));
    tersoff->one = 1.0;
#endif // TERSOFFMOD
    tersoff->init = 1;
  }

  // update only if the address has changed
  if (tersoff->A[0] != index) {
    // set the pair parameters
    for (int i = 0; i < g_calc.paircol; i++) {
      tersoff->A[i] = index++;
      tersoff->B[i] = index++;
      tersoff->lambda[i] = index++;
      tersoff->mu[i] = index++;
#if defined(TERSOFFMOD)
      tersoff->eta[i] = index++;
      tersoff->delta[i] = index++;
      tersoff->alpha[i] = index++;
      tersoff->beta[i] = index++;
      tersoff->c1[i] = index++;
      tersoff->c2[i] = index++;
      tersoff->c3[i] = index++;
      tersoff->c4[i] = index++;
      tersoff->c5[i] = index++;
      tersoff->h[i] = index++;
      tersoff->R1[i] = index++;
      tersoff->R2[i] = index++;
#else
      tersoff->gamma[i] = index++;
      tersoff->n[i] = index++;
      tersoff->c[i] = index++;
      tersoff->d[i] = index++;
      tersoff->h[i] = index++;
      tersoff->S[i] = index++;
      tersoff->R[i] = index++;
#endif
      index += 2;
    }
#if !defined(TERSOFFMOD)
    for (int i = 0; i < g_calc.paircol; i++) {
      if (0 == (i % g_param.ntypes)) {
        tersoff->chi[i] = &tersoff->one;
        tersoff->omega[i] = &tersoff->one;
      } else {
        tersoff->chi[i] = index++;
        tersoff->omega[i] = index++;
        index += 2;
      }
    }
  }

  // calculate c2 and d2
  for (int i = 0; i < g_calc.paircol; i++) {
    tersoff->c2[i] = dsquare(tersoff->c[i][0]);
    tersoff->d2[i] = dsquare(tersoff->d[i][0]);
#endif
  }
}
