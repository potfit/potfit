/****************************************************************
 *
 * force_stiweb.c: Routines used for calculating Stillinger-Weber
 * 	forces/energies
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

#if !defined(STIWEB)
#error force_stiweb.c compiled without STIWEB support
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
 *  compute forces using Stillinger-Weber potentials with spline interpolation
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

double calc_forces(double* xi_opt, double* forces, int flag)
{
  const sw_t* sw = &g_pot.apot_table.sw;

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

    update_stiweb_pointers(xi_opt);

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

      // first loop over atoms: reset forces
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

      // 2nd loop: calculate pair forces and energies
      for (int atom_idx = 0; atom_idx < g_config.inconf[config_idx]; atom_idx++) {
        atom_t* atom = g_config.conf_atoms + atom_idx + g_config.cnfstart[config_idx] - g_mpi.firstatom;
        int n_i = 3 * (g_config.cnfstart[config_idx] + atom_idx);
        // loop over neighbors
        for (int neigh_idx = 0; neigh_idx < atom->num_neigh; neigh_idx++) {
          neigh_t* neigh_j = atom->neigh + neigh_idx;
          // In small cells, an atom might interact with itself
          int self = (neigh_j->nr == atom_idx + g_config.cnfstart[config_idx]) ? 1 : 0;

          // pair potential part
          int col = neigh_j->col[0];
          if (neigh_j->r < *(sw->a1[col])) {
            // fn value and grad are calculated in the same step
            double x[2] = { neigh_j->r, x[0] };
            double y[2] = { -sw->p[col][0], -sw->q[col][0] };
            double power[2] = { 0.0, 0.0 };
            power_m(2, power, x, y);
            double phi_r = sw->A[col][0] * power[0];
            double phi_a = -sw->B[col][0] * power[1];
            double inv_c = 1.0 / (neigh_j->r - sw->a1[col][0]);
            double f_cut = exp(sw->delta[col][0] * inv_c);
            double v2_val = (phi_r + phi_a) * f_cut;
            double v2_grad = 0.0;
            if (uf) {
              v2_grad = -v2_val * *(sw->delta[col]) * inv_c * inv_c -
                        f_cut * neigh_j->inv_r *
                            (*(sw->p[col]) * phi_r + *(sw->q[col]) * phi_a);
            }
            // avoid double counting if atom is interacting with itself
            if (self) {
              v2_val *= 0.5;
              v2_grad *= 0.5;
            }

            // only half cohesive energy because of full neighbor list
            forces[g_calc.energy_p + config_idx] += 0.5 * v2_val;

            if (uf) {
              vector tmp_force;
              tmp_force.x = neigh_j->dist_r.x * v2_grad;
              tmp_force.y = neigh_j->dist_r.y * v2_grad;
              tmp_force.z = neigh_j->dist_r.z * v2_grad;
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
          }

          // calculate for later
          if (neigh_j->r < sw->a2[col][0]) {
            double tmp_r = neigh_j->r - sw->a2[col][0];
            if (tmp_r < -0.01 * sw->gamma[col][0]) {
              tmp_r = 1.0 / tmp_r;
              neigh_j->f = exp(sw->gamma[col][0] * tmp_r);
              neigh_j->df = -neigh_j->f * sw->gamma[col][0] * tmp_r * tmp_r / neigh_j->r;
            } else {
              neigh_j->f = 0.0;
              neigh_j->df = 0.0;
            }
          }
        } // loop over neighbors j

        // loop over all neighbors
        for (int neigh_j_idx = 0; neigh_j_idx < atom->num_neigh - 1; neigh_j_idx++) {
          // Get pointer to neighbor j
          neigh_t* neigh_j = atom->neigh + neigh_j_idx;
          int ijk = neigh_j->ijk_start;
          // Force location for atom j
          int n_j = 3 * neigh_j->nr;
          // check if we are inside the cutoff radius
          if (neigh_j->r < sw->a2[neigh_j->col[0]][0]) {
            // loop over remaining neighbors
            for (int neigh_k_idx = neigh_j_idx + 1; neigh_k_idx < atom->num_neigh; neigh_k_idx++) {
              // Store pointer to angular part (g)
              angle_t* angle = atom->angle_part + ijk++;
              // Get pointer to second neighbor
              neigh_t* neigh_k = atom->neigh + neigh_k_idx;
              // store lambda for atom triple i,j,k
              double lambda = sw->lambda[atom->type][neigh_j->type][neigh_k->type][0];
              // shortcut for types without threebody interaction
              if (0.0 == lambda)
                continue;
              // Force location for atom k
              int n_k = 3 * neigh_k->nr;
              // check if we are inside the cutoff radius
              if (neigh_k->r < sw->a2[neigh_k->col[0]][0]) {
                // potential term
                double tmp = angle->cos + 1.0 / 3.0;
                double v3_val = lambda * neigh_j->f * neigh_k->f * tmp * tmp;

                // total potential
                forces[g_calc.energy_p + config_idx] += v3_val;

                // forces
                double tmp_grad1 = lambda * neigh_j->f * neigh_k->f * 2.0 * tmp;
                double tmp_grad2 = lambda * tmp * tmp;

                double tmp_jj = 1.0 / (neigh_j->r2);
                double tmp_jk = 1.0 / (neigh_j->r * neigh_k->r);
                double tmp_kk = 1.0 / (neigh_k->r2);
                double tmp_1 = tmp_grad2 * neigh_j->df * neigh_k->f - tmp_grad1 * angle->cos * tmp_jj;
                double tmp_2 = tmp_grad1 * tmp_jk;

                vector force_j;
                force_j.x = tmp_1 * neigh_j->dist.x + tmp_2 * neigh_k->dist.x;
                force_j.y = tmp_1 * neigh_j->dist.y + tmp_2 * neigh_k->dist.y;
                force_j.z = tmp_1 * neigh_j->dist.z + tmp_2 * neigh_k->dist.z;

                tmp_1 = tmp_grad2 * neigh_k->df * neigh_j->f - tmp_grad1 * angle->cos * tmp_kk;

                vector force_k;
                force_k.x = tmp_1 * neigh_k->dist.x + tmp_2 * neigh_j->dist.x;
                force_k.y = tmp_1 * neigh_k->dist.y + tmp_2 * neigh_j->dist.y;
                force_k.z = tmp_1 * neigh_k->dist.z + tmp_2 * neigh_j->dist.z;

                // update force on particle i
                forces[n_i + 0] += force_j.x + force_k.x;
                forces[n_i + 1] += force_j.y + force_k.y;
                forces[n_i + 2] += force_j.z + force_k.z;

                // update force on particle j
                forces[n_j + 0] -= force_j.x;
                forces[n_j + 1] -= force_j.y;
                forces[n_j + 2] -= force_j.z;

                // update force on particle k
                forces[n_k + 0] -= force_k.x;
                forces[n_k + 1] -= force_k.y;
                forces[n_k + 2] -= force_k.z;

#if defined(STRESS)
                // Distribute stress among atoms
                if (us) {
                  forces[stress_idx + 0] -= force_j.x * neigh_j->dist.x + force_k.x * neigh_k->dist.x;
                  forces[stress_idx + 1] -= force_j.y * neigh_j->dist.y + force_k.y * neigh_k->dist.y;
                  forces[stress_idx + 2] -= force_j.z * neigh_j->dist.z + force_k.z * neigh_k->dist.z;
                  forces[stress_idx + 3] -= 0.5 * (force_j.x * neigh_j->dist.y + force_k.x * neigh_k->dist.y + force_j.y * neigh_j->dist.x + force_k.y * neigh_k->dist.x);
                  forces[stress_idx + 4] -= 0.5 * (force_j.y * neigh_j->dist.z + force_k.y * neigh_k->dist.z + force_j.z * neigh_j->dist.y + force_k.z * neigh_k->dist.y);
                  forces[stress_idx + 5] -= 0.5 * (force_j.z * neigh_j->dist.x + force_k.z * neigh_k->dist.x + force_j.x * neigh_j->dist.z + force_k.x * neigh_k->dist.z);
                }
#endif  // STRESS
              }
            } // neigh_k_idx
          }
        } // neigh_j_idx
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

/* sum up forces */
#if defined(CONTRIB)
          if (atom->contrib)
#endif  // CONTRIB
            error_sum += g_config.conf_weight[config_idx] *
                      (dsquare(forces[n_i + 0]) + dsquare(forces[n_i + 1]) +
                        dsquare(forces[n_i + 2]));
        }
      }
      // end third loop over all atoms

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
          error_sum += g_config.conf_weight[config_idx] * g_param.sweight *
                    dsquare(forces[stress_idx + i]);
        }
      }
#endif  // STRESS
      // limiting constraints per configuration
    } // loop over configurations

    // dummy constraints (global)
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
        return 10e30;
      } else
        return error_sum;
    }
  } // infinite while loop

  // once a non-root process arrives here, all is done
  return -1.0;
}

/****************************************************************
 *
 *  update_stiweb_pointers
 *
 ****************************************************************/

void update_stiweb_pointers(double* xi)
{
  double* index = xi + 2;
  sw_t* sw = &g_pot.apot_table.sw;

  // allocate if this has not been done
  if (sw->init == 0) {
    sw->A = (double**)Malloc(g_calc.paircol * sizeof(double*));
    sw->B = (double**)Malloc(g_calc.paircol * sizeof(double*));
    sw->p = (double**)Malloc(g_calc.paircol * sizeof(double*));
    sw->q = (double**)Malloc(g_calc.paircol * sizeof(double*));
    sw->delta = (double**)Malloc(g_calc.paircol * sizeof(double*));
    sw->a1 = (double**)Malloc(g_calc.paircol * sizeof(double*));
    sw->gamma = (double**)Malloc(g_calc.paircol * sizeof(double*));
    sw->a2 = (double**)Malloc(g_calc.paircol * sizeof(double*));
    sw->lambda = (double****)Malloc(g_calc.paircol * sizeof(double***));
    for (int i = 0; i < g_param.ntypes; i++) {
      sw->lambda[i] = (double***)Malloc(g_param.ntypes * sizeof(double**));
      for (int j = 0; j < g_param.ntypes; j++) {
        sw->lambda[i][j] = (double**)Malloc(g_param.ntypes * sizeof(double*));
      }
    }
    sw->init = 1;
  }

  // update only if the address has changed
  if (sw->A[0] != index) {
    // set the pair parameters (stiweb_2)
    for (int i = 0; i < g_calc.paircol; i++) {
      sw->A[i] = index++;
      sw->B[i] = index++;
      sw->p[i] = index++;
      sw->q[i] = index++;
      sw->delta[i] = index++;
      sw->a1[i] = index++;
      index += 2;
    }
    // set the threebody parameters (stiweb_3)
    for (int i = 0; i < g_calc.paircol; i++) {
      sw->gamma[i] = index++;
      sw->a2[i] = index++;
      index += 2;
    }
    // set the lambda pointer
    for (int i = 0; i < g_param.ntypes; i++)
      for (int j = 0; j < g_param.ntypes; j++)
        for (int k = j; k < g_param.ntypes; k++) {
          sw->lambda[i][j][k] = sw->lambda[i][k][j] = index++;
        }
  }
}
