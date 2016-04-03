/****************************************************************
 *
 * force_adp.c: Routine used for calculating adp forces/energies
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

#if !defined(ADP)
#error force_adp.c compiled without ADP support
#endif

#include "potfit.h"

#include "force.h"
#include "functions.h"
#include "memory.h"
#if defined(MPI)
#include "mpi_utils.h"
#endif
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

void init_force(int is_worker)
{
#if !defined(APOT)
  // set spline density corrections to 0
  g_calc.lambda = (double*)Malloc(g_param.ntypes * sizeof(double));
#endif  // APOT
}

/****************************************************************
 *
 *  compute forces using adp potentials with spline interpolation
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
 * flag is an integer controlling the behaviour of calc_forces_adp.
 *    flag == 1 will cause all processes to exit calc_forces_adp after
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
  double* xi = NULL;

  switch (g_pot.format_type) {
    case POTENTIAL_FORMAT_UNKNOWN:
      error(1, "Unknown potential format detected! (%s:%d)\n", __FILE__, __LINE__);
    case POTENTIAL_FORMAT_ANALYTIC:
      xi = g_pot.calc_pot.table;
      break;
    case POTENTIAL_FORMAT_TABULATED_EQ_DIST:
    case POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST:
      xi = xi_opt;
      break;
  }

#if !defined(MPI)
  g_mpi.myconf = g_config.nconf;
#endif  // MPI

  // This is the start of an infinite loop

  while (1) {
    // error_sum = Sum of all the forces, energies and constraints
    double error_sum = 0.0;
    // rho_sum = Sum of density, rho, for all atoms
    double rho_sum = 0.0;

#if defined APOT && !defined MPI
    if (g_pot.format_type == POTENTIAL_FORMAT_ANALYTIC) {
      apot_check_params(xi_opt);
      update_calc_table(xi_opt, xi, 0);
    }
#endif  // APOT && !MPI

#if defined(MPI)
#if !defined(APOT)
    // exchange potential and flag value
    MPI_Bcast(xi, g_pot.calc_pot.len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif  // APOT
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (flag == 1)
      break; // Exception: flag 1 means clean up

#if defined(APOT)
    if (g_mpi.myid == 0)
      apot_check_params(xi_opt);
    MPI_Bcast(xi_opt, g_calc.ndimtot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    update_calc_table(xi_opt, xi, 0);
#else
    // if flag==2 then the potential parameters have changed -> sync
    if (flag == 2)
      potsync();
#endif  // APOT
#endif  // MPI

    // init second derivatives for splines

    // pair potentials
    //   [0, ...,  paircol - 1]
    // transfer function
    //   [paircol, ..., paircol + ntypes - 1]
    update_splines(xi, 0, g_calc.paircol + g_param.ntypes, 1);

    // embedding function
    //   [paircol + ntypes, ..., paircol + 2 * ntypes - 1]
    update_splines(xi, g_calc.paircol + g_param.ntypes, g_param.ntypes, 3);

    // dipole function
    //   [paircol + 2 * ntypes, ..., 2 * paircol + 2 * ntypes - 1]
    // quadrupole function
    //   [2 * paircol + 2 * ntypes, ..., 3 * paircol + 2 * ntypes - 1]
    update_splines(xi, g_calc.paircol + 2 * g_param.ntypes, 2 * g_calc.paircol, 1);

    // loop over configurations
    for (int config_idx = g_mpi.firstconf; config_idx < g_mpi.firstconf + g_mpi.myconf; config_idx++) {
      int uf = g_config.conf_uf[config_idx - g_mpi.firstconf];
#if defined(STRESS)
      int us = g_config.conf_us[config_idx - g_mpi.firstconf];
#endif  // STRESS
      // reset energies and stresses
      forces[g_calc.energy_p + config_idx] = 0.0;
#if defined(STRESS)
      int stress_idx = g_calc.stress_p + 6 * config_idx;
      memset(forces + stress_idx, 0, 6 * sizeof(double));
#endif  // STRESS

      // set limiting constraints
      forces[g_calc.limit_p + config_idx] = -g_config.force_0[g_calc.limit_p + config_idx];

      /* first loop over atoms: reset forces, densities */
      for (int atom_idx = 0; atom_idx < g_config.inconf[config_idx]; atom_idx++) {
        int n_i = 3 * (g_config.cnfstart[config_idx] + atom_idx);
        if (uf) {
          forces[n_i + 0] = -g_config.force_0[n_i + 0];
          forces[n_i + 1] = -g_config.force_0[n_i + 1];
          forces[n_i + 2] = -g_config.force_0[n_i + 2];
        } else {
          memset(forces + 3 * n_i, 0, 3 * sizeof(double));
        }
        // reset atomic density, dipole and quadrupol distortions
        atom_t* atom = &g_config.conf_atoms[n_i / 3 - g_mpi.firstatom];
        atom->rho = 0.0;
        atom->mu.x = 0.0;
        atom->mu.y = 0.0;
        atom->mu.z = 0.0;
        atom->lambda.xx = 0.0;
        atom->lambda.yy = 0.0;
        atom->lambda.zz = 0.0;
        atom->lambda.xy = 0.0;
        atom->lambda.yz = 0.0;
        atom->lambda.zx = 0.0;
      }

      // second loop: calculate pair forces, energies and atomic densities
      for (int atom_idx = 0; atom_idx < g_config.inconf[config_idx]; atom_idx++) {
        atom_t* atom = g_config.conf_atoms + atom_idx + g_config.cnfstart[config_idx] - g_mpi.firstatom;
        int n_i = 3 * (g_config.cnfstart[config_idx] + atom_idx);
        // loop over all neighbors
        for (int neigh_idx = 0; neigh_idx < atom->num_neigh; neigh_idx++) {
          neigh_t* neigh = atom->neigh + neigh_idx;
          // In small cells, an atom might interact with itself
          int self = (neigh->nr == atom_idx + g_config.cnfstart[config_idx]) ? 1 : 0;

          // pair potential part
          if (neigh->r < g_pot.calc_pot.end[neigh->col[0]]) {
            double phi_val = 0.0;
            double phi_grad = 0.0;
            // potential value and gradient are calculated in the same step
            if (uf)
              phi_val = splint_comb_dir(&g_pot.calc_pot, xi, neigh->slot[0], neigh->shift[0], neigh->step[0], &phi_grad);
            else
              phi_val = splint_dir(&g_pot.calc_pot, xi, neigh->slot[0], neigh->shift[0], neigh->step[0]);

            // avoid double counting if atom is interacting with itself
            if (self) {
              phi_val *= 0.5;
              phi_grad *= 0.5;
            }

            // add cohesive energy
            forces[g_calc.energy_p + config_idx] += phi_val;

            // calculate forces
            if (uf) {
              vector tmp_force;
              tmp_force.x = neigh->dist_r.x * phi_grad;
              tmp_force.y = neigh->dist_r.y * phi_grad;
              tmp_force.z = neigh->dist_r.z * phi_grad;
              forces[n_i + 0] += tmp_force.x;
              forces[n_i + 1] += tmp_force.y;
              forces[n_i + 2] += tmp_force.z;
              // actio = reactio
              forces[3 * neigh->nr + 0] -= tmp_force.x;
              forces[3 * neigh->nr + 1] -= tmp_force.y;
              forces[3 * neigh->nr + 2] -= tmp_force.z;
#if defined(STRESS)
              // also calculate pair stresses
              if (us) {
                forces[stress_idx + 0] -= neigh->dist.x * tmp_force.x;
                forces[stress_idx + 1] -= neigh->dist.y * tmp_force.y;
                forces[stress_idx + 2] -= neigh->dist.z * tmp_force.z;
                forces[stress_idx + 3] -= neigh->dist.x * tmp_force.y;
                forces[stress_idx + 4] -= neigh->dist.y * tmp_force.z;
                forces[stress_idx + 5] -= neigh->dist.z * tmp_force.x;
              }
#endif  // STRESS
            } // uf
          } // neighbor in range

          // dipole distortion part
          if (neigh->r < g_pot.calc_pot.end[neigh->col[2]]) {
            // potential value and grad are calculated in the same step
            if (uf)
              neigh->u_val = splint_comb_dir(&g_pot.calc_pot, xi, neigh->slot[2], neigh->shift[2], neigh->step[2], &neigh->u_grad);
            else
              neigh->u_val = splint_dir(&g_pot.calc_pot, xi, neigh->slot[2], neigh->shift[2], neigh->step[2]);

            // avoid double counting if atom is interacting with itself
            if (self) {
              neigh->u_val *= 0.5;
              neigh->u_grad *= 0.5;
            }

            // sum up contribution for mu
            double tmp = neigh->u_val * neigh->dist.x;
            atom->mu.x += tmp;
            g_config.conf_atoms[neigh->nr - g_mpi.firstatom].mu.x -= tmp;
            tmp = neigh->u_val * neigh->dist.y;
            atom->mu.y += tmp;
            g_config.conf_atoms[neigh->nr - g_mpi.firstatom].mu.y -= tmp;
            tmp = neigh->u_val * neigh->dist.z;
            atom->mu.z += tmp;
            g_config.conf_atoms[neigh->nr - g_mpi.firstatom].mu.z -= tmp;
          }

          // quadrupole distortion part
          if (neigh->r < g_pot.calc_pot.end[neigh->col[3]]) {
            // potential value and grad are calculated in the same step
            if (uf)
              neigh->w_val = splint_comb_dir(&g_pot.calc_pot, xi, neigh->slot[3], neigh->shift[3], neigh->step[3], &neigh->w_grad);
            else
              neigh->w_val = splint_dir(&g_pot.calc_pot, xi, neigh->slot[3], neigh->shift[3], neigh->step[3]);

            // avoid double counting if atom is interacting with itself
            if (self) {
              neigh->w_val *= 0.5;
              neigh->w_grad *= 0.5;
            }

            // sum up contribution for lambda
            // diagonal elements
            double tmp = neigh->w_val * neigh->sqrdist.xx;
            atom->lambda.xx += tmp;
            g_config.conf_atoms[neigh->nr - g_mpi.firstatom].lambda.xx += tmp;
            tmp = neigh->w_val * neigh->sqrdist.yy;
            atom->lambda.yy += tmp;
            g_config.conf_atoms[neigh->nr - g_mpi.firstatom].lambda.yy += tmp;
            tmp = neigh->w_val * neigh->sqrdist.zz;
            atom->lambda.zz += tmp;
            g_config.conf_atoms[neigh->nr - g_mpi.firstatom].lambda.zz += tmp;
            /* offdiagonal elements */
            tmp = neigh->w_val * neigh->sqrdist.yz;
            atom->lambda.yz += tmp;
            g_config.conf_atoms[neigh->nr - g_mpi.firstatom].lambda.yz += tmp;
            tmp = neigh->w_val * neigh->sqrdist.zx;
            atom->lambda.zx += tmp;
            g_config.conf_atoms[neigh->nr - g_mpi.firstatom].lambda.zx += tmp;
            tmp = neigh->w_val * neigh->sqrdist.xy;
            atom->lambda.xy += tmp;
            g_config.conf_atoms[neigh->nr - g_mpi.firstatom].lambda.xy += tmp;
          }

          // calculate atomic densities
          if (atom->type == neigh->type) {
            // then transfer(a->b)==transfer(b->a)
            if (neigh->r < g_pot.calc_pot.end[neigh->col[1]]) {
              double rho_val = splint_dir(&g_pot.calc_pot, xi, neigh->slot[1], neigh->shift[1], neigh->step[1]);
              atom->rho += rho_val;
              // avoid double counting if atom is interacting with itself
              if (!self)
                g_config.conf_atoms[neigh->nr - g_mpi.firstatom].rho += rho_val;
            }
          } else {
            // transfer(a->b)!=transfer(b->a)
            if (neigh->r < g_pot.calc_pot.end[neigh->col[1]]) {
              atom->rho += splint_dir(&g_pot.calc_pot, xi, neigh->slot[1], neigh->shift[1], neigh->step[1]);
            }
            // cannot use slot/shift to access splines
            if (neigh->r < g_pot.calc_pot.end[g_calc.paircol + atom->type])
              g_config.conf_atoms[neigh->nr - g_mpi.firstatom].rho += (*g_splint)(&g_pot.calc_pot, xi, g_calc.paircol + atom->type, neigh->r);
          }
        } // loop over neighbors

        // column of F
        int col_F = g_calc.paircol + g_param.ntypes + atom->type;
#if defined(RESCALE)
        if (atom->rho > g_pot.calc_pot.end[col_F]) {
          // then punish target function -> bad potential
          forces[g_calc.limit_p + config_idx] += DUMMY_WEIGHT * 10.0 * dsquare(atom->rho - g_pot.calc_pot.end[col_F]);
          atom->rho = g_pot.calc_pot.end[col_F];
        }

        if (atom->rho < g_pot.calc_pot.begin[col_F]) {
          // then punish target function -> bad potential
          forces[g_calc.limit_p + config_idx] += DUMMY_WEIGHT * 10.0 * dsquare(g_pot.calc_pot.begin[col_F] - atom->rho);
          atom->rho = g_pot.calc_pot.begin[col_F];
        }
#endif  // RESCALE

// embedding energy, embedding gradient
// contribution to cohesive energy is F(n)

#if !defined(RESCALE)
#if defined(APOT)
        if ((atom->rho < g_pot.calc_pot.begin[col_F]) || (atom->rho > g_pot.calc_pot.end[col_F]) || (atom->rho < 0.1)) {
          // calculate analytic value explicitly
          double temp_eng = 0.0;
          g_pot.apot_table.fvalue[col_F](atom->rho, xi_opt + g_pot.opt_pot.first[col_F], &temp_eng);
          atom->gradF = apot_gradient(atom->rho, xi_opt + g_pot.opt_pot.first[col_F], g_pot.apot_table.fvalue[col_F]);
          forces[g_calc.energy_p + config_idx] += temp_eng;
        } else {
#else // APOT
        if (atom->rho < g_pot.calc_pot.begin[col_F]) {
          // linear extrapolation left
          double rho_val = g_splint_comb(&g_pot.calc_pot, xi, col_F, g_pot.calc_pot.begin[col_F], &atom->gradF);
          forces[g_calc.energy_p + config_idx] += rho_val + (atom->rho - g_pot.calc_pot.begin[col_F]) * atom->gradF;
        } else if (atom->rho > g_pot.calc_pot.end[col_F]) {
          // and right
          double rho_val = g_splint_comb(&g_pot.calc_pot, xi, col_F, g_pot.calc_pot.end[col_F] - 0.5 * g_pot.calc_pot.step[col_F], &atom->gradF);
          forces[g_calc.energy_p + config_idx] += rho_val + (atom->rho - g_pot.calc_pot.end[col_F]) * atom->gradF;
        } else { // and in-between
#endif // APOT
          forces[g_calc.energy_p + config_idx] += g_splint_comb(&g_pot.calc_pot, xi, col_F, atom->rho, &atom->gradF);
        }
#else // RESCALE
        forces[g_calc.energy_p + config_idx] += g_splint_comb(&g_pot.calc_pot, xi, col_F, atom->rho, &atom->gradF);
#endif  // !RESCALE

        // sum up rho
        rho_sum += atom->rho;

        double eng_store = 0.0;
        /* calculate ADP energy for atom i */
        eng_store += dsquare(atom->mu.x);
        eng_store += dsquare(atom->mu.y);
        eng_store += dsquare(atom->mu.z);
        atom->nu = atom->lambda.xx + atom->lambda.yy + atom->lambda.zz;
        double trace = atom->nu / 3.0;
        eng_store += dsquare(atom->lambda.xx - trace);
        eng_store += dsquare(atom->lambda.yy - trace);
        eng_store += dsquare(atom->lambda.zz - trace);
        eng_store += dsquare(atom->lambda.xy) * 2.0;
        eng_store += dsquare(atom->lambda.yz) * 2.0;
        eng_store += dsquare(atom->lambda.zx) * 2.0;
        eng_store *= 0.5;
        forces[g_calc.energy_p + config_idx] += eng_store;
      } // second loop over atoms

      // third loop over atom: ADP forces
      // only required if we calc forces
      if (uf) {
        for (int atom_idx = 0; atom_idx < g_config.inconf[config_idx]; atom_idx++) {
          atom_t* atom = g_config.conf_atoms + atom_idx + g_config.cnfstart[config_idx] - g_mpi.firstatom;
          int n_i = 3 * (g_config.cnfstart[config_idx] + atom_idx);
          for (int neigh_idx = 0; neigh_idx < atom->num_neigh; neigh_idx++) {
            // loop over all neighbors
            neigh_t* neigh = atom->neigh + neigh_idx;
            // In small cells, an atom might interact with itself
            int self = (neigh->nr == atom_idx + g_config.cnfstart[config_idx]) ? 1 : 0;
            // column of F
            int col_F = g_calc.paircol + g_param.ntypes + atom->type;

            // are we within reach?
            if ((neigh->r < g_pot.calc_pot.end[neigh->col[1]]) || (neigh->r < g_pot.calc_pot.end[col_F - g_param.ntypes]))
            {
              double rho_grad = 0.0;
              if (neigh->r < g_pot.calc_pot.end[neigh->col[1]])
                rho_grad = splint_grad_dir(&g_pot.calc_pot, xi, neigh->slot[1],
                                        neigh->shift[1], neigh->step[1]);
              double rho_grad_j = 0.0;
              // use actio = reactio
              if (atom->type == neigh->type)
                rho_grad_j = rho_grad;
              else
                if (neigh->r < g_pot.calc_pot.end[col_F - g_param.ntypes])
                  rho_grad_j = g_splint_grad(&g_pot.calc_pot, xi, col_F - g_param.ntypes, neigh->r);

              // now we know everything - calculate forces
              double eam_force = (rho_grad * atom->gradF + rho_grad_j * g_config.conf_atoms[(neigh->nr) - g_mpi.firstatom] .gradF);

              // avoid double counting if atom is interacting with itself
              if (self)
                eam_force *= 0.5;
              vector tmp_force;
              tmp_force.x = neigh->dist_r.x * eam_force;
              tmp_force.y = neigh->dist_r.y * eam_force;
              tmp_force.z = neigh->dist_r.z * eam_force;
              forces[n_i + 0] += tmp_force.x;
              forces[n_i + 1] += tmp_force.y;
              forces[n_i + 2] += tmp_force.z;
              // actio = reactio
              forces[3 * neigh->nr + 0] -= tmp_force.x;
              forces[3 * neigh->nr + 1] -= tmp_force.y;
              forces[3 * neigh->nr + 2] -= tmp_force.z;
#if defined(STRESS)
              // and stresses
              if (us) {
                forces[stress_idx + 0] -= neigh->dist.x * tmp_force.x;
                forces[stress_idx + 1] -= neigh->dist.y * tmp_force.y;
                forces[stress_idx + 2] -= neigh->dist.z * tmp_force.z;
                forces[stress_idx + 3] -= neigh->dist.x * tmp_force.y;
                forces[stress_idx + 4] -= neigh->dist.y * tmp_force.z;
                forces[stress_idx + 5] -= neigh->dist.z * tmp_force.x;
              }
#endif          // STRESS
            } // within reach

            if (neigh->r < g_pot.calc_pot.end[neigh->col[2]]) {
              vector u_force;
              u_force.x = (atom->mu.x - g_config.conf_atoms[(neigh->nr) - g_mpi.firstatom].mu.x);
              u_force.y = (atom->mu.y - g_config.conf_atoms[(neigh->nr) - g_mpi.firstatom].mu.y);
              u_force.z = (atom->mu.z - g_config.conf_atoms[(neigh->nr) - g_mpi.firstatom].mu.z);
              // avoid double counting if atom is interacting with itself
              if (self) {
                u_force.x *= 0.5;
                u_force.y *= 0.5;
                u_force.z *= 0.5;
              }
              double tmp = SPROD(u_force, neigh->dist) * neigh->u_grad;
              vector tmp_force;
              tmp_force.x = u_force.x * neigh->u_val + tmp * neigh->dist_r.x;
              tmp_force.y = u_force.y * neigh->u_val + tmp * neigh->dist_r.y;
              tmp_force.z = u_force.z * neigh->u_val + tmp * neigh->dist_r.z;
              forces[n_i + 0] += tmp_force.x;
              forces[n_i + 1] += tmp_force.y;
              forces[n_i + 2] += tmp_force.z;
              // actio = rectio
              forces[3 * neigh->nr + 0] -= tmp_force.x;
              forces[3 * neigh->nr + 1] -= tmp_force.y;
              forces[3 * neigh->nr + 2] -= tmp_force.z;
#if defined(STRESS)
              // and stresses
              if (us) {
                forces[stress_idx + 0] -= neigh->dist.x * tmp_force.x;
                forces[stress_idx + 1] -= neigh->dist.y * tmp_force.y;
                forces[stress_idx + 2] -= neigh->dist.z * tmp_force.z;
                forces[stress_idx + 3] -= neigh->dist.x * tmp_force.y;
                forces[stress_idx + 4] -= neigh->dist.y * tmp_force.z;
                forces[stress_idx + 5] -= neigh->dist.z * tmp_force.x;
              }
#endif  // STRESS
            }

            if (neigh->r < g_pot.calc_pot.end[neigh->col[3]]) {
              sym_tens w_force;
              w_force.xx = (atom->lambda.xx + g_config.conf_atoms[(neigh->nr) - g_mpi.firstatom] .lambda.xx);
              w_force.yy = (atom->lambda.yy + g_config.conf_atoms[(neigh->nr) - g_mpi.firstatom] .lambda.yy);
              w_force.zz = (atom->lambda.zz + g_config.conf_atoms[(neigh->nr) - g_mpi.firstatom] .lambda.zz);
              w_force.yz = (atom->lambda.yz + g_config.conf_atoms[(neigh->nr) - g_mpi.firstatom] .lambda.yz);
              w_force.zx = (atom->lambda.zx + g_config.conf_atoms[(neigh->nr) - g_mpi.firstatom] .lambda.zx);
              w_force.xy = (atom->lambda.xy + g_config.conf_atoms[(neigh->nr) - g_mpi.firstatom] .lambda.xy);
              // avoid double counting if atom is interacting with itself
              if (self) {
                w_force.xx *= 0.5;
                w_force.yy *= 0.5;
                w_force.zz *= 0.5;
                w_force.yz *= 0.5;
                w_force.zx *= 0.5;
                w_force.xy *= 0.5;
              }
              vector tmp_vect;
              tmp_vect.x = w_force.xx * neigh->dist.x + w_force.xy * neigh->dist.y + w_force.zx * neigh->dist.z;
              tmp_vect.y = w_force.xy * neigh->dist.x + w_force.yy * neigh->dist.y + w_force.yz * neigh->dist.z;
              tmp_vect.z = w_force.zx * neigh->dist.x + w_force.yz * neigh->dist.y + w_force.zz * neigh->dist.z;
              double nu = (atom->nu + g_config.conf_atoms[(neigh->nr) - g_mpi.firstatom].nu) / 3.0;
              double f1 = 2.0 * neigh->w_val;
              double f2 = (SPROD(tmp_vect, neigh->dist) - nu * neigh->r * neigh->r) * neigh->w_grad - nu * f1 * neigh->r;
              vector tmp_force;
              tmp_force.x = f1 * tmp_vect.x + f2 * neigh->dist_r.x;
              tmp_force.y = f1 * tmp_vect.y + f2 * neigh->dist_r.y;
              tmp_force.z = f1 * tmp_vect.z + f2 * neigh->dist_r.z;
              forces[n_i + 0] += tmp_force.x;
              forces[n_i + 1] += tmp_force.y;
              forces[n_i + 2] += tmp_force.z;
              // actio = reactio
              forces[3 * neigh->nr + 0] -= tmp_force.x;
              forces[3 * neigh->nr + 1] -= tmp_force.y;
              forces[3 * neigh->nr + 2] -= tmp_force.z;
#if defined(STRESS)
              /* and stresses */
              if (us) {
                forces[stress_idx + 0] -= neigh->dist.x * tmp_force.x;
                forces[stress_idx + 1] -= neigh->dist.y * tmp_force.y;
                forces[stress_idx + 2] -= neigh->dist.z * tmp_force.z;
                forces[stress_idx + 3] -= neigh->dist.x * tmp_force.y;
                forces[stress_idx + 4] -= neigh->dist.y * tmp_force.z;
                forces[stress_idx + 5] -= neigh->dist.z * tmp_force.x;
              }
#endif  // STRESS
            }
          } // loop over neighbors

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
        } // third loop over atoms
      } // uf

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

      // limiting constraints per configuration
      error_sum += g_config.conf_weight[config_idx] * dsquare(forces[g_calc.limit_p + config_idx]);

    } // loop over configurations

    gather_variable(&rho_sum);

// dummy constraints (global)
#if defined(APOT)
    // add punishment for out of bounds (mostly for powell_lsq)
    if (g_mpi.myid == 0)
      error_sum += apot_punish(xi_opt, forces);
#endif  // APOT

#if !defined(NOPUNISH)
    if (g_mpi.myid == 0) {
      for (int g = 0; g < g_param.ntypes; g++) {
#if !defined(RESCALE)
        // clear field
        forces[g_calc.dummy_p + g_param.ntypes + g] = 0.0; // Free end ...
        // NEW: Constraint on U': U'(1.0)=0.0;
        forces[g_calc.dummy_p + g] = DUMMY_WEIGHT * g_splint_grad(&g_pot.calc_pot, xi, g_calc.paircol + g_param.ntypes + g, 1.0);
#else   // !RESCALE
        forces[g_calc.dummy_p + g_param.ntypes + g] = 0.0; // Free end ...
        // constraints on U`(n)
        forces[g_calc.dummy_p + g] = DUMMY_WEIGHT * g_splint_grad(&g_pot.calc_pot, xi, g_calc.paircol + g_param.ntypes + g, 0.5 * (g_pot.calc_pot.begin[g_calc.paircol + g_param.ntypes + g] + g_pot.calc_pot.end[g_calc.paircol + g_param.ntypes + g])) -
            g_config.force_0[g_calc.dummy_p + g];
#endif  // !RESCALE

        // add punishments to total error sum
        error_sum += dsquare(forces[g_calc.dummy_p + g]);
        error_sum += dsquare(forces[g_calc.dummy_p + g_param.ntypes + g]);
      } /* loop over types */

#if !defined(RESCALE)
      // NEW: Constraint on n: <n>=1.0 ONE CONSTRAINT ONLY
      // Calculate averages
      rho_sum /= (double)g_config.natoms;
      // ATTN: if there are invariant potentials, things might be problematic
      forces[g_calc.dummy_p + g_param.ntypes] = DUMMY_WEIGHT * (rho_sum - 1.0);
      error_sum += dsquare(forces[g_calc.dummy_p + g_param.ntypes]);
#endif  // !RESCALE
    }   // only root process
#endif  // !NOPUNISH

    gather_forces(&error_sum, forces);

    // root process exits this function now
    if (g_mpi.myid == 0) {
      // increase function call counter
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
