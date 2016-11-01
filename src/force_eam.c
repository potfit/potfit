/****************************************************************
 *
 * force_eam.c: Routine used for calculating eam forces/energies
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

#if !defined(EAM)
#error force_eam.c compiled without EAM support
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
 *     g_pot.opt_pot.table - part of the struct g_pot.opt_pot, but it can also
 *be
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
    case POTENTIAL_FORMAT_KIM:
      error(1, "KIM format is not supported by EAM force routine!");
      break;
  }

#if !defined(MPI)
  g_mpi.myconf = g_config.nconf;
#endif  // MPI

  // This is the start of an infinite loop
  while (1) {
    // sum of squares of local process
    double error_sum = 0.0;
    double rho_sum = 0.0;
#if defined(TBEAM)
    double rho_s_sum = 0.0;
#endif  // TBEAM

#if defined(APOT) && !defined(MPI)
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
#else   // APOT
    // if flag == 2 then the potential parameters have changed -> sync
    if (flag == 2)
      potsync();
#endif  // APOT
#endif  // MPI

    // init second derivatives for splines

    // pair potential
    //   [0, ...,  paircol - 1]
    // transfer function
    //   [g_calc.paircol, ..., g_calc.paircol + g_param.ntypes - 1]
    update_splines(xi, 0, g_calc.paircol + g_param.ntypes, 1);

    // embedding function
    //   [g_calc.paircol + g_param.ntypes, ..., g_calc.paircol + 2 * g_param.ntypes - 1]
    update_splines(xi, g_calc.paircol + g_param.ntypes, g_param.ntypes, 3);

#if defined(TBEAM)
    // s-band transfer function
    //   [g_calc.paircol + 2 * g_param.ntypes, ..., g_calc.paircol + 3 * g_param.ntypes - 1]
    update_splines(xi, g_calc.paircol + 2 * g_param.ntypes, g_param.ntypes, 1);

    // s-band embedding function
    //   [g_calc.paircol + 3 * g_param.ntypes, ..., g_calc.paircol + 4 * g_param.ntypes - 1]
    update_splines(xi, g_calc.paircol + 3 * g_param.ntypes, g_param.ntypes, 3);
#endif  // TBEAM

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

#if defined(RESCALE)
      // set limiting constraints
      forces[g_calc.limit_p + config_idx] = -g_config.force_0[g_calc.limit_p + config_idx];
#endif  // RESCALE

      // first loop: reset forces and densities
      for (int atom_idx = 0; atom_idx < g_config.inconf[config_idx]; atom_idx++) {
        int n_i = g_config.cnfstart[config_idx] + atom_idx;
        if (uf) {
          forces[3 * n_i + 0] = -g_config.force_0[3 * n_i + 0];
          forces[3 * n_i + 1] = -g_config.force_0[3 * n_i + 1];
          forces[3 * n_i + 2] = -g_config.force_0[3 * n_i + 2];
        } else {
          memset(forces + 3 * n_i, 0, 3 * sizeof(double));
        }
        // reset atomic density
        g_config.conf_atoms[n_i - g_mpi.firstatom].rho = 0.0;
#if defined(TBEAM)
        g_config.conf_atoms[n_i - g_mpi.firstatom].rho_s = 0.0;
#endif  // TBEAM
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
#endif // STRESS
            } // uf
          } // neighbor in range

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
#if defined(TBEAM)
            if (neigh->r < g_pot.calc_pot.end[neigh->col[2]]) {
              double rho_s_val = splint_dir(&g_pot.calc_pot, xi, neigh->slot[2], neigh->shift[2], neigh->step[2]);
              atom->rho_s += rho_s_val;
              // avoid double counting if atom is interacting with itself
              if (!self)
                g_config.conf_atoms[neigh->nr - g_mpi.firstatom].rho_s += rho_s_val;
            }
#endif  // TBEAM
          } else {
            // transfer(a->b)!=transfer(b->a)
            if (neigh->r < g_pot.calc_pot.end[neigh->col[1]])
              atom->rho += splint_dir(&g_pot.calc_pot, xi, neigh->slot[1], neigh->shift[1], neigh->step[1]);
            // cannot use slot/shift to access splines
            if (neigh->r < g_pot.calc_pot.end[g_calc.paircol + atom->type])
              g_config.conf_atoms[neigh->nr - g_mpi.firstatom].rho +=
                g_splint(&g_pot.calc_pot, xi, g_calc.paircol + atom->type, neigh->r);
#if defined(TBEAM)
            if (neigh->r < g_pot.calc_pot.end[neigh->col[2]])
              atom->rho_s += splint_dir(&g_pot.calc_pot, xi, neigh->slot[2], neigh->shift[2], neigh->step[2]);
            if (neigh->r < g_pot.calc_pot .end[g_calc.paircol + 2 * g_param.ntypes + atom->type])
              g_config.conf_atoms[neigh->nr - g_mpi.firstatom].rho_s += g_splint(&g_pot.calc_pot, xi, g_calc.paircol + 2 * g_param.ntypes + atom->type, neigh->r);
#endif  // TBEAM
          }
        } // loop over all neighbors

        // column of F
        int col_F = g_calc.paircol + g_param.ntypes + atom->type;
#if defined(TBEAM)
        // column of F of the s-band
        int col_F_s = col_F + 2 * g_param.ntypes;
#endif  // TBEAM

#if defined(RESCALE)
        // we punish the potential for bad behavior:
        // if the density of one atom is smaller or greater than we have the
        // embedding function tabulated a punishment is added

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
#if defined(TBEAM)
        if (atom->rho_s > g_pot.calc_pot.end[col_F_s]) {
          // then punish target function -> bad potential
          forces[g_calc.limit_p + config_idx] += DUMMY_WEIGHT * 10.0 * dsquare(atom->rho_s - g_pot.calc_pot.end[col_F_s]);
          atom->rho_s = g_pot.calc_pot.end[col_F_s];
        }

        if (atom->rho_s < g_pot.calc_pot.begin[col_F_s]) {
          // then punish target function -> bad potential
          forces[g_calc.limit_p + config_idx] += DUMMY_WEIGHT * 10.0 * dsquare(g_pot.calc_pot.begin[col_F_s] - atom->rho_s);
          atom->rho_s = g_pot.calc_pot.begin[col_F_s];
        }
#endif  // TBEAM
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

#if defined(TBEAM)
#if !defined(RESCALE)
#if defined(APOT)
        if ((atom->rho_s < g_pot.calc_pot.begin[col_F_s]) || (atom->rho_s > g_pot.calc_pot.end[col_F_s]) || (atom->rho_s < 0.1)) {
          // calculate analytic value explicitly
          double temp_eng = 0.0;
          g_pot.apot_table.fvalue[col_F_s](atom->rho_s, xi_opt + g_pot.opt_pot.first[col_F_s], &temp_eng);
          atom->gradF_s = apot_gradient(atom->rho_s, xi_opt + g_pot.opt_pot.first[col_F_s], g_pot.apot_table.fvalue[col_F_s]);
          forces[g_calc.energy_p + config_idx] += temp_eng;
        } else {
#else // APOT
        if (atom->rho_s < g_pot.calc_pot.begin[col_F_s]) {
          // linear extrapolation left
          double rho_s_val = g_splint_comb(&g_pot.calc_pot, xi, col_F_s, g_pot.calc_pot.begin[col_F_s], &atom->gradF_s);
          forces[g_calc.energy_p + config_idx] += rho_s_val + (atom->rho_s - g_pot.calc_pot.begin[col_F_s]) * atom->gradF_s;
        } else if (atom->rho_s > g_pot.calc_pot.end[col_F_s]) {
          // and right
          double rho_s_val = g_splint_comb(&g_pot.calc_pot, xi, col_F_s, g_pot.calc_pot.end[col_F_s] - 0.5 * g_pot.calc_pot.step[col_F_s], &atom->gradF_s);
          forces[g_calc.energy_p + config_idx] += rho_s_val + (atom->rho_s - g_pot.calc_pot.end[col_F_s]) * atom->gradF_s;
        } else {
#endif // APOT
          forces[g_calc.energy_p + config_idx] += g_splint_comb(&g_pot.calc_pot, xi, col_F_s, atom->rho_s, &atom->gradF_s);
        }
#else // RESCALE
        forces[g_calc.energy_p + config_idx] += g_splint_comb(&g_pot.calc_pot, xi, col_F_s, atom->rho_s, &atom->gradF_s);
#endif // RESCALE

        // sum up rho_s
        rho_s_sum += atom->rho_s;
#endif // TBEAM

      } // second loop

      // third loop: EAM force
      // only required if we calculate forces
      if (uf) {
        for (int atom_idx = 0; atom_idx < g_config.inconf[config_idx]; atom_idx++) {
          atom_t* atom = g_config.conf_atoms + atom_idx + g_config.cnfstart[config_idx] -
                  g_mpi.firstatom;
          int n_i = 3 * (g_config.cnfstart[config_idx] + atom_idx);
          // loop over all neighbors
          for (int neigh_idx = 0; neigh_idx < atom->num_neigh; neigh_idx++) {
            neigh_t* neigh = atom->neigh + neigh_idx;
            // In small cells, an atom might interact with itself
            int self = (neigh->nr == atom_idx + g_config.cnfstart[config_idx]) ? 1 : 0;
            // column of F
            int col_F = g_calc.paircol + g_param.ntypes + atom->type;
#if defined(TBEAM)
            int col_F_s = col_F + 2 * g_param.ntypes;
#endif  // TBEAM
            double r = neigh->r;
            // are we within reach?
            if ((r < g_pot.calc_pot.end[neigh->col[1]]) || (r < g_pot.calc_pot.end[col_F - g_param.ntypes])) {
              double rho_grad = 0.0;
              if (r < g_pot.calc_pot.end[neigh->col[1]])
                rho_grad = splint_grad_dir(&g_pot.calc_pot, xi, neigh->slot[1], neigh->shift[1], neigh->step[1]);
              // use actio = reactio
              double rho_grad_j = 0.0;
              if (atom->type == neigh->type)
                rho_grad_j = rho_grad;
              else if (r < g_pot.calc_pot.end[col_F - g_param.ntypes])
                rho_grad_j =  g_splint_grad(&g_pot.calc_pot, xi, col_F - g_param.ntypes, r);
              // now we know everything - calculate forces
              double eam_force = (rho_grad * atom->gradF + rho_grad_j * g_config.conf_atoms[(neigh->nr) - g_mpi.firstatom] .gradF);

#if defined(TBEAM)
              // s-band contribution to force for TBEAM
              if ((r < g_pot.calc_pot.end[neigh->col[2]]) || (r < g_pot.calc_pot.end[col_F_s - g_param.ntypes])) {
                double rho_s_grad = 0.0;
                if (r < g_pot.calc_pot.end[neigh->col[2]])
                  rho_s_grad = splint_grad_dir(&g_pot.calc_pot, xi, neigh->slot[2], neigh->shift[2], neigh->step[2]);
                // use actio = reactio
                double rho_s_grad_j = 0.0;
                if (atom->type == neigh->type)
                  rho_s_grad_j = rho_s_grad;
                else if (r < g_pot.calc_pot.end[col_F_s - g_param.ntypes])
                  rho_s_grad_j = g_splint_grad(&g_pot.calc_pot, xi, col_F_s - g_param.ntypes, r);
                // now we know everything - calculate forces
                eam_force += (rho_s_grad * atom->gradF_s + rho_s_grad_j * g_config.conf_atoms[(neigh->nr) - g_mpi.firstatom] .gradF_s);
              }
#endif  // TBEAM

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
          }   // loop over neighbours

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
      } // use forces

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
          error_sum += g_config.conf_weight[config_idx] * g_param.sweight *
                    dsquare(forces[stress_idx + i]);
        }
      }
#endif  // STRESS

#if defined(RESCALE)
      // limiting constraints per configuration
      error_sum += g_config.conf_weight[config_idx] * dsquare(forces[g_calc.limit_p + config_idx]);
#endif  // RESCALE
    } // loop over configurations

    // dummy constraints (global)
#if defined(APOT)
    // add punishment for out of bounds (mostly for powell_lsq)
    if (g_mpi.myid == 0)
      error_sum += apot_punish(xi_opt, forces);
#endif  // APOT

    gather_variable(&rho_sum);
#if defined(TBEAM)
    gather_variable(&rho_s_sum);
#endif // TBEAM

    // rho_sum and rho_s_sum are now correct on root node

#if !defined(NOPUNISH)
    if (g_mpi.myid == 0) {
      for (int g = 0; g < g_param.ntypes; g++) {
#if !defined(RESCALE)
        // clear field
        forces[g_calc.dummy_p + g_param.ntypes + g] = 0.0; // Free end ...
        // NEW: Constraint on U': U'(1.0)=0.0;
        forces[g_calc.dummy_p + g] = DUMMY_WEIGHT * g_splint_grad(&g_pot.calc_pot, xi, g_calc.paircol + g_param.ntypes + g, 1.0);
#if defined(TBEAM)
        // clear field
        forces[g_calc.dummy_p + 3 * g_param.ntypes + g] = 0.0; // Free end ...
        // NEW: Constraint on U': U'(1.0)=0.0;
        forces[g_calc.dummy_p + 2 * g_param.ntypes + g] = DUMMY_WEIGHT * g_splint_grad(&g_pot.calc_pot, xi, g_calc.paircol + 3 * g_param.ntypes + g, 1.0);
#endif  // TBEAM
#else   // !RESCALE
        forces[g_calc.dummy_p + g_param.ntypes + g] = 0.0; // Free end ...
        // constraints on U`(n)
        forces[g_calc.dummy_p + g] = DUMMY_WEIGHT * g_splint_grad(&g_pot.calc_pot, xi, g_calc.paircol + g_param.ntypes + g, 0.5 * (g_pot.calc_pot .begin[g_calc.paircol + g_param.ntypes + g] + g_pot.calc_pot .end[g_calc.paircol + g_param.ntypes + g])) - g_config.force_0[g_calc.dummy_p + g];
#if defined(TBEAM)
        forces[g_calc.dummy_p + 3 * g_param.ntypes + g] = 0.0; // Free end ...
        // constraints on U`(n)
        forces[g_calc.dummy_p + 2 * g_param.ntypes + g] = DUMMY_WEIGHT * g_splint_grad(&g_pot.calc_pot, xi, g_calc.paircol + 3 * g_param.ntypes + g, 0.5 * (g_pot.calc_pot .begin[g_calc.paircol + 3 * g_param.ntypes + g] + g_pot.calc_pot .end[g_calc.paircol + 3 * g_param.ntypes + g])) - g_config.force_0[g_calc.dummy_p + 2 * g_param.ntypes + g];
#endif  // TBEAM
#endif  // !RESCALE

        // add punishments to total error sum
        error_sum += dsquare(forces[g_calc.dummy_p + g]);
        error_sum += dsquare(forces[g_calc.dummy_p + g_param.ntypes + g]);
#if defined(TBEAM)
        error_sum += dsquare(forces[g_calc.dummy_p + 2 * g_param.ntypes + g]);
        error_sum += dsquare(forces[g_calc.dummy_p + 3 * g_param.ntypes + g]);
#endif  // TBEAM
      } // loop over types

#if !defined(RESCALE)
      // NEW: Constraint on n: <n>=1.0 ONE CONSTRAINT ONLY
      if (rho_sum > 0.0) {
        // Calculate averages
        rho_sum /= (double)g_config.natoms;
        // ATTN: if there are invariant potentials, things might be problematic
        forces[g_calc.dummy_p + g_param.ntypes] = DUMMY_WEIGHT * (rho_sum - 1.0);
        error_sum += dsquare(forces[g_calc.dummy_p + g_param.ntypes]);
      }
#if defined(TBEAM)
      if (rho_s_sum > 0.0) {
        // Calculate averages
        rho_s_sum /= (double)g_config.natoms;
        // ATTN: if there are invariant potentials, things might be problematic
        forces[g_calc.dummy_p + 3 * g_param.ntypes] = DUMMY_WEIGHT * (rho_s_sum - 1.0);
        error_sum += dsquare(forces[g_calc.dummy_p + 3 * g_param.ntypes]);
      }
#endif  // TBEAM
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

  // once a non-root process arrives here, all is done.
  return -1.0;
}
