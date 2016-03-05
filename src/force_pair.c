/****************************************************************
 *
 * force_pair.c: Routines used for calculating pair forces/energies
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

#if !defined(PAIR)
#error force_pair.c compiled without PAIR support
#endif

#include "potfit.h"

#include "chempot.h"
#if defined(MPI)
#include "mpi_utils.h"
#endif
#include "force.h"
#include "functions.h"
#include "potential_input.h"
#include "splines.h"
#include "utils.h"

/****************************************************************
  init_force
    called after all parameters and potentials are read
    additional assignments and initializations can be made here
****************************************************************/

void init_force(int is_worker)
{
  (void)is_worker;
  // nothing to do here for pair potentials
}

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

double calc_forces(double* xi_opt, double* forces, int flag)
{
  int first = 0;
  int col = 0;
  int i = 0;

  double* xi = NULL;
  atom_t* atom = NULL;

  int h, j;
  int n_i, n_j;
  int self;
  int uf;
#if defined(STRESS)
  int us;
  int stresses;
#endif  // STRESS

  /* pointer for neighbor table */
  neigh_t* neigh = NULL;

  /* pair variables */
  double phi_val = 0.0;
  double phi_grad = 0.0;
  vector tmp_force;

  switch (g_pot.format_type) {
    case POTENTIAL_FORMAT_UNKNOWN:
      break;
    case POTENTIAL_FORMAT_ANALYTIC:
      xi = g_pot.calc_pot.table;
      break;
    case POTENTIAL_FORMAT_TABULATED_EQ_DIST:
    case POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST:
      xi = xi_opt;
      break;
  }

  /* This is the start of an infinite loop */
  while (1) {
    double error_sum = 0.0; /* sum of squares of local process */

#if defined(APOT) && !defined(MPI)
    if (g_pot.format_type == POTENTIAL_FORMAT_ANALYTIC) {
      apot_check_params(xi_opt);
      update_calc_table(xi_opt, xi, 0);
    }
#endif  // APOT && !MPI

#if defined(MPI)
#if !defined(APOT)
    /* exchange potential and flag value */
    MPI_Bcast(xi, g_pot.calc_pot.len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif  // !APOT
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (1 == flag)
      break; /* Exception: flag 1 means clean up */

#if defined(APOT)
    if (g_mpi.myid == 0)
      apot_check_params(xi_opt);
    MPI_Bcast(xi_opt, g_calc.ndimtot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    update_calc_table(xi_opt, xi, 0);
#else   // APOT
    /* if flag==2 then the potential parameters have changed -> sync */
    if (2 == flag)
      potsync();
#endif  // APOT
#endif  // MPI

    /* init second derivatives for splines */

    /* pair potentials */
    for (col = 0; col < g_calc.paircol; col++) {
      first = g_pot.calc_pot.first[col];

      switch (g_pot.format_type) {
        case POTENTIAL_FORMAT_UNKNOWN:
          error(1, "Unknown potential format detected! (%s:%d)", __FILE__,
                __LINE__);
        case POTENTIAL_FORMAT_ANALYTIC:
        case POTENTIAL_FORMAT_TABULATED_EQ_DIST: {
          spline_ed(g_pot.calc_pot.step[col], xi + first,
                    g_pot.calc_pot.last[col] - first + 1, *(xi + first - 2),
                    0.0, g_pot.calc_pot.d2tab + first);
          break;
        }
        case POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST: {
          spline_ne(g_pot.calc_pot.xcoord + first, xi + first,
                    g_pot.calc_pot.last[col] - first + 1, *(xi + first - 2),
                    0.0, g_pot.calc_pot.d2tab + first);
        }
      }
    }

#if !defined(MPI)
    g_mpi.myconf = g_config.nconf;
#endif  // !MPI

    /* region containing loop over configurations */
    {
      /* loop over configurations */
      for (h = g_mpi.firstconf; h < g_mpi.firstconf + g_mpi.myconf; h++) {
        uf = g_config.conf_uf[h - g_mpi.firstconf];
#if defined(STRESS)
        us = g_config.conf_us[h - g_mpi.firstconf];
#endif  // STRESS
        /* reset energies and stresses */
        forces[g_calc.energy_p + h] = 0.0;
#if defined(STRESS)
        stresses = g_calc.stress_p + 6 * h;
        for (i = 0; i < 6; i++)
          forces[stresses + i] = 0.0;
#endif  // STRESS

#if defined(APOT)
        if (g_param.enable_cp)
          forces[g_calc.energy_p + h] += chemical_potential(
              g_param.ntypes, g_config.na_type[h], xi_opt + g_pot.cp_start);
#endif  // APOT

        /* first loop over atoms: reset forces, densities */
        for (i = 0; i < g_config.inconf[h]; i++) {
          if (uf) {
            n_i = 3 * (g_config.cnfstart[h] + i);
            forces[n_i + 0] = -g_config.force_0[n_i + 0];
            forces[n_i + 1] = -g_config.force_0[n_i + 1];
            forces[n_i + 2] = -g_config.force_0[n_i + 2];
          } else {
            n_i = 3 * (g_config.cnfstart[h] + i);
            forces[n_i + 0] = 0.0;
            forces[n_i + 1] = 0.0;
            forces[n_i + 2] = 0.0;
          }
        }
        /* end first loop */

        /* 2nd loop: calculate pair forces and energies */
        for (i = 0; i < g_config.inconf[h]; i++) {
          atom =
              g_config.conf_atoms + i + g_config.cnfstart[h] - g_mpi.firstatom;
          n_i = 3 * (g_config.cnfstart[h] + i);
          /* loop over neighbors */
          for (j = 0; j < atom->num_neigh; j++) {
            neigh = atom->neigh + j;
            /* In small cells, an atom might interact with itself */
            self = (neigh->nr == i + g_config.cnfstart[h]) ? 1 : 0;

            /* pair potential part */
            if (neigh->r < g_pot.calc_pot.end[neigh->col[0]]) {
              /* fn value and grad are calculated in the same step */
              if (uf)
                phi_val =
                    splint_comb_dir(&g_pot.calc_pot, xi, neigh->slot[0],
                                    neigh->shift[0], neigh->step[0], &phi_grad);
              else
                phi_val = splint_dir(&g_pot.calc_pot, xi, neigh->slot[0],
                                     neigh->shift[0], neigh->step[0]);

              /* avoid double counting if atom is interacting with a copy of
               * itself */
              if (self) {
                phi_val *= 0.5;
                phi_grad *= 0.5;
              }

              /* add cohesive energy */
              forces[g_calc.energy_p + h] += phi_val;

              /* calculate forces */
              if (uf) {
                tmp_force.x = neigh->dist_r.x * phi_grad;
                tmp_force.y = neigh->dist_r.y * phi_grad;
                tmp_force.z = neigh->dist_r.z * phi_grad;
                forces[n_i + 0] += tmp_force.x;
                forces[n_i + 1] += tmp_force.y;
                forces[n_i + 2] += tmp_force.z;
                /* actio = reactio */
                n_j = 3 * neigh->nr;
                forces[n_j + 0] -= tmp_force.x;
                forces[n_j + 1] -= tmp_force.y;
                forces[n_j + 2] -= tmp_force.z;
#if defined(STRESS)
                /* also calculate pair stresses */
                if (us) {
                  forces[stresses + 0] -= neigh->dist.x * tmp_force.x;
                  forces[stresses + 1] -= neigh->dist.y * tmp_force.y;
                  forces[stresses + 2] -= neigh->dist.z * tmp_force.z;
                  forces[stresses + 3] -= neigh->dist.x * tmp_force.y;
                  forces[stresses + 4] -= neigh->dist.y * tmp_force.z;
                  forces[stresses + 5] -= neigh->dist.z * tmp_force.x;
                }
#endif  // STRESS
              }
            } /* neighbor in range */
          }   /* loop over all neighbors */

          /* then we can calculate contribution of forces right away */
          if (uf) {
#if defined(FWEIGHT)
            /* Weigh by absolute value of force */
            forces[n_i + 0] /= FORCE_EPS + atom->absforce;
            forces[n_i + 1] /= FORCE_EPS + atom->absforce;
            forces[n_i + 2] /= FORCE_EPS + atom->absforce;
#endif  // FWEIGHT

/* sum up forces */
#if defined(CONTRIB)
            if (atom->contrib)
#endif  // CONTRIB
              error_sum += g_config.conf_weight[h] * (dsquare(forces[n_i + 0]) +
                                                      dsquare(forces[n_i + 1]) +
                                                      dsquare(forces[n_i + 2]));
          }
        } /* second loop over atoms */

        /* energy contributions */
        forces[g_calc.energy_p + h] /= (double)g_config.inconf[h];
        forces[g_calc.energy_p + h] -= g_config.force_0[g_calc.energy_p + h];
        error_sum += g_config.conf_weight[h] * g_param.eweight *
                     dsquare(forces[g_calc.energy_p + h]);

#if defined(STRESS)
        /* stress contributions */
        if (uf && us) {
          for (i = 0; i < 6; i++) {
            forces[stresses + i] /= g_config.conf_vol[h - g_mpi.firstconf];
            forces[stresses + i] -= g_config.force_0[stresses + i];
            error_sum += g_config.conf_weight[h] * g_param.sweight *
                         dsquare(forces[stresses + i]);
          }
        }
#endif  // STRESS

      } /* loop over configurations */
    }   /* parallel region */

/* dummy constraints (global) */
#if defined(APOT)
    /* add punishment for out of bounds (mostly for powell_lsq) */
    if (g_mpi.myid == 0) {
      error_sum += apot_punish(xi_opt, forces);
    }
#endif  // APOT

    gather_forces(&error_sum, forces);

    /* root process exits this function now */
    if (g_mpi.myid == 0) {
      g_calc.fcalls++; /* Increase function call counter */
      if (isnan(error_sum)) {
#if defined(DEBUG)
        printf("\n--> Force is nan! <--\n\n");
#endif  // DEBUG
        return 10e10;
      } else
        return error_sum;
    }
  } /* end of infinite loop */

  /* once a non-root process arrives here, all is done. */
  return -1.0;
}
