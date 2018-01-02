/****************************************************************
 *
 * force_ang_elstat.c: Routine used for calculating pair and angular
 * interactions in combination with monopole forces/energies
 *
 ****************************************************************
 *
 * Copyright 2002-2018 - the potfit development team
 *
 * ANG elstat potential: Ariel Lozano
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

#if !defined(ANG) ||  !defined(COULOMB)
#error force_ang_elstat.c compiled without ANG or COULOMB support
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
#if defined(COULOMB)
  if (g_pot.apot_table.sw_kappa)  // FIXME is sw_kappa really correct here?
    init_tails(g_pot.apot_table.dp_kappa[0]);
#endif  // COULOMB
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
  int first, col, i = flag;
  double* xi = NULL;

  /* Some useful temp variables */
  double error_sum = 0.0;
  int ne = 0;
  int size = 0;
  double charge[g_param.ntypes];
  double sum_charges;
  double dp_kappa;
  double fnval, grad, fnval_tail, grad_tail;
  int type1, type2;


  /* Temp variables */
  atom_t* atom = NULL; /* atom pointer */
  int h, j, k;
  int n_i, n_j, n_k;
  int uf;
  double angener_sum;
#if defined(STRESS)
  int us, stresses;
#endif  // STRESS

  /* Some useful temp struct variable types */
  /* neighbor pointers */
  neigh_t *neigh_j, *neigh_k;

  /* Pair variables */
  double phi_val, phi_grad;
  vector tmp_force;

  /* Angular derivatives variables */
  double dV3j, dV3k, V3, vlj, vlk, vv3j, vv3k;
  vector dfj, dfk;
  angle_t* angle;

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
    case POTENTIAL_FORMAT_KIM:
      error(1, "KIM format is not supported by ANG elstat force routine!");
      break;
  }

  ne = g_pot.apot_table.total_ne_par;
  size = g_pot.apot_table.number;

  /* This is the start of an infinite loop */
  while (1) {
    /* Reset
       error_sum = Sum of all the forces, energies and constraints
        */
    error_sum = 0.0;

#if defined(APOT) && !defined(MPI)
    if (g_pot.format_type == POTENTIAL_FORMAT_ANALYTIC) {
      apot_check_params(xi_opt);
      update_calc_table(xi_opt, xi, 0);
    }
#endif  // APOT && !MPI

#if defined(MPI)
/* exchange potential and flag value */
#if !defined(APOT)
    MPI_Bcast(xi, g_pot.calc_pot.len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif  // APOT
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (1 == flag)
      break; /* Exception: flag 1 means clean up */

#if defined(APOT)
    if (g_mpi.myid == 0)
      apot_check_params(xi_opt);
    MPI_Bcast(xi_opt, g_calc.ndimtot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    update_calc_table(xi_opt, xi, 0);
#else
    /* if flag==2 then the potential parameters have changed -> sync */
    if (flag == 2)
      potsync();
#endif  // APOT
#endif  // MPI

    /* local arrays for electrostatic parameters */
    sum_charges = 0;
    for (i = 0; i < g_param.ntypes - 1; i++) {
      if (xi_opt[2 * size + ne + i]) {
        charge[i] = xi_opt[2 * size + ne + i];
        sum_charges += g_pot.apot_table.ratio[i] * charge[i];
      } else {
        charge[i] = 0.0;
      }
    }
    g_pot.apot_table.last_charge = -sum_charges / g_pot.apot_table.ratio[g_param.ntypes - 1];
    charge[g_param.ntypes - 1] = g_pot.apot_table.last_charge;
    if (xi_opt[2 * size + ne + g_param.ntypes - 1]) {
      dp_kappa = xi_opt[2 * size + ne + g_param.ntypes - 1];
    } else {
      dp_kappa = 0.0;
    }

    /* First step is to initialize 2nd derivatives for splines */

    /* Pair potential (phi)
     * where paircol is number of pair potential columns */
    for (col = 0; col < 2 * g_calc.paircol + g_param.ntypes; col++) {
      /* Pointer to first entry */
      first = g_pot.calc_pot.first[col];

      /* Initialize 2nd derivatives
         step = width of spline knots (known as h)
         xi+first = array with spline values
         g_pot.calc_pot.last[col1] - first + 1 = num of spline pts
         *(xi + first - 2) = value of endpoint gradient (default: 1e30)
         *(xi + first - 1) = value of other endpoint gradient
         (default: phi=0.0)
         g_pot.calc_pot.d2tab + first = array to hold 2nd deriv */
      spline_ed(g_pot.calc_pot.step[col], xi + first,
                g_pot.calc_pot.last[col] - first + 1, *(xi + first - 2),
                *(xi + first - 1), g_pot.calc_pot.d2tab + first);
    }

#if !defined(MPI)
    g_mpi.myconf = g_config.nconf;
#endif  // !MPI
    /* region containing loop over configurations */
    {
      /* Loop over configurations */
      for (h = g_mpi.firstconf; h < g_mpi.firstconf + g_mpi.myconf; h++) {
        uf = g_config.conf_uf[h - g_mpi.firstconf];
#if defined(STRESS)
        us = g_config.conf_us[h - g_mpi.firstconf];
#endif  // STRESS
        /* Reset energies */
        forces[g_calc.energy_p + h] = 0.0;
#if defined(STRESS)
        /* Reset stresses */
        stresses = g_calc.stress_p + 6 * h;
        for (i = 0; i < 6; ++i)
          forces[stresses + i] = 0.0;
#endif  // STRESS

        /* FIRST LOOP: Reset forces for each atom */
        for (i = 0; i < g_config.inconf[h]; i++) {
          /* Skip every 3 spots in force array starting from position of first
           * atom */
          n_i = 3 * (g_config.cnfstart[h] + i);
          if (uf) {
            /* Set initial forces to negative of user given forces so we can
             * take difference */
            forces[n_i + 0] = -g_config.force_0[n_i + 0];
            forces[n_i + 1] = -g_config.force_0[n_i + 1];
            forces[n_i + 2] = -g_config.force_0[n_i + 2];
          } else {
            /* Set initial forces to zero if not using forces */
            forces[n_i + 0] = 0.0;
            forces[n_i + 1] = 0.0;
            forces[n_i + 2] = 0.0;
          } /* uf */
        } /* i */
        /* END OF FIRST LOOP */

        /* SECOND LOOP: Calculate pair and monopole forces and energies */
        for (i = 0; i < g_config.inconf[h]; i++) {
          /* Set pointer to temp atom pointer */
          atom = g_config.conf_atoms +
                 (g_config.cnfstart[h] - g_mpi.firstatom + i);
          /* Skip every 3 spots for force array */
          n_i = 3 * (g_config.cnfstart[h] + i);
          type1 = atom->type;
          /* Loop over neighbors */
          for (j = 0; j < atom->num_neigh; j++) {
            /* Set pointer to temp neighbor pointer */
            neigh_j = atom->neigh + j;
            type2 = neigh_j->type;


            /* calculate short-range forces */
            /* Correct column in the potential table for pair
               potential phi_ij: neigh_j->col[0] */
            /* Check neighbor lies inside pair potential's radius */
            if (neigh_j->r < g_pot.calc_pot.end[neigh_j->col[0]]) {
              /* Compute phi and phi' value given radial distance
                 NOTE: slot = spline point index right below radial distance
                 shift = % distance from 'slot' spline pt
                 step = width of spline points (given as 'h' in books)
                 0 means the pair potential columns */
              /* fn value and grad are calculated in the same step */
              if (uf)
                phi_val = splint_comb_dir(&g_pot.calc_pot, xi, neigh_j->slot[0],
                                          neigh_j->shift[0], neigh_j->step[0],
                                          &phi_grad);
              else
                phi_val = splint_dir(&g_pot.calc_pot, xi, neigh_j->slot[0],
                                     neigh_j->shift[0], neigh_j->step[0]);

              /* Add in piece contributed by neighbor to energy */
              forces[g_calc.energy_p + h] += 0.5 * phi_val;

              if (uf) {
                /* Compute tmp force values */
                tmp_force.x = neigh_j->dist_r.x * phi_grad;
                tmp_force.y = neigh_j->dist_r.y * phi_grad;
                tmp_force.z = neigh_j->dist_r.z * phi_grad;
                /* Add in force on atom i from atom j */
                forces[n_i + 0] += tmp_force.x;
                forces[n_i + 1] += tmp_force.y;
                forces[n_i + 2] += tmp_force.z;
#if defined(STRESS)
                if (us) {
                  /* also calculate pair stresses */
                  forces[stresses + 0] -= 0.5 * neigh_j->dist.x * tmp_force.x;
                  forces[stresses + 1] -= 0.5 * neigh_j->dist.y * tmp_force.y;
                  forces[stresses + 2] -= 0.5 * neigh_j->dist.z * tmp_force.z;
                  forces[stresses + 3] -= 0.5 * neigh_j->dist.x * tmp_force.y;
                  forces[stresses + 4] -= 0.5 * neigh_j->dist.y * tmp_force.z;
                  forces[stresses + 5] -= 0.5 * neigh_j->dist.z * tmp_force.x;
                }
#endif  // STRESS
              }
            }

            /* calculate monopole forces */
            /* updating tail-functions - only necessary with variing kappa */
            if (!g_pot.apot_table.sw_kappa)
#if defined(DSF)
            elstat_dsf(neigh_j->r, dp_kappa, &neigh_j->fnval_el,
                            &neigh_j->grad_el, &neigh_j->ggrad_el);
#else
            elstat_shift(neigh_j->r, dp_kappa, &neigh_j->fnval_el,
                           &neigh_j->grad_el, &neigh_j->ggrad_el);
#endif // DSF

            if (neigh_j->r < g_config.dp_cut &&
                (charge[type1] || charge[type2])) {
              fnval_tail = neigh_j->fnval_el;
              grad_tail = neigh_j->grad_el;

              fnval = charge[type1] * charge[type2] * fnval_tail;
              grad = charge[type1] * charge[type2] * grad_tail;

              forces[g_calc.energy_p + h] += 0.5 * fnval;

              if (uf) {
                tmp_force.x = 0.5 * neigh_j->dist.x * grad;
                tmp_force.y = 0.5 * neigh_j->dist.y * grad;
                tmp_force.z = 0.5 * neigh_j->dist.z * grad;
                forces[n_i + 0] += tmp_force.x;
                forces[n_i + 1] += tmp_force.y;
                forces[n_i + 2] += tmp_force.z;
                /* actio = reactio */
                n_j = 3 * neigh_j->nr;
                forces[n_j + 0] -= tmp_force.x;
                forces[n_j + 1] -= tmp_force.y;
                forces[n_j + 2] -= tmp_force.z;
#if defined(STRESS)
                /* calculate coulomb stresses */
                if (us) {
                  forces[stresses + 0] -= neigh_j->dist.x * tmp_force.x;
                  forces[stresses + 1] -= neigh_j->dist.y * tmp_force.y;
                  forces[stresses + 2] -= neigh_j->dist.z * tmp_force.z;
                  forces[stresses + 3] -= neigh_j->dist.x * tmp_force.y;
                  forces[stresses + 4] -= neigh_j->dist.y * tmp_force.z;
                  forces[stresses + 5] -= neigh_j->dist.z * tmp_force.x;
                }
#endif  // STRESS
              }
             }

            /* Compute the f_ij values and store the fn and grad in each
             * neighbor struct for easy access later */

            /* Find the correct column in the potential table for "f": f_ij
               For Binary Alloy: 0 = f_AA, 1 = f_AB, f_BA, 2 = f_BB
               where typ = A = 0 and typ = B = 1
               Note: it is "paircol" spots away in the array */

            /* Check that atom j lies inside f_col2 */
            if (neigh_j->r < g_pot.calc_pot.end[neigh_j->col[1]]) {
              /* Store the f(r_ij) value and the gradient for future use */
              neigh_j->f = splint_comb_dir(&g_pot.calc_pot, xi,
                                           neigh_j->slot[1], neigh_j->shift[1],
                                           neigh_j->step[1], &neigh_j->df);
            } else {
              /* Store f and f' = 0 if doesn't lie in boundary to be used later
               * when calculating forces */
              neigh_j->f = 0.0;
              neigh_j->df = 0.0;
            }
          } /* END LOOP OVER NEIGHBORS */

          /* Compute angular energies and forces */

          /* Find the correct column in the potential table for angle part:
             g_ijk
             Binary Alloy: 0 = g_A, 1 = g_B
             where A, B are atom type for the main atom i
             Note: it is now "2*paircol" from beginning column
             to account for
             phi(paircol)+f(paircol)
             col2 = 2 * paircol + typ1; */

          /* Loop over every angle stored with neighbors */

          /* set angl pointer to angl_part of current atom */
          angle = atom->angle_part;

          /* reset sum of angular component */
          angener_sum = 0.0;

          for (j = 0; j < atom->num_neigh - 1; j++) {
            /* Get pointer to neighbor jj */
            neigh_j = atom->neigh + j;
            /* check that j lies inside f_ij */
            if (neigh_j->r < g_pot.calc_pot.end[neigh_j->col[1]]) {

              for (k = j + 1; k < atom->num_neigh; k++) {
                /* Get pointer to neighbor k */
                neigh_k = atom->neigh + k;
                /* check that k lies inside f_ik */
                if (neigh_k->r < g_pot.calc_pot.end[neigh_k->col[1]]) {

                  /* The cos(theta) should always lie inside -1 ... 1
                     So store the g and g' without checking bounds */
                  angle->g = splint_comb_dir(&g_pot.calc_pot, xi, angle->slot,
                                             angle->shift, angle->step,
					     &angle->dg);

                  /* Sum up angular contribution for atom i caused by j and k
                     f_ij * f_ik * m_ijk */
                  angener_sum += neigh_j->f * neigh_k->f * angle->g;

                  if (uf) {
                    /* Force location for atom j */
                    n_j = 3 * neigh_j->nr;
                    /* Force location for atom k */
                    n_k = 3 * neigh_k->nr;

                    /* Some tmp variables to clean up force fn below */
                    dV3j = angle->g * neigh_j->df * neigh_k->f;
                    dV3k = angle->g * neigh_j->f * neigh_k->df;
                    V3 = neigh_j->f * neigh_k->f * angle->dg;

                    vlj = V3 * neigh_j->inv_r;
                    vlk = V3 * neigh_k->inv_r;
                    vv3j = dV3j - vlj * angle->cos;
                    vv3k = dV3k - vlk * angle->cos;

                    dfj.x = vv3j * neigh_j->dist_r.x + vlj * neigh_k->dist_r.x;
                    dfj.y = vv3j * neigh_j->dist_r.y + vlj * neigh_k->dist_r.y;
                    dfj.z = vv3j * neigh_j->dist_r.z + vlj * neigh_k->dist_r.z;

                    dfk.x = vv3k * neigh_k->dist_r.x + vlk * neigh_j->dist_r.x;
                    dfk.y = vv3k * neigh_k->dist_r.y + vlk * neigh_j->dist_r.y;
                    dfk.z = vv3k * neigh_k->dist_r.z + vlk * neigh_j->dist_r.z;

                    /* Force on atom i from j and k */
                    forces[n_i + 0] += (dfj.x + dfk.x);
                    forces[n_i + 1] += (dfj.y + dfk.y);
                    forces[n_i + 2] += (dfj.z + dfk.z);

                    /* Reaction force on atom j from i and k */
                    forces[n_j + 0] -= dfj.x;
                    forces[n_j + 1] -= dfj.y;
                    forces[n_j + 2] -= dfj.z;

                    /* Reaction force on atom k from i and j */
                    forces[n_k + 0] -= dfk.x;
                    forces[n_k + 1] -= dfk.y;
                    forces[n_k + 2] -= dfk.z;

#if defined(STRESS)
                    if (us) {
                      /* Force from j on atom i */
                      tmp_force.x = dfj.x;
                      tmp_force.y = dfj.y;
                      tmp_force.z = dfj.z;
                      forces[stresses + 0] -= neigh_j->dist.x * tmp_force.x;
                      forces[stresses + 1] -= neigh_j->dist.y * tmp_force.y;
                      forces[stresses + 2] -= neigh_j->dist.z * tmp_force.z;
                      forces[stresses + 3] -= neigh_j->dist.x * tmp_force.y;
                      forces[stresses + 4] -= neigh_j->dist.y * tmp_force.z;
                      forces[stresses + 5] -= neigh_j->dist.z * tmp_force.x;

                      /* Force from k on atom i */
                      tmp_force.x = dfk.x;
                      tmp_force.y = dfk.y;
                      tmp_force.z = dfk.z;
                      forces[stresses + 0] -= neigh_k->dist.x * tmp_force.x;
                      forces[stresses + 1] -= neigh_k->dist.y * tmp_force.y;
                      forces[stresses + 2] -= neigh_k->dist.z * tmp_force.z;
                      forces[stresses + 3] -= neigh_k->dist.x * tmp_force.y;
                      forces[stresses + 4] -= neigh_k->dist.y * tmp_force.z;
                      forces[stresses + 5] -= neigh_k->dist.z * tmp_force.x;
                    }
#endif  // STRESS
                  }     /* uf */
                  /* Increase angl pointer */
                  angle++;
                }
	      }  /* End inner loop over angles (neighbor atom k) */
	    }
          }  /* End outer loop over angles (neighbor atom j) */

	  /* add angular contribution for atom i */
          forces[g_calc.energy_p + h] += angener_sum;

        } /* END OF SECOND LOOP OVER ATOM i */


        /* THIRD loop: self energy contributions and sum-up force
         * contributions */
        double qq;
#if defined(DSF)
        double fnval_cut, gtail_cut, ggrad_cut;
        elstat_value(g_config.dp_cut, dp_kappa, &fnval_cut, &gtail_cut, &ggrad_cut);
#endif // DSF
        for (i = 0; i < g_config.inconf[h]; i++) { /* atoms */
          atom =
              g_config.conf_atoms + i + g_config.cnfstart[h] - g_mpi.firstatom;
          type1 = atom->type;
          n_i = 3 * (g_config.cnfstart[h] + i);

          /* self energy contributions */
          if (charge[type1]) {
            qq = charge[type1] * charge[type1];
#if defined(DSF)
            fnval = qq * ( DP_EPS * dp_kappa / sqrt(M_PI) +
             (fnval_cut - gtail_cut * g_config.dp_cut * g_config.dp_cut )*0.5 );
#else
            fnval = DP_EPS * dp_kappa * qq / sqrt(M_PI);
#endif // DSF
            forces[g_calc.energy_p + h] -= fnval;
          }

	  /* sum-up: whole force contributions flow into tmpsum */
          if (uf) {
#if defined(FWEIGHT)
            /* Weigh by absolute value of force */
            forces[n_i + 0] /= FORCE_EPS + atom->absforce;
            forces[n_i + 1] /= FORCE_EPS + atom->absforce;
            forces[n_i + 2] /= FORCE_EPS + atom->absforce;
#endif  // FWEIGHT
#if defined(CONTRIB)
            if (atom->contrib)
#endif  // CONTRIB
            error_sum += g_config.conf_weight[h] *
                        (dsquare(forces[n_i + 0]) + dsquare(forces[n_i + 1]) +
                         dsquare(forces[n_i + 2]));
          }
        } /* end THIRD loop over atoms */

        /* Add in the energy per atom and its weight to the sum */
        /* First divide by num atoms */
        forces[g_calc.energy_p + h] /= (double)g_config.inconf[h];

        /* Then subtract off the cohesive energy given to use by user */
        forces[g_calc.energy_p + h] -= g_config.force_0[g_calc.energy_p + h];

        /* Sum up square of this new energy term for each config
           multiplied by its respective weight */
        error_sum += g_config.conf_weight[h] * g_param.eweight *
                     dsquare(forces[g_calc.energy_p + h]);

#if defined(STRESS)
        /* LOOP OVER STRESSES */
        for (i = 0; i < 6; ++i) {
          /* Multiply weight to stresses and divide by volume */
          forces[stresses + i] /= g_config.conf_vol[h - g_mpi.firstconf];
          /* Subtract off user supplied stresses */
          forces[stresses + i] -= g_config.force_0[stresses + i];
          /* Sum in the square of each stress component with config weight */
          error_sum += g_config.conf_weight[h] * g_param.sweight *
                       dsquare(forces[stresses + i]);
        }
#endif  // STRESS
      }
    } /* END MAIN LOOP OVER CONFIGURATIONS */

/* dummy constraints (global) */
#if defined(APOT)
    /* add punishment for out of bounds (mostly for powell_lsq) */
    if (g_mpi.myid == 0) {
      error_sum += apot_punish(xi_opt, forces);
    }
#endif  // APOT

    gather_forces(&error_sum, forces);

    /* Root process only */
    if (g_mpi.myid == 0) {
      /* Increment function calls */
      g_calc.fcalls++;
      /* If total error_sum is NAN return large number instead */
      if (isnan(error_sum)) {
#if defined(DEBUG)
        printf("\n--> Force is nan! <--\n\n");
#endif  // DEBUG
        return 10e10;
      } else
        return error_sum;
    }
  } /* END OF INFINITE LOOP */

  /* Kill off other procs */
  return -1.0;
}
