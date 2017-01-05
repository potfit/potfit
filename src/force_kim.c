/****************************************************************
 *
 * force_kim.c: Routines used for calculating forces/energies via KIM
 *
 ****************************************************************
 *
 * Copyright 2002-2017 - the potfit development team
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

#include "potfit.h"

#include "chempot.h"
#include "kim.h"
#include "utils.h"
#include "functions.h"

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
  /* local variables */
  double tmpsum = 0.0, sum = 0.0;
  int   h;
  int   n_i;
  int   uf;
#if defined(STRESS)
  int   us, stresses;
#endif /* STRESS */
#if defined(FWEIGHT) || defined(CONTRIB)
  atom_t *atom;
#endif // FWEIGHT || CONTRIB
  int status;
  double* kimenergy;
  double* kimforce;
#if defined(STRESS)
  double* kimvirial;
#endif
  double weight;

  // publish KIM parameters
  for (int i = 0; i < g_config.nconf; i++) {
    int idx = 0;
    for (int j = 0; j < g_kim.num_opt_param; ++j) {
      double* data = ((double**)g_kim.param_value[i])[j];
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
    tmpsum = 0.0; // sum of squares of local process

#if !defined(MPI)
    g_mpi.myconf = g_config.nconf;
#endif /* MPI */

    // region containing loop over configurations
    {
      // loop over configurations
      for (h = g_mpi.firstconf; h < g_mpi.firstconf + g_mpi.myconf; h++) {
        uf = g_config.conf_uf[h - g_mpi.firstconf];
#if defined(STRESS)
        us = g_config.conf_us[h - g_mpi.firstconf];
#endif /* STRESS */

/* Calculate forces from KIM (general forces, including forces, virial and energy) */

#if defined(STRESS)
          /* get data */
          KIM_API_getm_data(g_kim.pkim[h], &status, 3*3,
                            "energy", &kimenergy, 1,
                            "forces", &kimforce,  uf,
                            "virial", &kimvirial, us);
          if (KIM_STATUS_OK > status)
            error(1, "");
#else
          /* get data */
          KIM_API_getm_data(g_kim.pkim[h], &status, 2*3,
                            "energy", &kimenergy, 1,
                            "forces", &kimforce,  uf);
          if (KIM_STATUS_OK > status)
            error(1, "");
#endif

          /* Call model compute */
          status = KIM_API_model_compute(g_kim.pkim[h]);
          if (KIM_STATUS_OK > status)
            error(1, "");

        /* forces contributation */
        weight = sqrt(g_config.conf_weight[h]);
        for (int i = 0; i < g_config.inconf[h]; i++) {
#if defined(FWEIGHT) || defined(CONTRIB)
          atom = g_config.conf_atoms + i + g_config.cnfstart[h] - g_mpi.firstatom;
#endif // FWEIGHT || CONTRIB
          n_i = DIM * (g_config.cnfstart[h] + i);
          if (uf) {
            forces[n_i + 0] = weight * (kimforce[DIM * i + 0] - g_config.force_0[n_i + 0]);
            forces[n_i + 1] = weight * (kimforce[DIM * i + 1] - g_config.force_0[n_i + 1]);
            forces[n_i + 2] = weight * (kimforce[DIM * i + 2] - g_config.force_0[n_i + 2]);
          } else {
            forces[n_i + 0] = 0.0;
            forces[n_i + 1] = 0.0;
            forces[n_i + 2] = 0.0;
          }

#if defined(FWEIGHT)
          /* weigh by absolute value of force */
          forces[n_i + 0] /= FORCE_EPS + atom->absforce;
          forces[n_i + 1] /= FORCE_EPS + atom->absforce;
          forces[n_i + 2] /= FORCE_EPS + atom->absforce;
#endif /* FWEIGHT */

#if defined(CONTRIB)
          if (atom->contrib)
#endif /* CONTRIB */
            tmpsum += (dsquare(forces[n_i + 0])
                    + dsquare(forces[n_i + 1])
                    + dsquare(forces[n_i + 2]) );

        }

        /* reset energies and stresses */
        forces[g_calc.energy_p + h] = 0.0;
#if defined(STRESS)
        stresses = g_calc.stress_p + 6 * h;
        for (int i = 0; i < 6; i++)
          forces[stresses + i] = 0.0;
#endif /* STRESS */

        /* energy contributation */
#if defined(APOT)
        if (g_param.enable_cp)
          forces[g_calc.energy_p + h] = chemical_potential(g_param.ntypes, g_config.na_type[h], xi_opt + g_pot.cp_start);
#endif /* APOT */

		    weight = sqrt(g_config.conf_weight[h] * g_param.eweight);
        forces[g_calc.energy_p + h] = weight * (*kimenergy) / (double)g_config.inconf[h];
        forces[g_calc.energy_p + h] -=  weight * g_config.force_0[g_calc.energy_p + h];
        tmpsum += dsquare(forces[g_calc.energy_p + h]);


#if defined(STRESS)
        /* stress contributions */
        if (uf && us) {
		      weight = sqrt(g_config.conf_weight[h] * g_param.sweight);
          for (int i = 0; i < 6; i++) {
            forces[stresses + i]  = weight*kimvirial[i];
            forces[stresses + i] /= g_config.conf_vol[h - g_mpi.firstconf];
            forces[stresses + i] -= weight*g_config.force_0[stresses + i];
            tmpsum += dsquare(forces[stresses + i]);
          }
        }
#endif /* STRESS */

      }       /* loop over configurations */
    }       /* parallel region */

#if defined(MPI)
    /* reduce global sum */
    sum = 0.0;
    MPI_Reduce(&tmpsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    /* gather forces, energies, stresses */
    if (0 == g_mpi.myid) {    /* root node already has data in place */
      /* forces */
      MPI_Gatherv(MPI_IN_PLACE, myatoms, MPI_VECTOR, forces,
          atom_len, atom_dist, MPI_VECTOR, 0, MPI_COMM_WORLD);
      /* energies */
      MPI_Gatherv(MPI_IN_PLACE, g_mpi.myconf, MPI_DOUBLE, forces + g_calc.energy_p,
          conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if defined(STRESS)
      /* stresses */
      MPI_Gatherv(MPI_IN_PLACE, g_mpi.myconf, MPI_STENS, forces + g_calc.stress_p,
          conf_len, conf_dist, MPI_STENS, 0, MPI_COMM_WORLD);
#endif /* STRESS */
    } else {
      /* forces */
      MPI_Gatherv(forces + g_mpi.firstatom * 3, myatoms, MPI_VECTOR,
          forces, atom_len, atom_dist, MPI_VECTOR, 0, MPI_COMM_WORLD);
      /* energies */
      MPI_Gatherv(forces + g_calc.energy_p + g_mpi.firstconf, g_mpi.myconf, MPI_DOUBLE,
          forces + g_calc.energy_p, conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if defined(STRESS)
      /* stresses */
      MPI_Gatherv(forces + g_calc.stress_p + 6 * g_mpi.firstconf, g_mpi.myconf, MPI_STENS,
          forces + g_calc.stress_p, conf_len, conf_dist, MPI_STENS, 0, MPI_COMM_WORLD);
#endif /* STRESS */
    }
#else
    sum = tmpsum;   /* global sum = local sum  */
#endif /* MPI */

    /* root process exits this function now */
    if (0 == g_mpi.myid) {
      g_calc.fcalls++;     /* Increase function call counter */
      if (isnan(sum)) {
#if defined(DEBUG)
        printf("\n--> Force is nan! <--\n\n");
#endif /* DEBUG */
        return 10e10;
      } else
        return sum;
    }
  }       /* end of infinite loop */

  /* once a non-root process arrives here, all is done. */
  return -1.0;
}
