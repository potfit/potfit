/*******************************************************************************
 *
 * force_kim.c: Routines used for calculating forces/energies via KIM
 *
 *******************************************************************************
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

#include "potfit.h"
#include "kim.h"


#include "force.h"
#include "functions.h"
#include "memory.h"
#if defined(MPI)
#include "mpi_utils.h"
#endif
#include "potential_input.h"
#include "potential_output.h"
#include "utils.h"


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
 *     opt_pot.table - part of the struct opt_pot, but it can also be
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
  /* local variables */
  int h;
  int n_i;
  int uf;
#if defined(STRESS)
  int us, stresses;
#endif // STRESS 
  atom_t* atom;
  int status; 
  double* kimenergy;
  double* kimforce;       
  double* kimvirial;  
  int kim_us = 0;
  int i;

  /* publish KIM parameters */
  for (i = 0; i < nconf; i++) {
    status = publish_param(pkimObj[i], &FreeParamAllConfig[i], xi_opt);  
    if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__, "publish_parameters", status);
      exit(1);
    }
  }

  /*
    KIM_API_print(pkimObj[0],&status);
  */

  /* This is the start of an infinite loop */
  while (1) {
    
#if !defined(MPI)
    g_mpi.myconf = g_config.nconf;
#endif  // MPI
    
    /* region containing loop over configurations */
    {
      /* loop over configurations */
      for (h = g_mpi.firstconf; h < g_mpi.firstconf + g_mpi.myconf; h++) {
        uf = g_config.conf_uf[h - g_mpi.firstconf];
#if defined(STRESS)
        us = g_config.conf_us[h - g_mpi.firstconf];
#endif // STRESS 

#if defined(APOT)
        if (g_param.enable_cp)
          forces[g_calc.energy_p + h] += chemical_potential(
              g_param.ntypes, g_config.na_type[h], xi_opt + g_pot.cp_start);
#endif  // APOT

/* Calculate forces from KIM (general forces, including forces, virial and energy) */
#if defined(STRESS) 
        kim_us = us;
#endif // STRESS


        /* compute forces */
        status = calc_force_KIM(pkimObj[h], &kimenergy,  &kimforce,  &kimvirial, uf, kim_us);
        if (KIM_STATUS_OK > status) {
          KIM_API_report_error(__LINE__, __FILE__, "KIM: compute forces failed", status);
	  error(1,"KIM Status error in calc_force_KIM");
        }
        

/* 
double* coords;
KIM_API_getm_data(pkimObj[h], &status, 1*3,
                              "coordinates",         &coords,              1);
        if (KIM_STATUS_OK > status) {
              KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);
              return status;
            }



       FILE* fp = fopen("potfit_force", "w");
       for (i = 0; i < inconf[h]; i++){
           fprintf(fp, "%18.10e %18.10e %18.10e %18.10e %18.10e %18.10e\n",
             coords[DIM*i+0], coords[DIM*i+1], coords[DIM*i+2], kimforce[DIM*i+0], 
	     kimforce[DIM*i+1], kimforce[DIM*i+2]);
         }
         fflush(fp);
         fclose(fp);

         fp = fopen("potfit_energy", "w");
         fprintf(fp, "total energy%18.10e\n", *kimenergy);
         fprintf(fp, "average energy%18.10e", *kimenergy/(double)inconf[h]);
         fflush(fp);
         fclose(fp);

         printf("finished job, let's exit. \n");
         exit(1);


*/



        /* forces contribution */
        for (i = 0; i < g_config.inconf[h]; i++) {
          atom = conf_atoms + i + cnfstart[h] - firstatom;
          n_i = DIM * (g_config.cnfstart[h] + i);  
          if (uf) {
            forces[n_i + 0] = kimforce[DIM*i + 0] - g_config.force_0[n_i + 0];
            forces[n_i + 1] = kimforce[DIM*i + 1] - g_config.force_0[n_i + 1];
            forces[n_i + 2] = kimforce[DIM*i + 2] - g_config.force_0[n_i + 2];
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
#endif // FWEIGHT 

#if defined(CONTRIB)
          if (atom->contrib) 
#endif // CONTRIB 
	    error_sum += g_config.conf_weight[h] * (dsquare(forces[n_i + 0]) +
						    dsquare(forces[n_i + 1]) +
						    dsquare(forces[n_i + 2]));
	  
        }

        /* reset energies and stresses */
        forces[g_calc.energy_p + h] = 0.0;
#if defined(STRESS)
        stresses = g_calc.stress_p + 6 * h;
        for (i = 0; i < 6; i++)
          forces[stresses + i] = 0.0;
#endif // STRESS 

        /* energy contributation */
#if defined(APOT)
        if (g_param.enable_cp)
          forces[g_calc.energy_p + h] = chemical_potential(ntypes, na_type[h], xi_opt + cp_start);
#endif // APOT 

        forces[g_calc.energy_p + h] = *kimenergy / (double)inconf[h];
        forces[g_calc.energy_p + h] -= g_config.force_0[g_calc.energy_p + h];
        error_sum += g_config.conf_weight[h] * g_param.eweight *
                     dsquare(forces[g_calc.energy_p + h]);

#if defined(STRESS)
        /* stress contributions */
        if (uf && us) {
          for (i = 0; i < 6; i++) {
            forces[stresses + i]  = kimvirial[i]; 
            forces[stresses + i] /= g_config.conf_vol[h - firstconf];
            forces[stresses + i] -= g_config.force_0[stresses + i];
            error_sum += g_config.conf_weight[h] * g_param.sweight *
                         dsquare(forces[stresses + i]);
          } 
        } 
#endif // STRESS 

      }       /* loop over configurations */
    }       /* parallel region */


    /* dummy constraints (global) */
#if defined(APOT)
    /* add punishment for out of bounds (mostly for powell_lsq) */
    if (g_mpi.myid == 0) {
      error_sum += apot_punish(xi_opt, forces);
    }
#endif // APOT 

    gather_forces(&error_sum, forces);
    
    /* root process exits this function now */
    if (g_mpi.myid == 0) {
      g_calc.fcalls++; /* Increase function call counter */
      if (isnan(error_sum)) {
#if defined(DEBUG)
        printf("\n--> Force is nan! <--\n\n");
#endif // DEBUG 
        return 10e10;
      } else
        return error_sum;
    }
  }       /* end of infinite loop */

  /* once a non-root process arrives here, all is done. */
  return -1.0;
}
