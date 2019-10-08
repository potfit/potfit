/****************************************************************
 *
 * potfit.c: Contains main potfit program
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

#include <time.h>

#include "potfit.h"

#include "config.h"
#include "errors.h"
#include "force.h"
#include "functions.h"
#include "kim.h"
#include "memory.h"
#include "mpi_utils.h"
#include "optimize.h"
#include "params.h"
#include "potential_input.h"
#include "potential_output.h"
#include "random.h"
#include "utils.h"
#include "uq.h"

// forward declarations of helper functions

void read_input_files(int argc, char** argv);
void start_mpi_worker(double* force);

// potfit global variables

potfit_calculation g_calc;
potfit_configurations g_config;
potfit_filenames g_files;
#if defined(KIM)
potfit_kim g_kim;
#endif // KIM
potfit_mpi_config g_mpi;
potfit_parameters g_param;
potfit_potentials g_pot;

/****************************************************************
    main potfit routine
****************************************************************/

int main(int argc, char** argv)
{
  initialize_global_variables();

  if (initialize_mpi(&argc, &argv) != POTFIT_SUCCESS) {
    shutdown_mpi();
    return EXIT_FAILURE;
  }

  read_input_files(argc, argv);

  g_mpi.init_done = 1;

#if defined(KIM)
  initialize_KIM();
#endif // KIM

  switch (broadcast_params_mpi()) {
    case POTFIT_ERROR_MPI_CLEAN_EXIT:
      shutdown_mpi();
      return EXIT_SUCCESS;
    case POTFIT_ERROR:
      shutdown_mpi();
      return EXIT_FAILURE;
  }

  g_calc.ndim = g_pot.opt_pot.idxlen;
  g_calc.ndimtot = g_pot.opt_pot.len;

  // main force vector, all forces, energies, stresses, etc. will be stored here
  g_calc.force = (double*)Malloc(g_calc.mdim * sizeof(double));

  // starting positions for the force vector
  set_force_vector_pointers();

#if defined(APOT)
#if defined(MPI)
  MPI_Bcast(g_pot.opt_pot.table, g_calc.ndimtot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif  // MPI
#if !defined(KIM)
  update_calc_table(g_pot.opt_pot.table, g_pot.calc_pot.table, 1);
#endif  // !KIM
#endif  // APOT

  if (g_mpi.myid > 0) {
    start_mpi_worker(g_calc.force);
  } else {
    time_t start_time;
    time_t end_time;

    time(&start_time);

    if (g_param.opt && g_calc.ndim > 0) {
      run_optimization();
    } else if (g_calc.ndim == 0) {
      printf(
          "\nOptimization disabled due to 0 free parameters. Calculating "
          "errors.\n");
    } else {
      printf("\nOptimization disabled. Calculating errors.\n\n");
    }

    time(&end_time);

#if defined(APOT) || defined(KIM)
    double tot = calc_forces(g_pot.opt_pot.table, g_calc.force, 0);
#else
    double tot = calc_forces(g_pot.calc_pot.table, g_calc.force, 0);
#endif  // APOT

    write_pot_table_potfit(g_files.endpot);

    {
      int format = -1;

      switch (g_pot.format_type) {
        case POTENTIAL_FORMAT_UNKNOWN:
          error(1, "Unknown potential format detected! (%s:%d)\n", __FILE__,
                __LINE__);
        case POTENTIAL_FORMAT_ANALYTIC:
          format = 0;
          break;
        case POTENTIAL_FORMAT_TABULATED_EQ_DIST:
          format = 3;
          break;
        case POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST:
          format = 4;
          break;
        case POTENTIAL_FORMAT_KIM:
          format = 5;
          break;
      }

      printf("\nPotential in format %d written to file \t%s\n", format,
             g_files.endpot);
    }

#if !defined(KIM)
    if (g_param.write_imd == 1)
      write_pot_table_imd(g_files.imdpot);

    if (g_param.plot == 1)
      write_plotpot_pair(&g_pot.calc_pot, g_files.plotfile);

    if (g_param.write_lammps == 1)
      write_pot_table_lammps();
#endif // !KIM

#if defined(BINDIST) && !defined(MPI)
    // bindist does not work with MPI
    write_bindist_file(&g_pot.opt_pot, g_files.bindistfile);
#endif // BINDIST && !MPI

    // write the error files for forces, energies, stresses, ...
    write_errors(g_calc.force, tot);

#if defined(UQ)
    for (int i=0;i<g_pot.opt_pot.idxlen;i++) {
      printf("Parameter  %d = %f ", i, g_pot.opt_pot.table[g_pot.opt_pot.idx[i]]);
    }
    printf("\n");

    ensemble_generation(tot);
#endif //UQ

    /* calculate total runtime */
    if (g_param.opt && g_mpi.myid == 0 && g_calc.ndim > 0) {
      printf("\nRuntime: %d hours, %d minutes and %d seconds.\n",
             (int)difftime(end_time, start_time) / 3600,
             ((int)difftime(end_time, start_time) % 3600) / 60,
             (int)difftime(end_time, start_time) % 60);
      printf("%d force calculations, each took %f seconds\n", g_calc.fcalls,
             (double)difftime(end_time, start_time) / g_calc.fcalls);
    }

#if defined(MPI)
    calc_forces(NULL, NULL, 1); // go wake up other threads
#endif // MPI
  } // myid == 0

  exit_potfit();
}

/****************************************************************
  read_input_files -- process all input files
****************************************************************/

void read_input_files(int argc, char** argv)
{
  // only root process reads input files

  if (g_mpi.myid == 0) {
    read_parameters(argc, argv);
#if defined(KIM)
    init_kim_model();
#endif
    read_pot_table(g_files.startpot);
    read_config(g_files.config);

    printf("Global energy weight: %f\n", g_param.eweight);
#if defined(STRESS)
    printf("Global stress weight: %f\n", g_param.sweight);
#endif  // STRESS

    /* initialize additional force variables and parameters */
    init_force_common(0);
    init_force(0);

    init_rng(g_param.rng_seed);
  }
}

/****************************************************************
  start_mpi_worker
****************************************************************/

void start_mpi_worker(double* force)
{
  /* Select correct spline interpolation and other functions */
  /* Root process has done this earlier */

  init_force_common(1);
  init_force(1);

#if defined(APOT) || defined(KIM)
  calc_forces(g_pot.opt_pot.table, force, 0);
#else
  calc_forces(g_pot.calc_pot.table, force, 0);
#endif  // APOT
}

/****************************************************************
  exit_potfit
****************************************************************/

void exit_potfit()
{
  // do some cleanups before exiting

#if defined(MPI)
  // kill MPI
  shutdown_mpi();
#endif  // MPI

#if defined(KIM)
  shutdown_KIM();
#endif  // KIM

  free_allocated_memory();

  exit(EXIT_SUCCESS);
}

/****************************************************************
  error -- complain and abort
****************************************************************/

void error(int done, const char* msg, ...)
{
  va_list ap;

  fflush(stderr);
  fprintf(stderr, "[ERROR] ");

  va_start(ap, msg);
  vfprintf(stderr, msg, ap);
  va_end(ap);

  fflush(stderr);

  if (done == 1) {
#if defined(MPI)
    if (g_mpi.init_done == 1) {
      /* go wake up other threads */
      calc_forces(NULL, NULL, 1);
      shutdown_mpi();
    }
#endif  // MPI
#if defined(KIM)
    shutdown_KIM();
#endif  // KIM
    free_allocated_memory();
    exit(EXIT_FAILURE);
  }
}

/****************************************************************
  warning -- just complain, don't abort
****************************************************************/

void warning(const char* msg, ...)
{
  va_list ap;

  fflush(stdout);
  fprintf(stderr, "[WARNING] ");

  va_start(ap, msg);
  vfprintf(stderr, msg, ap);
  va_end(ap);

  fflush(stderr);
}
