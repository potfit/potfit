/****************************************************************
 *
 * potfit.c: Contains main potfit program
 *
 ****************************************************************
 *
 * Copyright 2002-2014
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.sourceforge.net/
 *
 ****************************************************************
 *
 *   This file is part of potfit.
 *
 *   potfit is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   potfit is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

#define MAIN
#include "potfit.h"
#undef MAIN

#include <time.h>

#include "config.h"
#include "functions.h"
#include "optimize.h"
#include "potential.h"
#include "splines.h"
#include "utils.h"

/****************************************************************
 *
 *  main potfit routine
 *
 ****************************************************************/

int main(int argc, char **argv)
{
  int   i;
  double *force, tot;
  time_t t_begin, t_end;

#ifdef MPI
  /* initialize the MPI communication */
  init_mpi(argc, argv);
#endif /* MPI */

  /* print version and compile date and time */
  if (myid == 0) {
    printf("This is %s compiled on %s, %s.\n\n", POTFIT_VERSION, __DATE__, __TIME__);
#ifdef MPI
    printf("Starting up MPI with %d processes.\n", num_cpus);
#endif /* MPI */
  }

  /* read the parameters and the potential file */
  if (myid == 0) {
    read_parameters(argc, argv);
    read_pot_table(&opt_pot, startpot);
    read_config(config);
    printf("Global energy weight: %f\n", eweight);
#ifdef STRESS
    printf("Global stress weight: %f\n", sweight);
#endif /* STRESS */

    /* Select correct spline interpolation and other functions */
#ifdef APOT
    if (format == 0) {
      splint = splint_ed;
      splint_comb = splint_comb_ed;
      splint_grad = splint_grad_ed;
      write_pot_table = write_pot_table0;
    }
#else /* APOT */
    if (format == 3) {
      splint = splint_ed;
      splint_comb = splint_comb_ed;
      splint_grad = splint_grad_ed;
      write_pot_table = write_pot_table3;
    } else if (format >= 4) {
      splint = splint_ne;
      splint_comb = splint_comb_ne;
      splint_grad = splint_grad_ne;
      write_pot_table = write_pot_table4;
    }
#endif /* APOT */

    /* initialize additional force variables and parameters */
    init_forces();

    init_done = 1;

    /* properly initialize random number generator */
#define R_SIZE 624
#define RAND_MAX 2147483647
    {
      uint32_t *array;
      array = (uint32_t *) malloc(R_SIZE * sizeof(uint32_t));
      srand(seed);
      for (i = 0; i < R_SIZE; i++)
	array[i] = rand();

      dsfmt_init_by_array(&dsfmt, array, R_SIZE);
      for (i = 0; i < 10e5; i++)
	eqdist();
      free(array);
    }
#undef R_SIZE
#undef RAND_MAX
  }
  /* myid == 0 */
#ifdef MPI
  if (0 == myid) {
    printf("Broadcasting data to MPI workers ... ");
    fflush(stdout);
  }
  broadcast_params();		/* let the others know what's going on */
  if (0 == myid) {
    printf("done\n");
    fflush(stdout);
  }
#else
  /* Identify subset of atoms/volumes belonging to individual process
     with complete set of atoms/volumes */
  conf_atoms = atoms;
  conf_vol = volume;
  conf_uf = useforce;
#ifdef STRESS
  conf_us = usestress;
#endif /* STRESS */
#endif /* MPI */

  ndim = opt_pot.idxlen;
  ndimtot = opt_pot.len;
  idx = opt_pot.idx;

  /* main force vector, all forces, energies, stresses, ... will be stored here */
  force = (double *)malloc((mdim) * sizeof(double));
  if (NULL == force)
    error(1, "Could not allocate memory for main force vector.");
  for (i = 0; i < mdim; i++)
    force[i] = 0.0;
  reg_for_free(force, "force vector");

  /* starting positions for the force vector */
  set_force_vector_pointers();

#ifdef APOT
#ifdef MPI
  MPI_Bcast(opt_pot.table, ndimtot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif /* MPI */
  update_calc_table(opt_pot.table, calc_pot.table, 1);
#endif /* APOT */

  /* Select correct spline interpolation and other functions */
  /* Root process has done this earlier */
  if (myid > 0) {
#ifdef APOT
    if (format == 0) {
      splint = splint_ed;
      splint_comb = splint_comb_ed;
      splint_grad = splint_grad_ed;
    }
#else /* APOT */
    if (format == 3) {
      splint = splint_ed;
      splint_comb = splint_comb_ed;
      splint_grad = splint_grad_ed;
    } else if (format >= 4) {
      splint = splint_ne;
      splint_comb = splint_comb_ne;
      splint_grad = splint_grad_ne;
    }
#endif /* !APOT */

    /* all but root go to calc_forces */
#ifdef APOT
    calc_forces(opt_pot.table, force, 0);
#else
    calc_forces(calc_pot.table, force, 0);
#endif /* APOT */
  } else {			/* root thread does minimization */
#ifdef MPI
    if (num_cpus > nconf) {
      warning("You are using more CPUs than you have configurations!\n");
      warning("While this will not do any harm, you are wasting %d CPUs.\n", num_cpus - nconf);
    }
#endif /* MPI */
    time(&t_begin);
    if (opt && ndim > 0) {
      printf("\nStarting optimization with %d parameters.\n", ndim);
      fflush(stdout);
#ifndef EVO
      anneal(opt_pot.table);
#else /* EVO */
      diff_evo(opt_pot.table);
#endif /* EVO */
      printf("\nStarting powell minimization ...\n");
      powell_lsq(opt_pot.table);
      printf("\nFinished powell minimization, calculating errors ...\n");
    } else if (0 == ndim) {
      printf("\nOptimization disabled due to 0 free parameters. Calculating errors.\n");
    } else {
      printf("\nOptimization disabled. Calculating errors.\n\n");
    }
    time(&t_end);

#ifdef APOT
    tot = calc_forces(opt_pot.table, force, 0);
    if (opt) {
      write_pot_table(&apot_table, endpot);
#else /* APOT */
    tot = calc_forces(calc_pot.table, force, 0);
    if (opt) {
      write_pot_table(&opt_pot, endpot);
#endif /* APOT */
      printf("\nPotential in format %d written to file \t%s\n", format, endpot);
    }
#ifndef APOT
    /* then we can also write format 4 */
    if (format == 3) {
      sprintf(endpot, "%s_4", endpot);
      write_pot_table4(&opt_pot, endpot);
      printf("Potential in format 4 written to file \t%s\n", endpot);
    }
#endif /* !APOT */

    if (1 == writeimd)
      write_pot_table_imd(&calc_pot, imdpot);

    if (1 == plot)
      write_plotpot_pair(&calc_pot, plotfile);

#ifdef APOT
    if (1 == write_lammps)
      write_pot_table_lammps(&calc_pot);
#endif /* APOT */

    /* will not work with MPI */
#if defined PDIST && !defined MPI
    write_pairdist(&opt_pot, distfile);
#endif /* PDIST && !MPI */

    /* write the error files for forces, energies, stresses, ... */
    write_errors(force, tot);

#ifdef MPI
    calc_forces(calc_pot.table, force, 1);	/* go wake up other threads */
#endif /* MPI */
  }				/* myid == 0 */

  /* calculate total runtime */
  if (opt && myid == 0 && ndim > 0) {
    printf("\nRuntime: %d hours, %d minutes and %d seconds.\n",
      (int)difftime(t_end, t_begin) / 3600, ((int)difftime(t_end,
	  t_begin) % 3600) / 60, (int)difftime(t_end, t_begin) % 60);
    printf("%d force calculations, each took %f seconds\n", fcalls, (double)difftime(t_end,
	t_begin) / fcalls);
  }

  /* do some cleanups before exiting */
#ifdef MPI
  /* kill MPI */
  shutdown_mpi();
#endif /* MPI */

  free(u_address);
  free_all_pointers();


/* added */
/* free KIM stuff */

/*	FreeKIM();
*/
/* added ends */


  return 0;
}

/****************************************************************
 *
 *  error -- complain and abort
 *
 ****************************************************************/

void error(int done, char *msg, ...)
{
  va_list ap;

  fflush(stderr);
  fprintf(stderr, "[ERROR] ");

  va_start(ap, msg);
  vfprintf(stderr, msg, ap);
  va_end(ap);

  fflush(stderr);
  if (done == 1) {
#ifdef MPI
    double *force = NULL;
    /* go wake up other threads */
    calc_forces(calc_pot.table, force, 1);
    fprintf(stderr, "\n");
    shutdown_mpi();
#endif /* MPI */
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
  }
}

/****************************************************************
 *
 *  warning -- just complain, don't abort
 *
 ****************************************************************/

void warning(char *msg, ...)
{
  va_list ap;

  fflush(stdout);
  fprintf(stderr, "[WARNING] ");

  va_start(ap, msg);
  vfprintf(stderr, msg, ap);
  va_end(ap);

  fflush(stderr);
}
