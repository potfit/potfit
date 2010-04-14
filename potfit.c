/****************************************************************
*
*  potfit.c: Contains main potfit programme.
*
*****************************************************************/
/*
*   Copyright 2002-2010 Peter Brommer, Franz G"ahler, Daniel Schopf
*             Institute for Theoretical and Applied Physics
*             University of Stuttgart, D-70550 Stuttgart, Germany
*             http://www.itap.physik.uni-stuttgart.de/
*
*****************************************************************/
/*
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
*   along with potfit; if not, write to the Free Software
*   Foundation, Inc., 51 Franklin St, Fifth Floor,
*   Boston, MA  02110-1301  USA
*/
/****************************************************************
* $Revision: 1.77 $
* $Date: 2010/04/14 10:14:17 $
*****************************************************************/

#define MAIN

#include "potfit.h"
#include "utils.h"

/******************************************************************************
*
*  error -- complain and abort
*
******************************************************************************/

void error(char *msg)
{
  fprintf(stderr, "Error: %s\n", msg);
  fflush(stderr);
#ifdef MPI
  real *force = NULL;
  calc_forces(calc_pot.table, force, 1);	/* go wake up other threads */
  shutdown_mpi();
#endif
  exit(2);
}

/******************************************************************************
 *
 *  warning -- just complain, don't abort
 *
 *****************************************************************************/

void warning(char *msg)
{
  fprintf(stderr, "Warning: %s\n", msg);
  fflush(stderr);
  return;
}

/******************************************************************************
*
*  main
*
******************************************************************************/

int main(int argc, char **argv)
{
  real *force;
  real  tot, min, max, sqr, *totdens = NULL;
  int   i, j, *ntyp = NULL;
  char  msg[255], file[255];
  FILE *outfile;

#ifdef MPI
  init_mpi(&argc, argv);
#endif

  /* assign correct force routine */
#ifdef PAIR
  calc_forces = calc_forces_pair;
  strcpy(interaction, "PAIR");
#elif defined EAM
  calc_forces = calc_forces_eam;
  strcpy(interaction, "EAM");
#elif defined MEAM
  calc_forces = calc_forces_meam;
  strcpy(interaction, "MEAM");
#elif defined ADP
  calc_forces = calc_forces_adp;
  strcpy(interaction, "ADP");
#endif

  /* read the parameters and the potential file */
  if (myid == 0) {
    read_parameters(argc, argv);
    read_pot_table(&opt_pot, startpot, ntypes * (ntypes + 1) / 2);
    read_config(config);
    printf("Global energy weight: %f\n", eweight);
#ifdef STRESS
    printf("Global stress weight: %f\n", sweight);
#endif
    /* Select correct spline interpolation and other functions */
    if (format == 0) {
#ifndef APOT
      sprintf(msg,
	      "potfit binary compiled without analytic potential support\n");
      error(msg);
#else
      splint = splint_ed;
      splint_comb = splint_comb_ed;
      splint_grad = splint_grad_ed;
      write_pot_table = write_apot_table;
#endif
    } else if (format == 3) {
#ifdef APOT
      sprintf(msg,
	      "potfit binary compiled without tabulated potential support\n");
      error(msg);
#else
      splint = splint_ed;
      splint_comb = splint_comb_ed;
      splint_grad = splint_grad_ed;
      write_pot_table = write_pot_table3;
#ifdef PARABEL
      parab = parab_ed;
      parab_comb = parab_comb_ed;
      parab_grad = parab_grad_ed;
#endif /* PARABEL */
#endif /* APOT */
    } else if (format >= 4) {	/*format == 4 ! */
#ifdef APOT
      sprintf(msg,
	      "potfit binary compiled without tabulated potential support\n");
      error(msg);
#else
      splint = splint_ne;
      splint_comb = splint_comb_ne;
      splint_grad = splint_grad_ne;
      write_pot_table = write_pot_table4;
#ifdef PARABEL
      parab = parab_ne;
      parab_comb = parab_comb_ne;
      parab_grad = parab_grad_ne;
#endif /* PARABEL */
#endif /* APOT */
    }

    /* set spline density corrections to 0 */
#if defined EAM || defined MEAM || defined ADP
    lambda = (real *)malloc(ntypes * sizeof(real));
    reg_for_free(lambda, "lambda");

    totdens = (real *)malloc(ntypes * sizeof(real));
    reg_for_free(totdens, "totdens");

    ntyp = (int *)malloc(ntypes * sizeof(int));
    reg_for_free(ntyp, "ntyp");

    for (i = 0; i < ntypes; i++)
      lambda[i] = 0.;
#ifndef NORESCALE
    rescale(&opt_pot, 1., 1);	/* rescale now... */
#ifdef WZERO
//    embed_shift(&opt_pot);    /* and shift  - diabolical */
#endif /* WZERO */
#endif /* NORESCALE */
#endif /* EAM */
    init_done = 1;
  }

  /* initialize random number generator */
  srandom(seed + myid);

  /* initialize the remaining parameters and assign the atoms */
#ifdef MPI
  MPI_Bcast(&init_done, 1, MPI_INT, 0, MPI_COMM_WORLD);
  broadcast_params();		/* let the others know what's going on */
#else /* MPI */
  /* Identify subset of atoms/volumes belonging to individual process
     with complete set of atoms/volumes */
  conf_atoms = atoms;
  conf_vol = volumen;
  conf_uf = useforce;
  conf_us = usestress;
  myatoms = natoms;
#endif /* MPI */
  /*   mdim=3*natoms+nconf; */
  ndim = opt_pot.idxlen;
  ndimtot = opt_pot.len;
  paircol = (ntypes * (ntypes + 1)) / 2;
  idx = opt_pot.idx;

  /* main force vector, all forces, energies, ... will be stored here */
  force = (real *)malloc((mdim) * sizeof(real));
  reg_for_free(force, "force");

  /* starting positions for the force vector */
  energy_p = 3 * natoms;
  stress_p = energy_p + nconf;
#if defined EAM || defined MEAM || defined ADP
  limit_p = stress_p + 6 * nconf;
  dummy_p = limit_p + nconf;
#ifdef APOT
  punish_par_p = dummy_p + 2 * ntypes;
  punish_pot_p = punish_par_p + apot_table.total_par;
#endif
#else /* EAM */
#ifdef APOT
  punish_par_p = stress_p + 6 * nconf;
  punish_pot_p = punish_par_p + apot_table.total_par;
#endif
#endif /* EAM MEAM */
  rms = (real *)malloc(3 * sizeof(real));
  reg_for_free(rms, "rms");

#ifdef APOT
#ifdef MPI
  MPI_Bcast(opt_pot.table, ndimtot, REAL, 0, MPI_COMM_WORLD);
#endif
  update_calc_table(opt_pot.table, calc_pot.table, 1);
  for (i = 0; i < paircol; i++)
    new_slots(i, 1);
#endif

  /* Select correct spline interpolation and other functions */
  /* Root process has done this earlier */
  if (myid > 0) {
    if (format == 0) {
#ifdef APOT
      splint = splint_ed;
      splint_comb = splint_comb_ed;
      splint_grad = splint_grad_ed;
      write_pot_table = write_apot_table;
#endif
    } else if (format == 3) {
#ifndef APOT
      splint = splint_ed;
      splint_comb = splint_comb_ed;
      splint_grad = splint_grad_ed;
      write_pot_table = write_pot_table3;
#ifdef PARABEL
      parab = parab_ed;
      parab_comb = parab_comb_ed;
      parab_grad = parab_grad_ed;
#endif /* PARABEL */
#endif /* APOT */
    } else if (format == 4) {	/*format == 4 ! */
#ifndef APOT
      splint = splint_ne;
      splint_comb = splint_comb_ne;
      splint_grad = splint_grad_ne;
      write_pot_table = write_pot_table4;
#ifdef PARABEL
      parab = parab_ne;
      parab_comb = parab_comb_ne;
      parab_grad = parab_grad_ne;
#endif /* PARABEL */
#endif /* APOT */
    }

    /* all but root go to calc_forces */
#ifndef APOT
    calc_forces(calc_pot.table, force, 0);
#else
    calc_forces(opt_pot.table, force, 0);
#endif
  } else {			/* root thread does minimization */
    if (opt) {
      printf("\nStarting optimization with %d parameters ...\n", ndim);
      fflush(stdout);
#ifdef EVO
      diff_evo(opt_pot.table);
#else
      anneal(opt_pot.table);
#endif
      powell_lsq(opt_pot.table);
      printf("\nFinished powell minimization, calculating errors ...\n");
    } else {
      printf("\nOptimization disabled. Calculating errors.\n\n");
    }

#ifndef APOT
    tot = calc_forces(calc_pot.table, force, 0);
    if (opt) {
      write_pot_table(&opt_pot, endpot);
#else
    tot = calc_forces(opt_pot.table, force, 0);
    if (opt) {
      write_pot_table(&apot_table, endpot);
#endif
      printf("\nPotential in format %d written to file %s\n", format, endpot);
    }
/*#ifndef APOT*/
/*    printf("Plotpoint file written to file %s\n", plotpointfile);*/
/*#endif*/
    if (writeimd)
      write_pot_table_imd(&calc_pot, imdpot);
    if (plot)
      write_plotpot_pair(&calc_pot, plotfile);

#ifdef PDIST
#ifndef MPI			/* will not work with MPI */
    write_pairdist(&opt_pot, distfile);
#endif
#endif
    if (format == 3) {		/* then we can also write format 4 */
      sprintf(endpot, "%s_4", endpot);
      write_pot_table4(&opt_pot, endpot);
      printf("Potential in format 4 written to file %s\n", endpot);
    }
#if defined EAM || defined MEAM
#ifndef MPI
/* Not much sense in printing rho when not communicated... */
    if (write_output_files) {
      strcpy(file, output_prefix);
      strcat(file, ".rho_loc");
      outfile = fopen(file, "w");
      if (NULL == outfile) {
	sprintf(msg, "Could not open file %s\n", file);
	error(msg);
      }
    } else {
      outfile = stdout;
      printf("Local electron density rho\n");
    }
    for (i = 0; i < ntypes; i++) {
      totdens[i] = 0.;
      ntyp[i] = 0;
    }
    fprintf(outfile, "#    atomtype\trho\n");
    for (i = 0; i < natoms; i++) {
      fprintf(outfile, "%d\t%d\t%f\n", i, atoms[i].typ, atoms[i].rho);
      totdens[atoms[i].typ] += atoms[i].rho;
      ntyp[atoms[i].typ]++;
    }
    if (write_output_files) {
      printf("Local electron density data written to %s\n", file);
      fclose(outfile);
    }
    for (i = 0; i < ntypes; i++) {
      totdens[i] /= (real)ntyp[i];
      printf("Average local electron density at atom sites type %d: %f\n",
	     i, totdens[i]);
    }
#ifdef NEWSCALE
    for (i = 0; i < ntypes; i++) {
#ifdef NORESCALE
      /* U'(1.) = 0. */
      lambda[i] = splint_grad(&calc_pot, calc_pot.table,
			      paircol + ntypes + i, 1.0);
#else /* NORESCALE */
      /* U'(<n>)=0; */
      lambda[i] = splint_grad(&calc_pot, calc_pot.table,
			      paircol + ntypes + i, totdens[i]);
#endif /* NORESCALE */
      printf("lambda[%d] = %f \n", i, lambda[i]);
    }
    sprintf(plotfile, "%s_new", plotfile);
    sprintf(imdpot, "%s_new", imdpot);
    /* write new potential plotting table */
    if (plot)
      write_altplot_pair(&opt_pot, plotfile);
    /* write NEW imd potentials */
    if (writeimd)
      write_pot_table_imd(&opt_pot, imdpot);
#endif /* NEWSCALE */
#endif /* MPI */
#endif /* EAM MEAM */

    /* prepare for error calculations */
    max = 0.0;
    min = 100000.0;
    real  f_sum = 0., e_sum = 0., s_sum = 0.;

    /* write force deviations */
    if (write_output_files) {
      strcpy(file, output_prefix);
      strcat(file, ".force");
      outfile = fopen(file, "w");
      if (NULL == outfile) {
	sprintf(msg, "Could not open file %s\n", file);
	error(msg);
      }
    } else {
      outfile = stdout;
      printf("Forces:\n");
    }
    for (i = 0; i < 3 * natoms; i++) {
      sqr = conf_weight[atoms[i / 3].conf] * SQR(force[i]);
      f_sum += sqr;
      max = MAX(max, sqr);
      min = MIN(min, sqr);
#ifdef FWEIGHT
      if (i == 0)
	fprintf(outfile,
		"conf:atom\ttype\t(w*df)^2\t\tf\t\tf0\t\tdf/f0\t\t|f|\n");
      fprintf(outfile,
	      "%3d:%5d\t%4s\t%11.8f\t%11.8f\t%11.8f\t%11.8f\t%11.8f\n",
	      atoms[i / 3].conf, i / 3, elements[atoms[i / 3].typ], sqr,
	      force[i] * (FORCE_EPS + atoms[i / 3].absforce) + force_0[i],
	      force_0[i],
	      (force[i] * (FORCE_EPS + atoms[i / 3].absforce)) / force_0[i],
	      atoms[i / 3].absforce);
#else /* FWEIGHT */
      if (i == 0)
	fprintf(outfile, "conf:atom\ttype\t(w*df)^2\t\tf\t\tf0\t\tdf/f0\n");
      fprintf(outfile, "%3d:%5d\t%4s\t%.8f\t%11.8f\t%11.8f\t%11.8f\n",
	      atoms[i / 3].conf, i / 3, elements[atoms[i / 3].typ], sqr,
	      force[i] + force_0[i], force_0[i], force[i] / force_0[i]);
#endif /* FWEIGHT */
    }
    if (write_output_files) {
      printf("Force data written to %s\n", file);
      fclose(outfile);
    }

    /* write energy deviations */
    if (write_output_files) {
      strcpy(file, output_prefix);
      strcat(file, ".energy");
      outfile = fopen(file, "w");
      if (NULL == outfile) {
	sprintf(msg, "Could not open file %s\n", file);
	error(msg);
      }
    } else {
      outfile = stdout;
      printf("Cohesive Energies\n");
    }

    if (write_output_files) {
      fprintf(outfile, "# global energy weight w is %f\n", eweight);
      fprintf(outfile,
	      "# nr.\tconf_w\t(w*de)^2\te\t\te0\t\t|e-e0|\t\te-e0\t\tde/e0\n");
    } else {
      fprintf(outfile, "energy weight is %f\n", eweight);
      fprintf(outfile, "conf\tconf_w\t(w*de)^2\te\t\te0\t\tde/e0\n");
    }

    for (i = 0; i < nconf; i++) {
      sqr = conf_weight[i] * SQR(force[energy_p + i]);
      e_sum += sqr;
      max = MAX(max, sqr);
      min = MIN(min, sqr);
      if (write_output_files) {
	fprintf(outfile, "%3d\t%.4f\t%f\t%.10f\t%.10f\t%f\t%f\t%f\n", i,
		conf_weight[i], sqr,
		(force[energy_p + i] + force_0[energy_p + i]) / eweight,
		force_0[energy_p + i] / eweight,
		fabs(force[energy_p + i]) / eweight,
		force[energy_p + i] / eweight,
		force[energy_p + i] / force_0[energy_p + i]);
      } else
	fprintf(outfile, "%d\t%.4f\t%f\t%f\t%f\t%f\n", i, conf_weight[i], sqr,
		(force[energy_p + i] + force_0[energy_p + i]) / eweight,
		force_0[energy_p + i] / eweight,
		force[energy_p + i] / force_0[energy_p + i]);
    }
    if (write_output_files) {
      printf("Energy data written to %s\n", file);
      fclose(outfile);
    }
#ifdef STRESS
    /* write stress deviations */
    if (write_output_files) {
      strcpy(file, output_prefix);
      strcat(file, ".stress");
      outfile = fopen(file, "w");
      if (NULL == outfile) {
	sprintf(msg, "Could not open file %s\n", file);
	error(msg);
      }
      fprintf(outfile, "# global stress weight w is %f\n", sweight);
    } else {
      outfile = stdout;
      fprintf(outfile, "Stresses on unit cell\n");
    }
    fprintf(outfile, "#\tconf_w\t\t(w*ds)^2\ts\t\ts0\t\tds/s0\n");
    for (i = stress_p; i < stress_p + 6 * nconf; i++) {
      sqr = SQR(force[i]);
      s_sum += sqr;
      max = MAX(max, sqr);
      min = MIN(min, sqr);
      fprintf(outfile, "%d\t%f\t%f\t%f\t%f\t%f\n",
	      (i - stress_p) / 6, conf_weight[(i - stress_p) / 6], sqr,
	      (force[i] + force_0[i]) / sweight, force_0[i] / sweight,
	      force[i] / force_0[i]);
    }
    if (write_output_files) {
      printf("Stress data written to %s\n", file);
      fclose(outfile);
    }
#endif
#if defined EAM || defined MEAM
/*    if (opt) {*/
    /* write EAM punishments */
    if (write_output_files) {
      strcpy(file, output_prefix);
      strcat(file, ".punish");
      outfile = fopen(file, "w");
      if (NULL == outfile) {
	sprintf(msg, "Could not open file %s\n", file);
	error(msg);
      }
      fprintf(outfile, "Limiting constraints\n");
      fprintf(outfile, "#conf\tp^2\t\tpunishment\n");
    } else {
      outfile = stdout;
      printf("Punishment Constraints\n");
    }
    for (i = limit_p; i < dummy_p; i++) {
      sqr = SQR(force[i]);
      max = MAX(max, sqr);
      min = MIN(min, sqr);
      if (write_output_files)
	fprintf(outfile, "%d\t%f\t%f\n", i - limit_p, sqr,
		force[i] + force_0[i]);
      else
	fprintf(outfile, "%d %f %f %f %f\n", i - limit_p,
		sqr, force[i] + force_0[i], force_0[i],
		force[i] / force_0[i]);
    }
    if (write_output_files) {
      real  zero = 0;
      fprintf(outfile, "\nDummy Constraints\n");
      fprintf(outfile, "element\tU^2\t\tU'^2\t\tU\t\tU'\n");
      for (i = dummy_p; i < dummy_p + ntypes; i++) {
#ifdef NORESCALE
	sqr = SQR(force[i]);
	max = MAX(max, sqr);
	min = MIN(min, sqr);
	fprintf(outfile, "%s\t%f\t%f\t%f\t%g\n", elements[i - dummy_p],
		zero, sqr, zero, force[i]);
#else
	sqr = SQR(force[i]);
	max = MAX(max, sqr);
	min = MIN(min, sqr);
	sqr = SQR(force[i + ntypes]);
	max = MAX(max, sqr);
	min = MIN(min, sqr);
	fprintf(outfile, "%s\t%f\t%f\t%f\t%f\n", elements[i - dummy_p], sqr,
		SQR(force[i]), force[i + ntypes], force[i]);
#endif
      }
#ifdef NORESCALE
      fprintf(outfile, "\nNORESCALE: <n>!=1\n");
      fprintf(outfile, "<n>=%f\n",
	      force[dummy_p + ntypes] / DUMMY_WEIGHT + 1);
      fprintf(outfile, "Additional punishment of %f added.\n",
	      SQR(force[dummy_p + ntypes]));
#endif
      printf("Punishment constraints data written to %s\n", file);
      fclose(outfile);
    } else {
      fprintf(outfile, "Dummy Constraints\n");
      for (i = 0; i < 2 * ntypes; i++) {
	sqr = SQR(force[i]);
	max = MAX(max, sqr);
	min = MIN(min, sqr);
	fprintf(outfile, "%d %f %f %f %f\n", i - dummy_p, sqr,
		force[dummy_p + i] + force_0[dummy_p + i],
		force_0[dummy_p + i],
		force[dummy_p + i] / force_0[dummy_p + i]);
      }
    }
/*    }*/
#endif

    /* final error report */
    printf("\n###### error report ######\n");
    if (write_output_files) {
      strcpy(file, output_prefix);
      strcat(file, ".error");
      outfile = fopen(file, "w");
      if (NULL == outfile) {
	sprintf(msg, "Could not open file %s\n", file);
	error(msg);
      }
#ifndef STRESS
      fprintf(outfile,
	      "total error sum %f, count %d (%d forces, %d energies)\n", tot,
	      mdim - 6 * nconf, 3 * natoms, nconf);
#else
      fprintf
	(outfile,
	 "total error sum %f, count %d (%d forces, %d energies, %d stresses)\n",
	 tot, mdim, 3 * natoms, nconf, 6 * nconf);
#endif
    }
#ifndef STRESS
    printf("total error sum %f, count %d (%d forces, %d energies)\n", tot,
	   mdim - 6 * nconf, 3 * natoms, nconf);
#else
    printf
      ("total error sum %f, count %d (%d forces, %d energies, %d stresses)\n",
       tot, mdim, 3 * natoms, nconf, 6 * nconf);
#endif

    /* calculate the rms errors for forces, energies, stress */
    rms[0] = 0;			/* rms rms for forces */
    rms[1] = 0;			/* energies */
    rms[2] = 0;			/* stresses */

    /* forces */
    for (i = 0; i < 3 * natoms; i++)
      rms[0] += SQR(force[i]);
    rms[0] = sqrt(rms[0] / natoms);

    /* energies */
    for (i = 0; i < nconf; i++)
      rms[1] += SQR(force[3 * natoms + i] / eweight);
    if (isnan(rms[1]))
      rms[1] = 0;
    rms[1] = sqrt(rms[1] / nconf);

    /* stresses */
    for (i = 0; i < nconf; i++)
      for (j = 0; j < 6; j++)
	rms[2] += SQR(force[3 * natoms + nconf + 6 * i + j] / sweight);
    if (isnan(rms[2]))
      rms[2] = 0;
    rms[2] = sqrt(rms[2] / (6 * nconf));

    if (write_output_files) {
      fprintf(outfile, "sum of force-errors = %f\t\t( %.3f%% - av: %f)\n",
	      f_sum, f_sum / tot * 100, f_sum / (3 * natoms));
      fprintf(outfile, "sum of energy-errors = %f\t\t( %.3f%% )\n", e_sum,
	      e_sum / tot * 100);
#ifdef STRESS
      fprintf(outfile, "sum of stress-errors = %f\t\t( %.3f%% )\n", s_sum,
	      s_sum / tot * 100);
#endif
      if ((tot - f_sum - e_sum - s_sum) > 0.01 && opt == 1) {
	fprintf
	  (outfile,
	   "\n --> Warning <--\nThis sum contains punishments! Check your results.\n");
	fprintf(outfile, "sum of punishments = %f\t\t( %.3f%% )\n\n",
		tot - f_sum - e_sum - s_sum,
		(tot - f_sum - e_sum - s_sum) / tot * 100);
      }
      fprintf(outfile, "min: %e - max: %e\n", min, max);
      fprintf(outfile, "rms-errors:\n");
      fprintf(outfile, "force \t%f\t(%f meV/A)\n", rms[0], rms[0] * 1000);
      fprintf(outfile, "energy \t%f\t(%f meV)\n", rms[1], rms[1] * 1000);
#ifdef STRESS
      fprintf(outfile, "stress \t%f\t(%f MPa)\n", rms[2],
	      rms[2] / 160.2 * 1000);
#endif
      fprintf(outfile, "\n");
      fprintf(outfile,
	      "\tforce [meV/A]\tenergy [meV]\tstress [MPa]\terror sum\n");
      fprintf(outfile, "RMS:\t%f\t%f\t%f\t%f\n", rms[0] * 1000, rms[1] * 1000,
	      rms[2] / 160.2 * 1000, tot);

    }
    printf("sum of force-errors = %f\t\t( %.3f%% - av: %f)\n",
	   f_sum, f_sum / tot * 100, f_sum / (3 * natoms));
    printf("sum of energy-errors = %f\t\t( %.3f%% )\n", e_sum,
	   e_sum / tot * 100);
#ifdef STRESS
    printf("sum of stress-errors = %f\t\t( %.3f%% )\n", s_sum,
	   s_sum / tot * 100);
#endif
    if ((tot - f_sum - e_sum - s_sum) > 0.01 && opt == 1) {
      printf
	("\n --> Warning <--\nThis sum contains punishments! Check your results.\n");
      printf("sum of punishments = %f\t\t( %.3f%% )\n\n",
	     tot - f_sum - e_sum - s_sum,
	     (tot - f_sum - e_sum - s_sum) / tot * 100);
    }
    printf("min: %e - max: %e\n", min, max);
    printf("rms-errors:\n");
    printf("force \t%f\t(%f meV/A)\n", rms[0], rms[0] * 1000);
    printf("energy \t%f\t(%f meV)\n", rms[1], rms[1] * 1000);
#ifdef STRESS
    printf("stress \t%f\t(%f MPa)\n", rms[2], rms[2] / 160.2 * 1000);
#endif
#ifdef MPI
    calc_forces(calc_pot.table, force, 1);	/* go wake up other threads */
#endif /* MPI */
  }
  /* do some cleanups before exiting */
#ifdef MPI
  /* kill MPI */
  shutdown_mpi();
#else
  free_all_pointers();
#endif
  return 0;
}
