/****************************************************************
* 
*  potfit.c: Contains main potfit programme.
*
*****************************************************************/
/*
*   Copyright 2002-2009 Peter Brommer, Franz G"ahler, Daniel Schopf
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
* $Revision: 1.53 $
* $Date: 2009/04/21 13:48:09 $
*****************************************************************/

#define MAIN

#include "potfit.h"

/******************************************************************************
*
*  error -- complain and abort
*
******************************************************************************/

void error(char *msg)
{
  real *force = NULL;
  fprintf(stderr, "Error: %s\n", msg);
  fflush(stderr);
#ifdef MPI
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
  real  tot, min, max, sqr, *totdens;
  int   i, j, diff, *ntyp;
  char  msg[255];
  FILE *outforce, *outenergy, *outstress;

  pi = 4.0 * atan(1.);
#ifdef MPI
  init_mpi(&argc, argv);
#endif
  srandom(seed + myid);
  random();
  random();
  random();
  random();
  calc_forces = calc_forces_pair;
  if (myid == 0) {
    read_parameters(argc, argv);
#ifdef APOT
    read_pot_table(&opt_pot, &apot_table, startpot,
		   ntypes * (ntypes + 1) / 2);
#else
    read_pot_table(&opt_pot, startpot, ntypes * (ntypes + 1) / 2);
#endif
    read_config(config);
    printf("Energy weight: %f\n", eweight);
#ifdef STRESS
    printf("Stress weight: %f\n", sweight);
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
#ifdef EAM
    lambda = (real *)malloc(ntypes * sizeof(real));
    totdens = (real *)malloc(ntypes * sizeof(real));
    ntyp = (int *)malloc(ntypes * sizeof(int));
    for (i = 0; i < ntypes; i++)
      lambda[i] = 0.;
#ifndef NORESCALE
    rescale(&opt_pot, 1., 1);	/* rescale now... */
#ifdef WZERO
//    embed_shift(&opt_pot);    /* and shift  - diabolical */
#endif /* WZERO */
#endif /* NORESCALE */
#endif /* EAM */
  }
#ifdef MPI
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

  force = (real *)malloc((mdim) * sizeof(real));
  rms = (real *)malloc(3 * sizeof(real));

#ifdef APOT
#ifdef MPI
  MPI_Bcast(opt_pot.table, ndimtot, REAL, 0, MPI_COMM_WORLD);
#endif
  update_calc_table(opt_pot.table, calc_pot.table, 1);
#endif

  if (myid > 0) {
    /* Select correct spline interpolation and other functions */
    /* Root process has done this earlier */
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
      anneal(opt_pot.table);
      powell_lsq(opt_pot.table);
    }
/*  for (i=0; i<opt_pot.ncols; i++) 
      spline_ed(opt_pot.step[i],opt_pot.table+opt_pot.first[i],
      opt_pot.last[i]-opt_pot.first[i]+1,
      1e30,0,opt_pot.d2tab+opt_pot.first[i]);*/

//    rescale(&opt_pot,1.);
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
#ifndef APOT
    printf("Plotpoint file written to file %s\n", plotpointfile);
#endif
    if (writeimd)
      write_pot_table_imd(&calc_pot, imdpot);
    if (plot)
      write_plotpot_pair(&calc_pot, plotfile);
//    if (plot) write_altplot_pair(&opt_pot, plotfile);


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
#ifdef EAM
#ifndef MPI
/* Not much sense in printing rho when not communicated... */
    printf("Local electron density rho\n");
    for (i = 0; i < ntypes; i++) {
      totdens[i] = 0.;
      ntyp[i] = 0;
    }
    for (i = 0; i < natoms; i++) {
      printf("%d %d %f\n", i, atoms[i].typ, atoms[i].rho);
      totdens[atoms[i].typ] += atoms[i].rho;
      ntyp[atoms[i].typ]++;
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
#endif /* EAM */

    max = 0.0;
    min = 100000.0;
    if (strcmp(endforce, "stdout") != 0) {
      outforce = fopen(endforce, "w");
      if (NULL == outforce) {
	sprintf(msg, "Could not open file %s\n", endforce);
	error(msg);
      }
    } else
      outforce = stdout;
    fprintf(outforce, "conf-atom   df^2\t f\t\t f0\t\t df/f0\t\t |f|\n");
    for (i = 0; i < 3 * natoms; i++) {
      sqr = SQR(force[i]);
      max = MAX(max, sqr);
      min = MIN(min, sqr);
#ifdef FWEIGHT
      fprintf(outforce, "%d-%d\t%f\t%f\t%f\t%f\t%f\n", atoms[i / 3].conf,
	      i / 3, sqr,
	      force[i] * (FORCE_EPS + atoms[i / 3].absforce) + force_0[i],
	      force_0[i],
	      (force[i] * (FORCE_EPS + atoms[i / 3].absforce)) / force_0[i],
	      atoms[i / 3].absforce);
#else /* FWEIGHT */
      fprintf(outforce, "%d-%d\t%f\t%f\t%f\t%f\n", atoms[i / 3].conf, i / 3,
	      sqr, force[i] + force_0[i], force_0[i], force[i] / force_0[i]);
#endif /* FWEIGHT */
    }
    if (strcmp(endforce, "stdout") != 0) {
      printf("Force data written to %s\n", endforce);
      fclose(outforce);
    }

    if (strcmp(endenergy, "stdout") != 0) {
      outenergy = fopen(endenergy, "w");
      if (NULL == outenergy) {
	sprintf(msg, "Could not open file %s\n", endenergy);
	error(msg);
      }
    } else {
      outenergy = stdout;
      printf("Cohesive Energies\n");
    }

    if (strcmp(endenergy, "stdout") != 0) {
      strcpy(msg, "");
      for (j = 4; j < config_name_max; j++)
	strcat(msg, " ");
      fprintf(outenergy,
	      "#\tconf%s\t(w*de)^2\te\t\te0\t\t|e-e0|\t\te-e0\t\tde/e0\n",
	      msg);
    } else
      fprintf(outenergy, "#\tconf\t(w*de)^2\te\t\te0\t\tde/e0\n");

    real  e_sum = 0, s_sum = 0;

    for (i = 0; i < nconf; i++) {
      sqr = SQR(force[3 * natoms + i]);
      e_sum += sqr;
      max = MAX(max, sqr);
      min = MIN(min, sqr);
      strcpy(msg, "");
      if (strcmp(endenergy, "stdout") != 0) {
	if (strlen(config_name[i]) < config_name_max)
	  for (j = 0; j < (config_name_max - strlen(config_name[i])); j++)
	    strcat(msg, " ");
	fprintf(outenergy, "%d\t%s%s\t%f\t%f\t%f\t%f\t%f\t%f\n", i,
		config_name[i], msg, sqr,
		(force[3 * natoms + i] + force_0[3 * natoms + i]) / eweight,
		force_0[3 * natoms + i] / eweight,
		fabs(force[3 * natoms + i]) / eweight,
		force[3 * natoms + i] / eweight,
		force[3 * natoms + i] / force_0[3 * natoms + i]);
      } else
	fprintf(outenergy, "%d\t%f\t%f\t%f\t%f\n", i, sqr,
		(force[3 * natoms + i] + force_0[3 * natoms + i]) / eweight,
		force_0[3 * natoms + i] / eweight,
		force[3 * natoms + i] / force_0[3 * natoms + i]);
    }
    if (strcmp(endenergy, "stdout") != 0) {
      printf("Energy data written to %s\n", endenergy);
      fclose(outenergy);
    }
#ifdef STRESS
    if (strcmp(endstress, "stdout") != 0) {
      outstress = fopen(endstress, "w");
      if (NULL == outstress) {
	sprintf(msg, "Could not open file %s\n", endstress);
	error(msg);
      }
    } else
      outstress = stdout;
    fprintf(outstress, "Stresses on unit cell\n");
    fprintf(outstress, "conf\t(w*ds)^2\ts\t\ts0\t\tds/s0\n");
    for (i = 3 * natoms + nconf; i < 3 * natoms + 7 * nconf; i++) {
      sqr = SQR(force[i]);
      s_sum += sqr;
      max = MAX(max, sqr);
      min = MIN(min, sqr);
      fprintf(outstress, "%d\t%f\t%f\t%f\t%f\n",
	      (i - (3 * natoms + nconf)) / 6, sqr,
	      (force[i] + force_0[i]) / sweight, force_0[i] / sweight,
	      force[i] / force_0[i]);
    }
#endif
#ifdef EAM
    printf("Punishment Constraints\n");
//    printf("conf dp p p0 dp/p0");
#ifdef STRESS
    diff = 6 * nconf;
#else
    diff = 0;
#endif
    for (i = 3 * natoms + nconf + diff; i < 3 * natoms + 2 * nconf + diff;
	 i++) {
      sqr = SQR(force[i]);
      max = MAX(max, sqr);
      min = MIN(min, sqr);
      printf("%d %f %f %f %f\n", i - (3 * natoms + nconf + diff), sqr,
	     force[i] + force_0[i], force_0[i], force[i] / force_0[i]);
    }
    printf("Dummy Constraints\n");
    for (i = 2 * ntypes; i > 0; i--) {
      sqr = SQR(force[mdim - i]);
      max = MAX(max, sqr);
      min = MIN(min, sqr);
      printf("%d %f %f %f %f\n", ntypes - i, sqr,
	     force[mdim - i] + force_0[mdim - i],
	     force_0[mdim - i], force[mdim - i] / force_0[mdim - i]);
    }
#endif

    printf("\n###### error report #####\n");
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

    for (i = 0; i < 3 * natoms; i++)
      rms[0] += SQR(force[i]);
    rms[0] = sqrt(rms[0] / natoms);

    for (i = 0; i < nconf; i++)
      rms[1] += SQR(force[3 * natoms + i] / eweight);
    if (isnan(rms[1]))
      rms[1] = 0;
    rms[1] = sqrt(rms[1] / nconf);

    for (i = 0; i < nconf; i++)
      for (j = 0; j < 6; j++)
	rms[2] += SQR(force[3 * natoms + nconf + 6 * i + j] / sweight);
    if (isnan(rms[2]))
      rms[2] = 0;
    rms[2] = sqrt(rms[2] / (6 * nconf));

#ifndef STRESS
    printf("sum of force-errors = %f\t\t( %.2f%% - av: %e)\n", tot - e_sum,
	   (tot - e_sum) / tot * 100, (tot - e_sum) / (3 * natoms));
    printf("sum of energy-errors = %f\t\t( %.2f%% )\n", e_sum,
	   e_sum / tot * 100);
    printf("min: %e - max: %e\n", min, max);
    printf("rms-errors:\nforce %f\nenergy %f\n", rms[0], rms[1]);
#else
    printf("sum of force-errors = %f\t\t( %.2f%% - av: %e)\n",
	   tot - e_sum - s_sum, (tot - e_sum - s_sum) / tot * 100,
	   (tot - e_sum - s_sum) / (3 * natoms));
    printf("sum of energy-errors = %f\t\t( %.2f%% )\n", e_sum,
	   e_sum / tot * 100);
    printf("sum of stress-errors = %f\t\t( %.2f%% )\n", s_sum,
	   s_sum / tot * 100);
    printf("min: %e - max: %e\n", min, max);
    printf("rms-errors:\nforce %f\nenergy %f\nstress %f\n", rms[0], rms[1],
	   rms[2]);
#endif

#ifdef MPI
    calc_forces(calc_pot.table, force, 1);	/* go wake up other threads */
#endif /* MPI */
#ifndef APOT
    free(ntyp);
    free(totdens);
#endif
    free(lambda);
  }
  free(force);
#ifdef MPI
  /* kill MPI */
  shutdown_mpi();
#endif
  return 0;
}
