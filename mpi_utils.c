/****************************************************************
 *
 * mpi_utils.c: Contains utilities to be used with MPI
 *
 ****************************************************************
 *
 * Copyright 2002-2013
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

#ifdef MPI

#include "potfit.h"

#include "utils.h"

/****************************************************************
 *
 * set up mpi
 *
 ****************************************************************/

void init_mpi(int argc, char **argv)
{
  /* Initialize MPI */
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS && myid == 0)
    fprintf(stderr, "MPI_Init failed!\n");
  MPI_Comm_size(MPI_COMM_WORLD, &num_cpus);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}


/****************************************************************
 *
 * shut down mpi
 *
 ****************************************************************/

void shutdown_mpi(void)
{
  if (!init_done) {
    fprintf(stderr, "MPI will be killed, because the initialization is not yet complete.\n");
    fprintf(stderr, "This is not a bug!\n\n");
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  MPI_Barrier(MPI_COMM_WORLD);	/* Wait for all processes to arrive */
  MPI_Finalize();		/* Shutdown */
}

/****************************************************************
 *
 * broadcast_param: Broadcast parameters etc to other nodes
 *
 ****************************************************************/

#define MAX_MPI_COMPONENTS 20

void broadcast_params()
{
  int   blklens[MAX_MPI_COMPONENTS];
  MPI_Aint displs[MAX_MPI_COMPONENTS];
  MPI_Datatype typen[MAX_MPI_COMPONENTS];
  neigh_t testneigh;
#ifdef THREEBODY
  angl  testangl;
#endif /* THREEBODY */
  atom_t testatom;
  int   calclen, size, i, j, each, odd, count;

  /* Define Structures */
  /* first the easy ones: */
  /* MPI_VECTOR */
  MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_VECTOR);
  MPI_Type_commit(&MPI_VECTOR);
  /* MPI_STENS */
  MPI_Type_contiguous(6, MPI_DOUBLE, &MPI_STENS);
  MPI_Type_commit(&MPI_STENS);

  /* MPI_NEIGH */
  /* *INDENT-OFF* */
  size = 0;
  blklens[size] = 1;         	typen[size++] = MPI_INT;     	/* type */
  blklens[size] = 1;         	typen[size++] = MPI_INT;     	/* nr */
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;    	/* r */
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;    	/* r2 */
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;   	/* inv_r */
  blklens[size] = 1;         	typen[size++] = MPI_VECTOR;  	/* dist_r */
  blklens[size] = 1;         	typen[size++] = MPI_VECTOR;  	/* rdist */
  blklens[size] = SLOTS;     	typen[size++] = MPI_INT;    	/* slot */
  blklens[size] = SLOTS;     	typen[size++] = MPI_DOUBLE;     /* shift */
  blklens[size] = SLOTS;     	typen[size++] = MPI_DOUBLE;     /* step */
  blklens[size] = SLOTS;     	typen[size++] = MPI_INT;     	/* col */
#ifdef ADP
  blklens[size] = 1;         	typen[size++] = MPI_STENS;   	/* sqrdist */
  blklens[size] = 1;        	typen[size++] = MPI_DOUBLE;     /* u_val */
  blklens[size] = 1;        	typen[size++] = MPI_DOUBLE;     /* u_grad */
  blklens[size] = 1;        	typen[size++] = MPI_DOUBLE;     /* w_val */
  blklens[size] = 1;        	typen[size++] = MPI_DOUBLE;     /* w_grad */
#endif /* ADP */
#ifdef COULOMB
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;     /* fnval_el */
  blklens[size] = 1;        	typen[size++] = MPI_DOUBLE;     /* grad_el */
  blklens[size] = 1;        	typen[size++] = MPI_DOUBLE;     /* ggrad_el */
#endif /* COULOMB */
#ifdef THREEBODY
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;     /* f */
  blklens[size] = 1;        	typen[size++] = MPI_DOUBLE;     /* df */
  blklens[size] = 1;        	typen[size++] = MPI_INT;     	/* ijk_start */
#endif /* THREEBODY */
#ifdef MEAM
  blklens[size] = 1; 		typen[size++] = MPI_DOUBLE; 	/* drho */
#endif /* MEAM */
#ifdef TERSOFF
  blklens[size] = 1; 		typen[size++] = MPI_VECTOR; 	/* dzeta */
#endif /* TERSOFF */

  count = 0;
  MPI_Address(&testneigh.type, 		&displs[count++]);
  MPI_Address(&testneigh.nr, 		&displs[count++]);
  MPI_Address(&testneigh.r, 		&displs[count++]);
  MPI_Address(&testneigh.r2, 		&displs[count++]);
  MPI_Address(&testneigh.inv_r, 	&displs[count++]);
  MPI_Address(&testneigh.dist_r,	&displs[count++]);
  MPI_Address(&testneigh.rdist, 	&displs[count++]);
  MPI_Address(testneigh.slot, 		&displs[count++]);
  MPI_Address(testneigh.shift, 		&displs[count++]);
  MPI_Address(testneigh.step, 		&displs[count++]);
  MPI_Address(testneigh.col, 		&displs[count++]);
#ifdef ADP
  MPI_Address(&testneigh.sqrdist, 	&displs[count++]);
  MPI_Address(&testneigh.u_val, 	&displs[count++]);
  MPI_Address(&testneigh.u_grad, 	&displs[count++]);
  MPI_Address(&testneigh.w_val, 	&displs[count++]);
  MPI_Address(&testneigh.w_grad, 	&displs[count++]);
#endif /* ADP */
#ifdef COULOMB
  MPI_Address(&testneigh.fnval_el, 	&displs[count++]);
  MPI_Address(&testneigh.grad_el, 	&displs[count++]);
  MPI_Address(&testneigh.ggrad_el, 	&displs[count++]);
#endif /* COULOMB */
#ifdef THREEBODY
  MPI_Address(&testneigh.f, 		&displs[count++]);
  MPI_Address(&testneigh.df, 		&displs[count++]);
  MPI_Address(&testneigh.ijk_start, 	&displs[count++]);
#endif /* THREEBODY */
#ifdef MEAM
  MPI_Address(&testneigh.drho, 		&displs[count++]);
#endif /* MEAM */
#ifdef TERSOFF
  MPI_Address(&testneigh.dzeta, 	&displs[count++]);
#endif /* MEAM */

  /* *INDENT-ON* */

  /* set displacements */
  for (i = 1; i < count; i++) {
    displs[i] -= displs[0];
  }
  displs[0] = 0;

  MPI_Type_create_struct(size, blklens, displs, typen, &MPI_NEIGH);
  MPI_Type_commit(&MPI_NEIGH);

#ifdef THREEBODY
  /* MPI_ANGL */
  /* *INDENT-OFF* */
  size = 0;
  blklens[size] = 1; 		typen[size++] = MPI_DOUBLE;    	/* cos */
#ifdef MEAM
  blklens[size] = 1;         	typen[size++] = MPI_INT;     	/* slot */
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;    	/* shift */
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;    	/* step */
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;   	/* g */
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;    	/* dg */
#endif /* MEAM */

  count = 0;
  MPI_Address(&testangl.cos, 		&displs[count++]);
#ifdef MEAM
  MPI_Address(&testangl.slot, 		&displs[count++]);
  MPI_Address(&testangl.shift, 		&displs[count++]);
  MPI_Address(&testangl.step, 		&displs[count++]);
  MPI_Address(&testangl.g, 		&displs[count++]);
  MPI_Address(&testangl.dg, 		&displs[count++]);
#endif /* MEAM */
  /* *INDENT-ON* */

  /* set displacements */
  for (i = 1; i < count; i++) {
    displs[i] -= displs[0];
  }
  displs[0] = 0;

  MPI_Type_struct(size, blklens, displs, typen, &MPI_ANGL);
  MPI_Type_commit(&MPI_ANGL);
#endif /* THREEBODY */

  /* MPI_ATOM */
  /* *INDENT-OFF* */
  size = 0;
  blklens[size] = 1;         	typen[size++] = MPI_INT;     	/* type */
  blklens[size] = 1;         	typen[size++] = MPI_INT;     	/* num_neigh */
  blklens[size] = 1;         	typen[size++] = MPI_VECTOR;  	/* pos */
  blklens[size] = 1;         	typen[size++] = MPI_VECTOR;  	/* force */
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;    	/* absforce */
  blklens[size] = 1;         	typen[size++] = MPI_INT;     	/* conf */
#ifdef CONTRIB
  blklens[size] = 1;         	typen[size++] = MPI_INT;     	/* contrib */
#endif /* CONTRIB */
#if defined EAM || defined ADP || defined MEAM
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;    	/* rho */
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;    	/* gradF */
#endif /* EAM || ADP || MEAM */
#ifdef ADP
  blklens[size] = 1;         	typen[size++] = MPI_VECTOR;  	/* mu */
  blklens[size] = 1;         	typen[size++] = MPI_STENS;   	/* lambda */
  blklens[size] = 1;        	typen[size++] = MPI_DOUBLE;    	/* nu */
#endif /* ADP */
#if defined DIPOLE
  blklens[size] = 1;         	typen[size++] = MPI_VECTOR;    	/* E_stat */
  blklens[size] = 1;         	typen[size++] = MPI_VECTOR;   	/* p_sr */
  blklens[size] = 1;         	typen[size++] = MPI_VECTOR;   	/* E_ind */
  blklens[size] = 1;         	typen[size++] = MPI_VECTOR;   	/* p_ind */
  blklens[size] = 1;        	typen[size++] = MPI_VECTOR;   	/* E_old */
  blklens[size] = 1;        	typen[size++] = MPI_VECTOR;   	/* E_tot */
#endif /* DIPOLE */
#ifdef THREEBODY
  blklens[size] = 1;         	typen[size++] = MPI_INT;    	/* num_angl */
#ifdef MEAM
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;    	/* rho_eam */
#endif /* MEAM */
#endif /* THREEBODY */

  /* DO NOT BROADCAST NEIGHBORS !!! DYNAMIC ALLOCATION */
  /* DO NOT BROADCAST ANGLES !!! DYNAMIC ALLOCATION */

  count = 0;
  MPI_Address(&testatom.type, 		&displs[count++]);
  MPI_Address(&testatom.num_neigh, 	&displs[count++]);
  MPI_Address(&testatom.pos, 		&displs[count++]);
  MPI_Address(&testatom.force, 		&displs[count++]);
  MPI_Address(&testatom.absforce, 	&displs[count++]);
  MPI_Address(&testatom.conf, 		&displs[count++]);
#ifdef CONTRIB
  MPI_Address(&testatom.contrib, 	&displs[count++]);
#endif /* CONTRIB */
#if defined EAM || defined ADP || defined MEAM
  MPI_Address(&testatom.rho, 		&displs[count++]);
  MPI_Address(&testatom.gradF, 		&displs[count++]);
#endif /* EAM || ADP */
#ifdef ADP
  MPI_Address(&testatom.mu, 		&displs[count++]);
  MPI_Address(&testatom.lambda, 	&displs[count++]);
  MPI_Address(&testatom.nu, 		&displs[count++]);
#endif /* ADP */
#ifdef DIPOLE
  MPI_Address(&testatom.E_stat, 	&displs[count++]);
  MPI_Address(&testatom.p_sr, 		&displs[count++]);
  MPI_Address(&testatom.E_ind, 		&displs[count++]);
  MPI_Address(&testatom.p_ind, 		&displs[count++]);
  MPI_Address(&testatom.E_old, 		&displs[count++]);
  MPI_Address(&testatom.E_tot, 		&displs[count++]);
#endif /* DIPOLE */
#ifdef THREEBODY
  MPI_Address(&testatom.num_angl, 	&displs[count++]);
#ifdef MEAM
  MPI_Address(&testatom.rho_eam,	&displs[count++]);
#endif /* MEAM */
#endif /* THREEBODY */

  /* *INDENT-ON* */

  /* set displacements */
  for (i = 1; i < count; i++) {
    displs[i] -= displs[0];
  }
  displs[0] = 0;

  MPI_Type_create_struct(size, blklens, displs, typen, &MPI_ATOM);
  MPI_Type_commit(&MPI_ATOM);

  /* Distribute fundamental parameters */
  MPI_Bcast(&mdim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ntypes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nconf, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opt, 1, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef COULOMB
  MPI_Bcast(&dp_cut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif /* COULOMB */
#ifdef DIPOLE
  MPI_Bcast(&dp_tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dp_mix, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif /* DIPOLE */
  if (myid > 0) {
    inconf = (int *)malloc(nconf * sizeof(int));
    cnfstart = (int *)malloc(nconf * sizeof(int));
    force_0 = (double *)malloc(mdim * sizeof(double));
    conf_weight = (double *)malloc(nconf * sizeof(double));
    reg_for_free(inconf, "inconf");
    reg_for_free(cnfstart, "cnfstart");
    reg_for_free(force_0, "force_0");
    reg_for_free(conf_weight, "conf_weight");
  }
  MPI_Bcast(inconf, nconf, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(cnfstart, nconf, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(force_0, mdim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(conf_weight, nconf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&maxneigh, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Broadcast weights... */
  MPI_Bcast(&eweight, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sweight, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* Broadcast the potential... */
  MPI_Bcast(&format, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&calc_pot.len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&calc_pot.ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);
  size = calc_pot.ncols;
  calclen = calc_pot.len;
  if (myid > 0) {
    calc_pot.begin = (double *)malloc(size * sizeof(double));
    calc_pot.end = (double *)malloc(size * sizeof(double));
    calc_pot.step = (double *)malloc(size * sizeof(double));
    calc_pot.invstep = (double *)malloc(size * sizeof(double));
    calc_pot.first = (int *)malloc(size * sizeof(int));
    calc_pot.last = (int *)malloc(size * sizeof(int));
    calc_pot.table = (double *)malloc(calclen * sizeof(double));
    calc_pot.xcoord = (double *)malloc(calclen * sizeof(double));
    calc_pot.d2tab = (double *)malloc(calclen * sizeof(double));
    reg_for_free(calc_pot.begin, "calc_pot.begin");
    reg_for_free(calc_pot.end, "calc_pot.end");
    reg_for_free(calc_pot.step, "calc_pot.step");
    reg_for_free(calc_pot.invstep, "calc_pot.invstep");
    reg_for_free(calc_pot.first, "calc_pot.first");
    reg_for_free(calc_pot.last, "calc_pot.last");
    reg_for_free(calc_pot.table, "calc_pot.table");
    reg_for_free(calc_pot.xcoord, "calc_pot.xcoord");
    reg_for_free(calc_pot.d2tab, "calc_pot.d2tab");
  }
  MPI_Bcast(calc_pot.begin, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.end, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.step, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.invstep, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.first, size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.last, size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.table, calclen, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.d2tab, calclen, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.xcoord, calclen, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef APOT
  MPI_Bcast(&enable_cp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opt_pot.len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&apot_table.number, 1, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef COULOMB
  MPI_Bcast(&apot_table.total_ne_par, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif /* COULOMB */
  if (enable_cp) {
    if (myid > 0) {
      na_type = (int **)malloc((nconf + 1) * sizeof(int *));
      for (i = 0; i < (nconf + 1); i++) {
	na_type[i] = (int *)malloc(ntypes * sizeof(int));
	reg_for_free(na_type[i], "na_type[%d]", i);
      }
      reg_for_free(na_type, "na_type");
    }
    for (i = 0; i < (nconf + 1); i++)
      MPI_Bcast(na_type[i], ntypes, MPI_INT, 0, MPI_COMM_WORLD);
  }
  if (myid > 0) {
    calc_list = (double *)malloc(opt_pot.len * sizeof(double));
    apot_table.n_par = (int *)malloc(apot_table.number * sizeof(int));
    apot_table.begin = (double *)malloc(apot_table.number * sizeof(double));
    apot_table.end = (double *)malloc(apot_table.number * sizeof(double));
    apot_table.idxpot = (int *)malloc(apot_table.number * sizeof(int));
#ifdef COULOMB
    apot_table.ratio = (double *)malloc(ntypes * sizeof(double));
#endif /* COULOMB */
    smooth_pot = (int *)malloc(apot_table.number * sizeof(int));
    invar_pot = (int *)malloc(apot_table.number * sizeof(int));
    rcut = (double *)malloc(ntypes * ntypes * sizeof(double));
    rmin = (double *)malloc(ntypes * ntypes * sizeof(double));
    apot_table.fvalue = (fvalue_pointer *) malloc(apot_table.number * sizeof(fvalue_pointer));
    opt_pot.table = (double *)malloc(opt_pot.len * sizeof(double));
    opt_pot.first = (int *)malloc(apot_table.number * sizeof(int));
    reg_for_free(calc_list, "calc_list");
    reg_for_free(apot_table.n_par, "apot_table.n_par");
    reg_for_free(apot_table.begin, "apot_table.begin");
    reg_for_free(apot_table.end, "apot_table.end");
    reg_for_free(apot_table.idxpot, "apot_table.idxpot");
#ifdef COULOMB
    reg_for_free(apot_table.ratio, "apot_table.ratio");
#endif
    reg_for_free(smooth_pot, "smooth_pot");
    reg_for_free(invar_pot, "invar_pot");
    reg_for_free(rcut, "rcut");
    reg_for_free(rmin, "rmin");
    reg_for_free(apot_table.fvalue, "apot_table.fvalue");
    reg_for_free(opt_pot.table, "opt_pot.first");
    reg_for_free(opt_pot.first, "opt_pot.first");
  }
  MPI_Bcast(smooth_pot, apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(invar_pot, apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_list, opt_pot.len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(apot_table.n_par, apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(rcut, ntypes * ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(rmin, ntypes * ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(apot_table.fvalue, apot_table.number, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(apot_table.end, apot_table.number, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(apot_table.begin, apot_table.number, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(apot_table.idxpot, apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cp_start, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&have_globals, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&global_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&apot_table.globals, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(opt_pot.first, apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef COULOMB
  MPI_Bcast(&apot_table.last_charge, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(apot_table.ratio, ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif /* COULOMB */
  if (have_globals) {
    if (myid > 0) {
      apot_table.n_glob = (int *)malloc(apot_table.globals * sizeof(int));
      apot_table.global_idx = (int ***)malloc(apot_table.globals * sizeof(int **));
      reg_for_free(apot_table.n_glob, "apot_table.n_glob");
      reg_for_free(apot_table.global_idx, "apot_table.global_idx");
    }
    MPI_Bcast(apot_table.n_glob, apot_table.globals, MPI_INT, 0, MPI_COMM_WORLD);
    if (myid > 0) {
      for (i = 0; i < apot_table.globals; i++) {
	apot_table.global_idx[i] = (int **)malloc(apot_table.n_glob[i] * sizeof(int *));
	reg_for_free(apot_table.global_idx[i], "apot_table.global_idx[%d]", i);
      }
      for (i = 0; i < apot_table.globals; i++) {
	for (j = 0; j < apot_table.n_glob[i]; j++) {
	  apot_table.global_idx[i][j] = (int *)malloc(2 * sizeof(int));
	  reg_for_free(apot_table.global_idx[i][j], "apot_table.global_idx[%d][%d]", i, j);
	}
      }
    }
    for (i = 0; i < apot_table.globals; i++)
      for (j = 0; j < apot_table.n_glob[i]; j++)
	MPI_Bcast(apot_table.global_idx[i][j], 2, MPI_INT, 0, MPI_COMM_WORLD);
  }
#endif /* APOT */

  /* Distribute configurations */
  /* Each node: nconf/num_cpus configurations.
     Last nconf%num_cpus nodes: 1 additional config */
  each = (nconf / num_cpus);
  odd = (nconf % num_cpus) - num_cpus;
  if (myid == 0) {
    atom_len = (int *)malloc(num_cpus * sizeof(int));
    atom_dist = (int *)malloc(num_cpus * sizeof(int));
    conf_len = (int *)malloc(num_cpus * sizeof(int));
    conf_dist = (int *)malloc(num_cpus * sizeof(int));
    for (i = 0; i < num_cpus; i++)
      conf_dist[i] = i * each + (((i + odd) > 0) ? (i + odd) : 0);
    for (i = 0; i < num_cpus - 1; i++)
      conf_len[i] = conf_dist[i + 1] - conf_dist[i];
    conf_len[num_cpus - 1] = nconf - conf_dist[num_cpus - 1];
    for (i = 0; i < num_cpus; i++)
      atom_dist[i] = cnfstart[i * each + (((i + odd) > 0) ? (i + odd) : 0)];
    for (i = 0; i < num_cpus - 1; i++)
      atom_len[i] = atom_dist[i + 1] - atom_dist[i];
    atom_len[num_cpus - 1] = natoms - atom_dist[num_cpus - 1];
    reg_for_free(atom_len, "atom_len");
    reg_for_free(atom_dist, "atom_dist");
    reg_for_free(conf_len, "conf_len");
    reg_for_free(conf_dist, "conf_dist");
  }
  MPI_Scatter(atom_len, 1, MPI_INT, &myatoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(atom_dist, 1, MPI_INT, &firstatom, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(conf_len, 1, MPI_INT, &myconf, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(conf_dist, 1, MPI_INT, &firstconf, 1, MPI_INT, 0, MPI_COMM_WORLD);
  /* this broadcasts all atoms */
  conf_atoms = (atom_t *)malloc(myatoms * sizeof(atom_t));
  for (i = 0; i < natoms; i++) {
    if (myid == 0)
      testatom = atoms[i];
    MPI_Bcast(&testatom, 1, MPI_ATOM, 0, MPI_COMM_WORLD);
    if (i >= firstatom && i < (firstatom + myatoms)) {
      conf_atoms[i - firstatom] = testatom;
    }
  }
  broadcast_neighbors();
#ifdef THREEBODY
  broadcast_angles();
#endif /* THREEBODY */
  conf_vol = (double *)malloc(myconf * sizeof(double));
  conf_uf = (int *)malloc(myconf * sizeof(double));
  conf_us = (int *)malloc(myconf * sizeof(double));
  MPI_Scatterv(volumen, conf_len, conf_dist, MPI_DOUBLE, conf_vol, myconf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(useforce, conf_len, conf_dist, MPI_INT, conf_uf, myconf, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(usestress, conf_len, conf_dist, MPI_INT, conf_us, myconf, MPI_INT, 0, MPI_COMM_WORLD);

  reg_for_free(conf_vol, "conf_vol");
  reg_for_free(conf_uf, "conf_uf");
  reg_for_free(conf_us, "conf_us");
  reg_for_free(conf_atoms, "conf_atoms");
}

/****************************************************************
 *
 * scatter dynamic neighbor table
 *
 ****************************************************************/

void broadcast_neighbors()
{
  int   i, j, neighs;
  neigh_t neigh;
  atom_t *atom;

  for (i = 0; i < natoms; i++) {
    atom = conf_atoms + i - firstatom;
    if (myid == 0)
      neighs = atoms[i].num_neigh;
    MPI_Bcast(&neighs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (i >= firstatom && i < (firstatom + myatoms)) {
      atom->neigh = (neigh_t *)malloc(neighs * sizeof(neigh_t));
      reg_for_free(atom->neigh, "broadcast atom[%d]->neigh", i);
    }
    for (j = 0; j < neighs; j++) {
      if (myid == 0)
	neigh = atoms[i].neigh[j];
      MPI_Bcast(&neigh, 1, MPI_NEIGH, 0, MPI_COMM_WORLD);
      if (i >= firstatom && i < (firstatom + myatoms)) {
	atom->neigh[j] = neigh;
      }
    }
  }
}

#ifdef THREEBODY

/***************************************************************************
 *
 * scatter dynamic angle table
 *
 **************************************************************************/

void broadcast_angles()
{
  int   i, j, nangles;
  angl  angle;
  atom_t *atom;

  for (i = 0; i < natoms; ++i) {
    atom = conf_atoms + i - firstatom;
    if (myid == 0)
      nangles = atoms[i].num_angl;
    MPI_Bcast(&nangles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (i >= firstatom && i < (firstatom + myatoms)) {
      atom->angl_part = (angl *) malloc(nangles * sizeof(angl));
      reg_for_free(atom->angl_part, "broadcast atom[%d]->angle_part", i);
    }
    for (j = 0; j < nangles; ++j) {
      if (myid == 0)
	angle = atoms[i].angl_part[j];
      MPI_Bcast(&angle, 1, MPI_ANGL, 0, MPI_COMM_WORLD);
      if (i >= firstatom && i < (firstatom + myatoms)) {
	atom->angl_part[j] = angle;
      }
    }
  }
}

#endif /* THREEBODY */

#ifndef APOT

/****************************************************************
 *
 * potsync: Broadcast parameters etc to other nodes
 *
 ****************************************************************/

void potsync()
{
  int   firstcol, firstval, nvals;
  firstcol = paircol + ntypes;
  /* Memory is allocated - just bcast that changed potential... */
  /* bcast begin/end/step/invstep of embedding energy  */
  MPI_Bcast(calc_pot.begin + firstcol, ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.end + firstcol, ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.step + firstcol, ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.invstep + firstcol, ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.first + firstcol, ntypes, MPI_INT, 0, MPI_COMM_WORLD);
  /* bcast table values of transfer fn. and embedding energy */
  firstval = calc_pot.first[paircol];
  nvals = calc_pot.len - firstval;
  MPI_Bcast(calc_pot.table + firstval, nvals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

#endif /* !APOT */

#endif /* MPI */
