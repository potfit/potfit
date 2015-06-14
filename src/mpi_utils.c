/****************************************************************
 *
 * mpi_utils.c: Contains utilities to be used with MPI
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

#include "potfit.h"

#ifdef MPI

#include "config.h"
#include "utils.h"

/****************************************************************
 *
 * set up mpi
 *
 ****************************************************************/

int init_mpi(int* argc, char*** argv)
{
  /* initialize the MPI communication */
  int rval = MPI_Init(argc, argv);

  if (rval != MPI_SUCCESS) {
    printf("Error initializing MPI communication! (Error: %d)\n", rval);
    return -1;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &g_mpi.num_cpus);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_mpi.myid);

  if (g_mpi.myid == 0)
    printf("This is %s compiled on %s, %s.\n\n", POTFIT_VERSION, __DATE__, __TIME__);

  printf("Starting up MPI with %d processes.\n", g_mpi.num_cpus);

  return MPI_SUCCESS;
}


/****************************************************************
 *
 * shut down mpi
 *
 ****************************************************************/

void shutdown_mpi(void)
{
  if (!g_todo.init_done) {
    fprintf(stderr, "MPI will be killed, because the initialization is not yet complete.\n");
    fprintf(stderr, "This is not a bug!\n\n");
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
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

void broadcast_params_mpi()
{
  int   blklens[MAX_MPI_COMPONENTS];
  MPI_Aint displs[MAX_MPI_COMPONENTS];
  MPI_Datatype typen[MAX_MPI_COMPONENTS];
  neigh_t testneigh;
#ifdef THREEBODY
  angle_t testangl;
#endif /* THREEBODY */
  atom_t testatom;
  int   calclen, size, i, each, odd, count;
#ifdef APOT
  int   j;
#endif /* APOT */

  if (g_mpi.myid == 0) {
    printf("Broadcasting data to MPI workers ... ");
    fflush(stdout);
  }

  MPI_Bcast(&g_todo.init_done, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Define Structures */
  /* first the easy ones: */
  /* MPI_VECTOR */
  MPI_Type_contiguous(3, MPI_DOUBLE, &g_mpi.MPI_VECTOR);
  MPI_Type_commit(&g_mpi.MPI_VECTOR);
  /* MPI_STENS */
  MPI_Type_contiguous(6, MPI_DOUBLE, &g_mpi.MPI_STENS);
  MPI_Type_commit(&g_mpi.MPI_STENS);

  /* MPI_NEIGH */
  /* *INDENT-OFF* */
  size = 0;
  blklens[size] = 1;         	typen[size++] = MPI_INT;     	/* type */
  blklens[size] = 1;         	typen[size++] = MPI_INT;     	/* nr */
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;    	/* r */
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;    	/* r2 */
  blklens[size] = 1;         	typen[size++] = MPI_DOUBLE;   	/* inv_r */
  blklens[size] = 1;         	typen[size++] = g_mpi.MPI_VECTOR;  	/* dist */
  blklens[size] = 1;         	typen[size++] = g_mpi.MPI_VECTOR;  	/* dist_r */
  blklens[size] = SLOTS;     	typen[size++] = MPI_INT;    	/* slot */
  blklens[size] = SLOTS;     	typen[size++] = MPI_DOUBLE;     /* shift */
  blklens[size] = SLOTS;     	typen[size++] = MPI_DOUBLE;     /* step */
  blklens[size] = SLOTS;     	typen[size++] = MPI_INT;     	/* col */
#ifdef ADP
blklens[size] = 1;         	typen[size++] = g_mpi.MPI_STENS;   	/* sqrdist */
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
blklens[size] = 1; 		typen[size++] = g_mpi.MPI_VECTOR; 	/* dzeta */
#endif /* TERSOFF */

  count = 0;
  MPI_Address(&testneigh.type, 		&displs[count++]);
  MPI_Address(&testneigh.nr, 		&displs[count++]);
  MPI_Address(&testneigh.r, 		&displs[count++]);
  MPI_Address(&testneigh.r2, 		&displs[count++]);
  MPI_Address(&testneigh.inv_r, 	&displs[count++]);
  MPI_Address(&testneigh.dist, 		&displs[count++]);
  MPI_Address(&testneigh.dist_r,	&displs[count++]);
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

  MPI_Type_create_struct(size, blklens, displs, typen, &g_mpi.MPI_NEIGH);
  MPI_Type_commit(&g_mpi.MPI_NEIGH);

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
  blklens[size] = 1;         	typen[size++] = g_mpi.MPI_VECTOR;  	/* pos */
  blklens[size] = 1;         	typen[size++] = g_mpi.MPI_VECTOR;  	/* force */
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
blklens[size] = 1;         	typen[size++] = g_mpi.MPI_VECTOR;  	/* mu */
blklens[size] = 1;         	typen[size++] = g_mpi.MPI_STENS;   	/* lambda */
  blklens[size] = 1;        	typen[size++] = MPI_DOUBLE;    	/* nu */
#endif /* ADP */
#if defined DIPOLE
blklens[size] = 1;         	typen[size++] = g_mpi.MPI_VECTOR;    	/* E_stat */
blklens[size] = 1;         	typen[size++] = g_mpi.MPI_VECTOR;   	/* p_sr */
blklens[size] = 1;         	typen[size++] = g_mpi.MPI_VECTOR;   	/* E_ind */
blklens[size] = 1;         	typen[size++] = g_mpi.MPI_VECTOR;   	/* p_ind */
blklens[size] = 1;        	typen[size++] = g_mpi.MPI_VECTOR;   	/* E_old */
blklens[size] = 1;        	typen[size++] = g_mpi.MPI_VECTOR;   	/* E_tot */
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
  MPI_Address(&testatom.num_angles, 	&displs[count++]);
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

  MPI_Type_create_struct(size, blklens, displs, typen, &g_mpi.MPI_ATOM);
  MPI_Type_commit(&g_mpi.MPI_ATOM);

  /* Distribute fundamental parameters */
  MPI_Bcast(&g_calc.mdim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_param.ntypes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_config.natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_config.nconf, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_calc.paircol, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_param.opt, 1, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef COULOMB
  MPI_Bcast(&dp_cut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif /* COULOMB */
#ifdef DIPOLE
  MPI_Bcast(&dp_tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dp_mix, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif /* DIPOLE */
  if (g_mpi.myid > 0) {
    g_config.inconf = (int *)malloc(g_config.nconf * sizeof(int));
    g_config.cnfstart = (int *)malloc(g_config.nconf * sizeof(int));
    g_config.force_0 = (double *)malloc(g_calc.mdim * sizeof(double));
    g_config.conf_weight = (double *)malloc(g_config.nconf * sizeof(double));
    reg_for_free(g_config.inconf, "inconf");
    reg_for_free(g_config.cnfstart, "cnfstart");
    reg_for_free(g_config.force_0, "force_0");
    reg_for_free(g_config.conf_weight, "conf_weight");
  }
  MPI_Bcast(g_config.inconf, g_config.nconf, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_config.cnfstart, g_config.nconf, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_config.force_0, g_calc.mdim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_config.conf_weight, g_config.nconf, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* Broadcast weights... */
  MPI_Bcast(&g_param.eweight, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_param.sweight, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* Broadcast the potential... */
  MPI_Bcast(&g_pot.format, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_pot.calc_pot.len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_pot.calc_pot.ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);
  size = g_pot.calc_pot.ncols;
  calclen = g_pot.calc_pot.len;
  if (g_mpi.myid > 0) {
    g_pot.calc_pot.begin = (double *)malloc(size * sizeof(double));
    g_pot.calc_pot.end = (double *)malloc(size * sizeof(double));
    g_pot.calc_pot.step = (double *)malloc(size * sizeof(double));
    g_pot.calc_pot.invstep = (double *)malloc(size * sizeof(double));
    g_pot.calc_pot.first = (int *)malloc(size * sizeof(int));
    g_pot.calc_pot.last = (int *)malloc(size * sizeof(int));
    g_pot.calc_pot.table = (double *)malloc(calclen * sizeof(double));
    g_pot.calc_pot.xcoord = (double *)malloc(calclen * sizeof(double));
    g_pot.calc_pot.d2tab = (double *)malloc(calclen * sizeof(double));
    reg_for_free(g_pot.calc_pot.begin, "calc_pot.begin");
    reg_for_free(g_pot.calc_pot.end, "calc_pot.end");
    reg_for_free(g_pot.calc_pot.step, "calc_pot.step");
    reg_for_free(g_pot.calc_pot.invstep, "calc_pot.invstep");
    reg_for_free(g_pot.calc_pot.first, "calc_pot.first");
    reg_for_free(g_pot.calc_pot.last, "calc_pot.last");
    reg_for_free(g_pot.calc_pot.table, "calc_pot.table");
    reg_for_free(g_pot.calc_pot.xcoord, "calc_pot.xcoord");
    reg_for_free(g_pot.calc_pot.d2tab, "calc_pot.d2tab");
  }
  MPI_Bcast(g_pot.calc_pot.begin, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.calc_pot.end, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.calc_pot.step, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.calc_pot.invstep, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.calc_pot.first, size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.calc_pot.last, size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.calc_pot.table, calclen, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.calc_pot.d2tab, calclen, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.calc_pot.xcoord, calclen, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef APOT
  MPI_Bcast(&g_param.enable_cp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_pot.opt_pot.len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_pot.apot_table.number, 1, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef COULOMB
  MPI_Bcast(&g_pot.apot_table.total_ne_par, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif /* COULOMB */
  if (g_param.enable_cp) {
    if (g_mpi.myid > 0) {
      g_config.na_type = (int **)malloc((g_config.nconf + 1) * sizeof(int *));
      for (i = 0; i < (g_config.nconf + 1); i++) {
        g_config.na_type[i] = (int *)malloc(g_param.ntypes * sizeof(int));
        reg_for_free(g_config.na_type[i], "na_type[%d]", i);
      }
      reg_for_free(g_config.na_type, "na_type");
    }
    for (i = 0; i < (g_config.nconf + 1); i++)
      MPI_Bcast(g_config.na_type[i], g_param.ntypes, MPI_INT, 0, MPI_COMM_WORLD);
  }
  if (g_mpi.myid > 0) {
    g_pot.calc_list = (double *)malloc(g_pot.opt_pot.len * sizeof(double));
    g_pot.apot_table.n_par = (int *)malloc(g_pot.apot_table.number * sizeof(int));
    g_pot.apot_table.begin = (double *)malloc(g_pot.apot_table.number * sizeof(double));
    g_pot.apot_table.end = (double *)malloc(g_pot.apot_table.number * sizeof(double));
    g_pot.apot_table.idxpot = (int *)malloc(g_pot.apot_table.number * sizeof(int));
#ifdef COULOMB
    g_pot.apot_table.ratio = (double *)malloc(g_param.ntypes * sizeof(double));
#endif /* COULOMB */
    g_pot.smooth_pot = (int *)malloc(g_pot.apot_table.number * sizeof(int));
    g_pot.invar_pot = (int *)malloc(g_pot.apot_table.number * sizeof(int));
    g_config.rcut = (double *)malloc(g_param.ntypes * g_param.ntypes * sizeof(double));
    g_config.rmin = (double *)malloc(g_param.ntypes * g_param.ntypes * sizeof(double));
    g_pot.apot_table.fvalue = (fvalue_pointer *) malloc(g_pot.apot_table.number * sizeof(fvalue_pointer));
    g_pot.opt_pot.table = (double *)malloc(g_pot.opt_pot.len * sizeof(double));
    g_pot.opt_pot.first = (int *)malloc(g_pot.apot_table.number * sizeof(int));
    reg_for_free(g_pot.calc_list, "calc_list");
    reg_for_free(g_pot.apot_table.n_par, "apot_table.n_par");
    reg_for_free(g_pot.apot_table.begin, "apot_table.begin");
    reg_for_free(g_pot.apot_table.end, "apot_table.end");
    reg_for_free(g_pot.apot_table.idxpot, "apot_table.idxpot");
#ifdef COULOMB
    reg_for_free(g_pot.apot_table.ratio, "apot_table.ratio");
#endif
    reg_for_free(g_pot.smooth_pot, "smooth_pot");
    reg_for_free(g_pot.invar_pot, "invar_pot");
    reg_for_free(g_config.rcut, "rcut");
    reg_for_free(g_config.rmin, "rmin");
    reg_for_free(g_pot.apot_table.fvalue, "apot_table.fvalue");
    reg_for_free(g_pot.opt_pot.table, "opt_pot.first");
    reg_for_free(g_pot.opt_pot.first, "opt_pot.first");
  }
  MPI_Bcast(g_pot.smooth_pot, g_pot.apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.invar_pot, g_pot.apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.calc_list, g_pot.opt_pot.len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.apot_table.n_par, g_pot.apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_config.rcut, g_param.ntypes * g_param.ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_config.rmin, g_param.ntypes * g_param.ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.apot_table.fvalue, g_pot.apot_table.number, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.apot_table.end, g_pot.apot_table.number, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.apot_table.begin, g_pot.apot_table.number, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.apot_table.idxpot, g_pot.apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_pot.cp_start, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_pot.have_globals, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_pot.global_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_pot.apot_table.globals, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.opt_pot.first, g_pot.apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef COULOMB
  MPI_Bcast(&g_pot.apot_table.last_charge, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(g_pot.apot_table.ratio, g_param.ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif /* COULOMB */
  if (g_pot.have_globals) {
    if (g_mpi.myid > 0) {
      g_pot.apot_table.n_glob = (int *)malloc(g_pot.apot_table.globals * sizeof(int));
      g_pot.apot_table.global_idx = (int ***)malloc(g_pot.apot_table.globals * sizeof(int **));
      reg_for_free(g_pot.apot_table.n_glob, "apot_table.n_glob");
      reg_for_free(g_pot.apot_table.global_idx, "apot_table.global_idx");
    }
    MPI_Bcast(g_pot.apot_table.n_glob, g_pot.apot_table.globals, MPI_INT, 0, MPI_COMM_WORLD);
    if (g_mpi.myid > 0) {
      for (i = 0; i < g_pot.apot_table.globals; i++) {
        g_pot.apot_table.global_idx[i] = (int **)malloc(g_pot.apot_table.n_glob[i] * sizeof(int *));
        reg_for_free(g_pot.apot_table.global_idx[i], "apot_table.global_idx[%d]", i);
      }
      for (i = 0; i < g_pot.apot_table.globals; i++) {
        for (j = 0; j < g_pot.apot_table.n_glob[i]; j++) {
          g_pot.apot_table.global_idx[i][j] = (int *)malloc(2 * sizeof(int));
          reg_for_free(g_pot.apot_table.global_idx[i][j], "apot_table.global_idx[%d][%d]", i, j);
	}
      }
    }
    for (i = 0; i < g_pot.apot_table.globals; i++)
      for (j = 0; j < g_pot.apot_table.n_glob[i]; j++)
        MPI_Bcast(g_pot.apot_table.global_idx[i][j], 2, MPI_INT, 0, MPI_COMM_WORLD);
  }
#endif /* APOT */

  /* Distribute configurations */
  /* Each node: nconf/num_cpus configurations.
     Last nconf%num_cpus nodes: 1 additional config */
  each = (g_config.nconf / g_mpi.num_cpus);
  odd = (g_config.nconf % g_mpi.num_cpus) - g_mpi.num_cpus;
  if (g_mpi.myid == 0) {
    g_mpi.atom_len = (int *)malloc(g_mpi.num_cpus * sizeof(int));
    g_mpi.atom_dist = (int *)malloc(g_mpi.num_cpus * sizeof(int));
    g_mpi.conf_len = (int *)malloc(g_mpi.num_cpus * sizeof(int));
    g_mpi.conf_dist = (int *)malloc(g_mpi.num_cpus * sizeof(int));
    for (i = 0; i < g_mpi.num_cpus; i++)
      g_mpi.conf_dist[i] = i * each + (((i + odd) > 0) ? (i + odd) : 0);
    for (i = 0; i < g_mpi.num_cpus - 1; i++)
      g_mpi.conf_len[i] = g_mpi.conf_dist[i + 1] - g_mpi.conf_dist[i];
    g_mpi.conf_len[g_mpi.num_cpus - 1] = g_config.nconf - g_mpi.conf_dist[g_mpi.num_cpus - 1];
    for (i = 0; i < g_mpi.num_cpus; i++)
      g_mpi.atom_dist[i] = g_config.cnfstart[i * each + (((i + odd) > 0) ? (i + odd) : 0)];
    for (i = 0; i < g_mpi.num_cpus - 1; i++)
      g_mpi.atom_len[i] = g_mpi.atom_dist[i + 1] - g_mpi.atom_dist[i];
    g_mpi.atom_len[g_mpi.num_cpus - 1] = g_config.natoms - g_mpi.atom_dist[g_mpi.num_cpus - 1];
    reg_for_free(g_mpi.atom_len, "atom_len");
    reg_for_free(g_mpi.atom_dist, "atom_dist");
    reg_for_free(g_mpi.conf_len, "conf_len");
    reg_for_free(g_mpi.conf_dist, "conf_dist");
  }
  MPI_Scatter(g_mpi.atom_len, 1, MPI_INT, &g_mpi.myatoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(g_mpi.atom_dist, 1, MPI_INT, &g_mpi.firstatom, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(g_mpi.conf_len, 1, MPI_INT, &g_mpi.myconf, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(g_mpi.conf_dist, 1, MPI_INT, &g_mpi.firstconf, 1, MPI_INT, 0, MPI_COMM_WORLD);
  /* this broadcasts all atoms */
  g_config.conf_atoms = (atom_t *)malloc(g_mpi.myatoms * sizeof(atom_t));
  for (i = 0; i < g_config.natoms; i++) {
    if (g_mpi.myid == 0)
      testatom = g_config.atoms[i];
    MPI_Bcast(&testatom, 1, g_mpi.MPI_ATOM, 0, MPI_COMM_WORLD);
    if (i >= g_mpi.firstatom && i < (g_mpi.firstatom + g_mpi.myatoms)) {
      g_config.conf_atoms[i - g_mpi.firstatom] = testatom;
    }
  }
  broadcast_neighbors();
#ifdef THREEBODY
  broadcast_angles();
#endif /* THREEBODY */
  g_config.conf_vol = (double *)malloc(g_mpi.myconf * sizeof(double));
  g_config.conf_uf = (int *)malloc(g_mpi.myconf * sizeof(int));
#ifdef STRESS
  g_config.conf_us = (int *)malloc(g_mpi.myconf * sizeof(double));
#endif /* STRESS */
  MPI_Scatterv(g_config.volume, g_mpi.conf_len, g_mpi.conf_dist, MPI_DOUBLE, g_config.conf_vol, g_mpi.myconf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(g_config.useforce, g_mpi.conf_len, g_mpi.conf_dist, MPI_INT, g_config.conf_uf, g_mpi.myconf, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef STRESS
  MPI_Scatterv(g_config.usestress, g_mpi.conf_len, g_mpi.conf_dist, MPI_INT, g_config.conf_us, g_mpi.myconf, MPI_INT, 0, MPI_COMM_WORLD);
#endif /* STRESS */

  reg_for_free(g_config.conf_vol, "conf_vol");
  reg_for_free(g_config.conf_uf, "conf_uf");
#ifdef STRESS
  reg_for_free(g_config.conf_us, "conf_us");
#endif /* STRESS */
  reg_for_free(g_config.conf_atoms, "conf_atoms");

  if (g_mpi.myid == 0) {
    printf("done\n");
    fflush(stdout);
  }
}

/****************************************************************
 *
 * scatter dynamic neighbor table
 *
 ****************************************************************/

void broadcast_neighbors()
{
  int   i, j, neighs = 0;
  neigh_t neigh;
  atom_t *atom;

  init_neigh_memory(&neigh);

  for (i = 0; i < g_config.natoms; i++) {
    atom = g_config.conf_atoms + i - g_mpi.firstatom;
    if (g_mpi.myid == 0)
      neighs = g_config.atoms[i].num_neigh;
    MPI_Bcast(&neighs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (i >= g_mpi.firstatom && i < (g_mpi.firstatom + g_mpi.myatoms)) {
      atom->neigh = (neigh_t *)malloc(neighs * sizeof(neigh_t));
      for (j = 0; j < neighs; j++)
	init_neigh_memory(atom->neigh + j);
      reg_for_free(atom->neigh, "broadcast atom[%d]->neigh", i);
    }
    for (j = 0; j < neighs; j++) {
      if (g_mpi.myid == 0)
        neigh = g_config.atoms[i].neigh[j];
      MPI_Bcast(&neigh, 1, g_mpi.MPI_NEIGH, 0, MPI_COMM_WORLD);
      if (i >= g_mpi.firstatom && i < (g_mpi.firstatom + g_mpi.myatoms)) {
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
  int   i, j, nangles = 0;
  angle_t angle;
  atom_t *atom;

  init_angle(&angle);

  for (i = 0; i < natoms; ++i) {
    atom = conf_atoms + i - firstatom;
    if (myid == 0)
      nangles = atoms[i].num_angles;
    MPI_Bcast(&nangles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (i >= firstatom && i < (firstatom + myatoms)) {
      atom->angle_part = (angle_t *) malloc(nangles * sizeof(angle_t));
      for (j = 0; j < nangles; j++)
	init_angle(atom->angle_part + j);
      reg_for_free(atom->angle_part, "broadcast atom[%d]->angle_part", i);
    }
    for (j = 0; j < nangles; ++j) {
      if (myid == 0)
	angle = atoms[i].angle_part[j];
      MPI_Bcast(&angle, 1, MPI_ANGL, 0, MPI_COMM_WORLD);
      if (i >= firstatom && i < (firstatom + myatoms)) {
	atom->angle_part[j] = angle;
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
  firstcol = paircol + g_param.ntypes;
  /* Memory is allocated - just bcast that changed potential... */
  /* bcast begin/end/step/invstep of embedding energy  */
  MPI_Bcast(calc_pot.begin + firstcol, g_param.ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.end + firstcol, g_param.ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.step + firstcol, g_param.ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.invstep + firstcol, g_param.ntypes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.first + firstcol, g_param.ntypes, MPI_INT, 0, MPI_COMM_WORLD);
  /* bcast table values of transfer fn. and embedding energy */
  firstval = g_pot.calc_pot.first[paircol];
  nvals = g_pot.calc_pot.len - firstval;
  MPI_Bcast(calc_pot.table + firstval, nvals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

#endif /* !APOT */

#endif /* MPI */
