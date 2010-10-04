/****************************************************************
 *
 * mpi_utils.c: Contains utilities to be used with MPI
 *
 ****************************************************************
 *
 * Copyright 2004-2010 Peter Brommer, Franz G"ahler, Daniel Schopf
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://www.itap.physik.uni-stuttgart.de/
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
    fprintf(stderr,
      "MPI will be killed, because the initialization is not yet complete.\n");
    fprintf(stderr, "This is not a bug!\n\n");
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  MPI_Barrier(MPI_COMM_WORLD);	/* Wait for all processes to arrive */
  MPI_Finalize();		/* Shutdown */
}

void debug_mpi(int i)
{
/*   int j,k; */
/*   MPI_Status status; */
/*   if (myid==0)  */
/*     for (k=1;k<num_cpus;k++) { */
/*       MPI_Send(&i,1,MPI_INT,k,12,MPI_COMM_WORLD); */
/*       printf("Node 0, sent %d to node %d\n",i,k); */
/*     } */
/*   else {  */
/*     MPI_Recv(&j,1,MPI_INT,0,12,MPI_COMM_WORLD,&status); */
/*     printf("Node %d, received %d at point %d.\n",myid,j,i); */
/*   } */
  MPI_Barrier(MPI_COMM_WORLD);
  printf("I am node %d. Still running at point %d.\n", myid, i);
  fflush(stdout);
}

/****************************************************************
 *
 * broadcast_param: Broadcast parameters etc to other nodes
 *
 ****************************************************************/

#ifdef PAIR
#define MAX_MPI_COMPONENTS 8
#elif defined EAM
#define MAX_MPI_COMPONENTS 9
#elif defined ADP
#define MAX_MPI_COMPONENTS 14
#elif defined COULOMB
#define MAX_MPI_COMPONENTS 12
#endif /* PAIR */

void broadcast_params()
{
  int   blklens[MAX_MPI_COMPONENTS];
  MPI_Aint displs[MAX_MPI_COMPONENTS];
  MPI_Datatype typen[MAX_MPI_COMPONENTS];
  neigh_t testneigh;
  atom_t testatom;
  int   calclen, size, i, j, each, odd;

  /* Define Structures */
  /* first the easy ones: */
  /* MPI_VEKTOR */
  MPI_Type_contiguous(3, REAL, &MPI_VEKTOR);
  MPI_Type_commit(&MPI_VEKTOR);
  /* MPI_STENS */
  MPI_Type_contiguous(6, REAL, &MPI_STENS);
  MPI_Type_commit(&MPI_STENS);

  /* MPI_NEIGH */
  /* *INDENT-OFF* */
  blklens[0] = 1;         typen[0] = MPI_INT;     /* typ */
  blklens[1] = 1;         typen[1] = MPI_INT;     /* nr */
  blklens[2] = 1;         typen[2] = REAL;        /* r */
  blklens[3] = 1;         typen[3] = MPI_VEKTOR;  /* dist */
  blklens[4] = SLOTS;     typen[4] = MPI_INT;     /* slot */
  blklens[5] = SLOTS;     typen[5] = REAL;        /* shift */
  blklens[6] = SLOTS;     typen[6] = REAL;        /* step */
  blklens[7] = SLOTS;     typen[7] = REAL;        /* col */
  size = 8;
#ifdef ADP
  blklens[8] = 1;         typen[8] = MPI_VEKTOR;  /* rdist */
  blklens[9] = 1;         typen[9] = MPI_STENS;   /* sqrdist */
  blklens[10] = 1;        typen[10] = REAL;        /* u_val */
  blklens[11] = 1;        typen[11] = REAL;       /* u_grad */
  blklens[12] = 1;        typen[12] = REAL;       /* w_val */
  blklens[13] = 1;        typen[13] = REAL;       /* w_grad */
  size += 6;
#endif /* ADP */
#ifdef COULOMB
  blklens[8] = 1;         typen[8] = REAL;        /* r^2 */
  blklens[9] = 1;         typen[9] = REAL;        /* fnval_el */
  blklens[10] = 1;        typen[10] = REAL;       /* grad_el */
  blklens[11] = 1;        typen[11] = REAL;       /* ggrad_el */
  size += 4;
#endif /* COULOMB */
 
 /* *INDENT-ON* */
  MPI_Address(&testneigh.typ, displs);
  MPI_Address(&testneigh.nr, &displs[1]);
  MPI_Address(&testneigh.r, &displs[2]);
  MPI_Address(&testneigh.dist, &displs[3]);
  MPI_Address(testneigh.slot, &displs[4]);
  MPI_Address(testneigh.shift, &displs[5]);
  MPI_Address(testneigh.step, &displs[6]);
  MPI_Address(testneigh.col, &displs[7]);
#ifdef ADP
  MPI_Address(&testneigh.rdist, &displs[8]);
  MPI_Address(&testneigh.sqrdist, &displs[9]);
  MPI_Address(&testneigh.u_val, &displs[10]);
  MPI_Address(&testneigh.u_grad, &displs[11]);
  MPI_Address(&testneigh.w_val, &displs[12]);
  MPI_Address(&testneigh.w_grad, &displs[13]);
#endif
#ifdef COULOMB
  MPI_Address(&testneigh.r2, &displs[8]);
  MPI_Address(&testneigh.fnval_el, &displs[9]);
  MPI_Address(&testneigh.grad_el, &displs[10]);
  MPI_Address(&testneigh.ggrad_el, &displs[11]);
#endif

  for (i = 1; i < size; i++) {
    displs[i] -= displs[0];
  }
  displs[0] = 0;		/* set displacements */
  MPI_Type_struct(size, blklens, displs, typen, &MPI_NEIGH);
  MPI_Type_commit(&MPI_NEIGH);

  /* MPI_ATOM */
  /* *INDENT-OFF* */
  blklens[0] = 1;         typen[0] = MPI_INT;     /* typ */
  blklens[1] = 1;         typen[1] = MPI_INT;     /* n_neigh */
  blklens[2] = 1;         typen[2] = MPI_VEKTOR;  /* pos */
  blklens[3] = 1;         typen[3] = MPI_VEKTOR;  /* force */
  blklens[4] = 1;         typen[4] = REAL;        /* absforce */
  blklens[5] = 1;         typen[5] = MPI_INT;     /* conf */
  size=6;
#if defined EAM || defined ADP
  blklens[6] = 1;         typen[6] = REAL;        /* rho */
  blklens[7] = 1;         typen[7] = REAL;        /* gradF */
  size += 2;
#endif /* EAM || ADP */
#ifdef ADP
  blklens[8] = 1;         typen[8] = MPI_VEKTOR;  /* mu */
  blklens[9] = 1;         typen[9] = MPI_STENS;   /* lambda */
  blklens[10] = 1;        typen[10] = REAL;       /* nu */
  size += 3;
#endif /* ADP */
#ifdef DIPOLE
  blklens[6] = 1;         typen[6] = MPI_VEKTOR;     /* E_stat */
  blklens[7] = 1;         typen[7] =  MPI_VEKTOR;    /* p_sr */
  blklens[8] = 1;         typen[8] =  MPI_VEKTOR;   /* E_ind */
  blklens[9] = 1;         typen[9] =  MPI_VEKTOR;   /* p_ind */
  blklens[10] = 1;        typen[10] =  MPI_VEKTOR;   /* E_old */
  blklens[11] = 1;        typen[11] =  MPI_VEKTOR;   /* E_tot */
  size += 6;
#endif /* DIPOLE */

  /* DO NOT BROADCAST NEIGHBORS !!! DYNAMIC ALLOCATION */

  /* *INDENT-ON* */
  MPI_Address(&testatom.typ, &displs[0]);
  MPI_Address(&testatom.n_neigh, &displs[1]);
  MPI_Address(&testatom.pos, &displs[2]);
  MPI_Address(&testatom.force, &displs[3]);
  MPI_Address(&testatom.absforce, &displs[4]);
  MPI_Address(&testatom.conf, &displs[5]);
#if defined EAM || defined ADP
  MPI_Address(&testatom.rho, &displs[6]);
  MPI_Address(&testatom.gradF, &displs[7]);
#endif /* EAM || ADP */
#ifdef ADP
  MPI_Address(&testatom.mu, &displs[8]);
  MPI_Address(&testatom.lambda, &displs[9]);
  MPI_Address(&testatom.nu, &displs[10]);
#endif /* ADP */
#ifdef DIPOLE
  MPI_Address(&testatom.E_stat, &displs[6]);
  MPI_Address(&testatom.p_sr, &displs[7]);
  MPI_Address(&testatom.E_ind, &displs[8]);
  MPI_Address(&testatom.p_ind, &displs[9]);
  MPI_Address(&testatom.E_old, &displs[10]);
  MPI_Address(&testatom.E_tot, &displs[11]);
#endif

  for (i = 1; i < size; i++) {
    displs[i] -= displs[0];
  }
  displs[0] = 0;		/* set displacements */

  MPI_Type_struct(size, blklens, displs, typen, &MPI_ATOM);
  MPI_Type_commit(&MPI_ATOM);

  /* Distribute fundamental parameters */
  MPI_Bcast(&mdim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ntypes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nconf, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opt, 1, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef COULOMB
  MPI_Bcast(&dp_kappa, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dp_cut, 1, REAL, 0, MPI_COMM_WORLD);
#endif
#ifdef DIPOLE
  MPI_Bcast(&dp_tol, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dp_mix, 1, REAL, 0, MPI_COMM_WORLD);
#endif
  if (myid > 0) {
    inconf = (int *)malloc(nconf * sizeof(int));
    cnfstart = (int *)malloc(nconf * sizeof(int));
    force_0 = (real *)malloc(mdim * sizeof(real));
    conf_weight = (real *)malloc(nconf * sizeof(real));
  }
  MPI_Bcast(inconf, nconf, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(cnfstart, nconf, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(force_0, mdim, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(conf_weight, nconf, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&maxneigh, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Broadcast weights... */
  MPI_Bcast(&eweight, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sweight, 1, REAL, 0, MPI_COMM_WORLD);

  /* Broadcast the potential... */
  MPI_Bcast(&format, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&calc_pot.len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&calc_pot.ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);
  size = calc_pot.ncols;
  calclen = calc_pot.len;
  if (myid > 0) {
    calc_pot.begin = (real *)malloc(size * sizeof(real));
    calc_pot.end = (real *)malloc(size * sizeof(real));
    calc_pot.step = (real *)malloc(size * sizeof(real));
    calc_pot.invstep = (real *)malloc(size * sizeof(real));
    calc_pot.first = (int *)malloc(size * sizeof(int));
    calc_pot.last = (int *)malloc(size * sizeof(int));
    calc_pot.table = (real *)malloc(calclen * sizeof(real));
    calc_pot.xcoord = (real *)malloc(calclen * sizeof(real));
    calc_pot.d2tab = (real *)malloc(calclen * sizeof(real));
  }
  MPI_Bcast(calc_pot.begin, size, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.end, size, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.step, size, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.invstep, size, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.first, size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.last, size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.table, calclen, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.d2tab, calclen, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.xcoord, calclen, REAL, 0, MPI_COMM_WORLD);

#ifdef APOT
  MPI_Bcast(&do_smooth, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&enable_cp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opt_pot.len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&apot_table.number, 1, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef COULOMB
  MPI_Bcast(&apot_table.total_ne_par, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  if (enable_cp) {
    if (myid > 0) {
      na_type = (int **)malloc((nconf + 1) * sizeof(int *));
      for (i = 0; i < (nconf + 1); i++)
	na_type[i] = (int *)malloc(ntypes * sizeof(int));
    }
    for (i = 0; i < (nconf + 1); i++)
      MPI_Bcast(na_type[i], ntypes, MPI_INT, 0, MPI_COMM_WORLD);
  }
  if (myid > 0) {
    calc_list = (real *)malloc(opt_pot.len * sizeof(real));
    apot_table.n_par = (int *)malloc(apot_table.number * sizeof(int));
    apot_table.end = (real *)malloc(apot_table.number * sizeof(real));
    apot_table.begin = (real *)malloc(apot_table.number * sizeof(real));
#ifdef COULOMB
    apot_table.ratio = (real *)malloc(2 * sizeof(real));
#endif
    smooth_pot = (int *)malloc(apot_table.number * sizeof(int));
    invar_pot = (int *)malloc(apot_table.number * sizeof(int));
    rcut = (real *)malloc(ntypes * ntypes * sizeof(real));
    rmin = (real *)malloc(ntypes * ntypes * sizeof(real));
    apot_table.fvalue =
      (fvalue_pointer *) malloc(apot_table.number * sizeof(fvalue_pointer));
    opt_pot.table = (real *)malloc(opt_pot.len * sizeof(real));
  }
  MPI_Bcast(smooth_pot, apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(invar_pot, apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_list, opt_pot.len, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(apot_table.n_par, apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(rcut, ntypes * ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(rmin, ntypes * ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(apot_table.fvalue, apot_table.number, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(apot_table.end, apot_table.number, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(apot_table.begin, apot_table.number, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cp_start, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&have_globals, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&global_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&apot_table.globals, 1, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef COULOMB
  MPI_Bcast(&apot_table.last_charge, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(apot_table.ratio, 2, REAL, 0, MPI_COMM_WORLD);
#endif
  if (have_globals) {
    if (myid > 0) {
      apot_table.n_glob = (int *)malloc(apot_table.globals * sizeof(int));
      apot_table.global_idx =
	(int ***)malloc(apot_table.globals * sizeof(int **));
      opt_pot.first = (int *)malloc(apot_table.number * sizeof(int));
    }
    MPI_Bcast(apot_table.n_glob, apot_table.globals, MPI_INT, 0,
      MPI_COMM_WORLD);
    MPI_Bcast(opt_pot.first, apot_table.number, MPI_INT, 0, MPI_COMM_WORLD);
    if (myid > 0) {
      for (i = 0; i < apot_table.globals; i++)
	apot_table.global_idx[i] =
	  (int **)malloc(apot_table.n_glob[i] * sizeof(int *));
      for (i = 0; i < apot_table.globals; i++)
	for (j = 0; j < apot_table.n_glob[i]; j++)
	  apot_table.global_idx[i][j] = (int *)malloc(2 * sizeof(int));
    }
    for (i = 0; i < apot_table.globals; i++)
      for (j = 0; j < apot_table.n_glob[i]; j++)
	MPI_Bcast(apot_table.global_idx[i][j], 2, MPI_INT, 0, MPI_COMM_WORLD);
  }
#endif

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
  }
  MPI_Scatter(atom_len, 1, MPI_INT, &myatoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(atom_dist, 1, MPI_INT, &firstatom, 1, MPI_INT, 0,
    MPI_COMM_WORLD);
  MPI_Scatter(conf_len, 1, MPI_INT, &myconf, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(conf_dist, 1, MPI_INT, &firstconf, 1, MPI_INT, 0,
    MPI_COMM_WORLD);
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
  conf_vol = (real *)malloc(myconf * sizeof(real));
  conf_uf = (int *)malloc(myconf * sizeof(real));
  conf_us = (int *)malloc(myconf * sizeof(real));
  MPI_Scatterv(volumen, conf_len, conf_dist, REAL,
    conf_vol, myconf, REAL, 0, MPI_COMM_WORLD);
  MPI_Scatterv(useforce, conf_len, conf_dist, MPI_INT,
    conf_uf, myconf, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(usestress, conf_len, conf_dist, MPI_INT,
    conf_us, myconf, MPI_INT, 0, MPI_COMM_WORLD);
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
      neighs = atoms[i].n_neigh;
    MPI_Bcast(&neighs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (i >= firstatom && i < (firstatom + myatoms)) {
      atom->neigh = (neigh_t *)malloc(neighs * sizeof(neigh_t));
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
  MPI_Bcast(calc_pot.begin + firstcol, ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.end + firstcol, ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.step + firstcol, ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.invstep + firstcol, ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(calc_pot.first + firstcol, ntypes, MPI_INT, 0, MPI_COMM_WORLD);
  /* bcast table values of transfer fn. and embedding energy */
  firstval = calc_pot.first[paircol];
  nvals = calc_pot.len - firstval;
  MPI_Bcast(calc_pot.table + firstval, nvals, REAL, 0, MPI_COMM_WORLD);
}

#endif /* !APOT */
#endif /* MPI */
