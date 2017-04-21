/****************************************************************
 *
 * mpi_utils.c: Contains utilities to be used with MPI
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

#include "config.h"
#include "memory.h"
#include "mpi_utils.h"
#include "utils.h"

#define CHECK_RETURN(a)                                            \
  do {                                                             \
    int r = a;                                                     \
    if (r != MPI_SUCCESS) {                                        \
      printf("Error in %s on line %d: %d", __FILE__, __LINE__, r); \
      return POTFIT_ERROR;                                         \
    }                                                              \
  } while (0);

#if defined(MPI)
int create_custom_datatypes();
int broadcast_basic_data();
int broadcast_calcpot_table();
int broadcast_apot_table();
int broadcast_configurations();
int broadcast_atoms();
int broadcast_neighbors();
int broadcast_angles();
#endif  // MPI

/****************************************************************
  set up MPI communication
****************************************************************/

int initialize_mpi(int* argc, char*** argv)
{
#if defined(MPI)
  // initialize the MPI communication
  int rval = MPI_Init(argc, argv);

  if (rval != MPI_SUCCESS) {
    printf("Error initializing MPI communication! (Error: %d)\n", rval);
    return POTFIT_ERROR;
  }

  rval = MPI_Comm_size(MPI_COMM_WORLD, &g_mpi.num_cpus);

  if (rval != MPI_SUCCESS) {
    printf("Error getting MPI communicator size! (Error: %d)\n", rval);
    return POTFIT_ERROR;
  }

  rval = MPI_Comm_rank(MPI_COMM_WORLD, &g_mpi.myid);

  if (rval != MPI_SUCCESS) {
    printf("Error getting MPI communicator rank! (Error: %d)\n", rval);
    return POTFIT_ERROR;
  }
#endif  // MPI

  if (g_mpi.myid == 0) {
    printf("This is %s compiled on %s, %s.\n", POTFIT_VERSION, __DATE__,
           __TIME__);
#if defined(MPI)
    printf("Starting up MPI with %d processes.\n", g_mpi.num_cpus);
#endif  // MPI
    printf("\n");
  }

  return POTFIT_SUCCESS;
}

/****************************************************************
  shutdown_mpi
****************************************************************/

void shutdown_mpi(void)
{
#if defined(MPI)
  if (g_mpi.init_done == 0) {
    int abort = -1;
    MPI_Bcast(&abort, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  MPI_Finalize(); /* Shutdown */
#endif            // MPI
}

/****************************************************************
  broadcast_params_mpi: Broadcast parameters etc to other nodes
****************************************************************/

int broadcast_params_mpi()
{
#if defined(MPI)
  if (g_mpi.myid == 0) {
    printf("Broadcasting data to MPI workers ... ");
    fflush(stdout);
  }

  // non-root processes will wait here while root is reading input files
  CHECK_RETURN(MPI_Bcast(&g_mpi.init_done, 1, MPI_INT, 0, MPI_COMM_WORLD));

  // shutdown_mpi will broadcast value of -1 to signal failure
  if (g_mpi.init_done == -1)
    return POTFIT_ERROR_MPI_CLEAN_EXIT;

  CHECK_RETURN(create_custom_datatypes());
  CHECK_RETURN(broadcast_basic_data());
  CHECK_RETURN(broadcast_calcpot_table());
  CHECK_RETURN(broadcast_apot_table());
  CHECK_RETURN(broadcast_configurations());
  CHECK_RETURN(broadcast_atoms());
  CHECK_RETURN(broadcast_neighbors());
  CHECK_RETURN(broadcast_angles());

  if (g_mpi.myid == 0) {
    printf("done\n");
    fflush(stdout);
  }
#else
  // Identify subset of atoms/volumes belonging to individual
  // process with complete set of atoms/volumes
  g_config.conf_atoms = g_config.atoms;
  g_config.conf_vol = g_config.volume;
  g_config.conf_uf = g_config.useforce;
#if defined(STRESS)
  g_config.conf_us = g_config.usestress;
#endif  // STRESS
#endif  // MPI

  return POTFIT_SUCCESS;
}

#if !defined(APOT)

/****************************************************************
    potsync: Broadcast parameters etc to other nodes
****************************************************************/

void potsync()
{
#if defined(MPI)
  int firstcol = g_calc.paircol + g_param.ntypes;

  /* Memory is allocated - just bcast that changed potential... */
  /* bcast begin/end/step/invstep of embedding energy  */
  MPI_Bcast(g_pot.calc_pot.begin + firstcol, g_param.ntypes, MPI_DOUBLE, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(g_pot.calc_pot.end + firstcol, g_param.ntypes, MPI_DOUBLE, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(g_pot.calc_pot.step + firstcol, g_param.ntypes, MPI_DOUBLE, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(g_pot.calc_pot.invstep + firstcol, g_param.ntypes, MPI_DOUBLE, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(g_pot.calc_pot.first + firstcol, g_param.ntypes, MPI_INT, 0,
            MPI_COMM_WORLD);

  /* bcast table values of transfer fn. and embedding energy */
  int firstval = g_pot.calc_pot.first[g_calc.paircol];
  int nvals = g_pot.calc_pot.len - firstval;
  MPI_Bcast(g_pot.calc_pot.table + firstval, nvals, MPI_DOUBLE, 0,
            MPI_COMM_WORLD);
#endif  // MPI
}

#endif  // !APOT

#if defined(MPI)

/****************************************************************
    create_custom_datatypes
****************************************************************/

int create_custom_datatypes()
{
  int size_a = 0;
  int size_b = 0;

#define MAX_MPI_COMPONENTS 20
  int data_len[MAX_MPI_COMPONENTS];
  MPI_Datatype data_type[MAX_MPI_COMPONENTS];
  MPI_Aint data_size[MAX_MPI_COMPONENTS];
#undef MAX_MPI_COMPONENTS

  // define MPI structures
  // first the easy ones:

  // MPI_VECTOR
  CHECK_RETURN(MPI_Type_contiguous(3, MPI_DOUBLE, &g_mpi.MPI_VECTOR));
  CHECK_RETURN(MPI_Type_commit(&g_mpi.MPI_VECTOR));

  // MPI_STENS
  CHECK_RETURN(MPI_Type_contiguous(6, MPI_DOUBLE, &g_mpi.MPI_STENS));
  CHECK_RETURN(MPI_Type_commit(&g_mpi.MPI_STENS));

  // clang-format off

  // MPI_NEIGH
  data_len[size_a] = 1;       data_type[size_a++] = MPI_INT;            // type
  data_len[size_a] = 1;       data_type[size_a++] = MPI_INT;            // nr
  data_len[size_a] = 1;       data_type[size_a++] = MPI_DOUBLE;         // r
  data_len[size_a] = 1;       data_type[size_a++] = MPI_DOUBLE;         // r2
  data_len[size_a] = 1;       data_type[size_a++] = MPI_DOUBLE;         // inv_r
  data_len[size_a] = 1;       data_type[size_a++] = g_mpi.MPI_VECTOR;   // dist
  data_len[size_a] = 1;       data_type[size_a++] = g_mpi.MPI_VECTOR;   // dist_r
  data_len[size_a] = SLOTS;   data_type[size_a++] = MPI_INT;            // slot
  data_len[size_a] = SLOTS;   data_type[size_a++] = MPI_DOUBLE;         // shift
  data_len[size_a] = SLOTS;   data_type[size_a++] = MPI_DOUBLE;         // step
  data_len[size_a] = SLOTS;   data_type[size_a++] = MPI_INT;            // col
#if defined(ADP)
  data_len[size_a] = 1;       data_type[size_a++] = g_mpi.MPI_STENS;    // sqrdist
  data_len[size_a] = 1;       data_type[size_a++] = MPI_DOUBLE;         // u_val
  data_len[size_a] = 1;       data_type[size_a++] = MPI_DOUBLE;         // u_grad
  data_len[size_a] = 1;       data_type[size_a++] = MPI_DOUBLE;         // w_val
  data_len[size_a] = 1;       data_type[size_a++] = MPI_DOUBLE;         // w_grad
#endif // ADP
#if defined(COULOMB)
  data_len[size_a] = 1;       data_type[size_a++] = MPI_DOUBLE;         // fnval_el
  data_len[size_a] = 1;       data_type[size_a++] = MPI_DOUBLE;         // grad_el
  data_len[size_a] = 1;       data_type[size_a++] = MPI_DOUBLE;         // ggrad_el
#endif // COULOMB
#if defined(THREEBODY)
  data_len[size_a] = 1;       data_type[size_a++] = MPI_DOUBLE;         // f
  data_len[size_a] = 1;       data_type[size_a++] = MPI_DOUBLE;         // df
  data_len[size_a] = 1;       data_type[size_a++] = MPI_INT;            // ijk_start
#endif // THREEBODY
#if defined(MEAM)
  data_len[size_a] = 1;       data_type[size_a++] = MPI_DOUBLE;         // drho
#endif // MEAM
#if defined(TERSOFF)
  data_len[size_a] = 1;       data_type[size_a++] = g_mpi.MPI_VECTOR;   // dzeta
#endif // TERSOFF

  neigh_t neigh;
  CHECK_RETURN(MPI_Get_address(&neigh.type,      &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.nr,        &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.r,         &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.r2,        &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.inv_r,     &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.dist,      &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.dist_r,    &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(neigh.slot,       &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(neigh.shift,      &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(neigh.step,       &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(neigh.col,        &data_size[size_b++]));
#if defined(ADP)
  CHECK_RETURN(MPI_Get_address(&neigh.sqrdist,   &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.u_val,     &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.u_grad,    &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.w_val,     &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.w_grad,    &data_size[size_b++]));
#endif // ADP
#if defined(COULOMB)
  CHECK_RETURN(MPI_Get_address(&neigh.fnval_el,  &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.grad_el,   &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.ggrad_el,  &data_size[size_b++]));
#endif // COULOMB
#if defined(THREEBODY)
  CHECK_RETURN(MPI_Get_address(&neigh.f,         &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.df,        &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&neigh.ijk_start, &data_size[size_b++]));
#endif // THREEBODY
#if defined(MEAM)
  CHECK_RETURN(MPI_Get_address(&neigh.drho,      &data_size[size_b++]));
#endif // MEAM
#if defined(TERSOFF)
  CHECK_RETURN(MPI_Get_address(&neigh.dzeta,     &data_size[size_b++]));
#endif // TERSOFF

  // clang-format on

  if (size_a != size_b) {
    printf("Sizes of item array and address array differ in %s:%d", __FILE__,
           __LINE__);
    MPI_Finalize();
    exit(EXIT_SUCCESS);
  }

  // calculate displacements from addresses
  for (int i = 1; i < size_a; ++i)
    data_size[i] -= data_size[0];
  data_size[0] = 0;

  CHECK_RETURN(MPI_Type_create_struct(size_a, data_len, data_size, data_type,
                                      &g_mpi.MPI_NEIGH));
  CHECK_RETURN(MPI_Type_commit(&g_mpi.MPI_NEIGH));

#if defined(THREEBODY)
  // clang-format off

  // MPI_ANGL
  size_a = 0;
  size_b = 0;
  angle_t angl;

  data_len[size_a] = 1;         data_type[size_a++] = MPI_DOUBLE;   // cos
#if defined(ANG)
  data_len[size_a] = 1;         data_type[size_a++] = MPI_DOUBLE;   // theta
#endif
#if defined(MEAM) || defined(ANG)
  data_len[size_a] = 1;         data_type[size_a++] = MPI_INT;      // slot
  data_len[size_a] = 1;         data_type[size_a++] = MPI_DOUBLE;   // shift
  data_len[size_a] = 1;         data_type[size_a++] = MPI_DOUBLE;   // step
  data_len[size_a] = 1;         data_type[size_a++] = MPI_DOUBLE;   // g
  data_len[size_a] = 1;         data_type[size_a++] = MPI_DOUBLE;   // dg
#endif // MEAM || ANG

  CHECK_RETURN(MPI_Get_address(&angl.cos,    &data_size[size_b++]));
#if defined(ANG)
  CHECK_RETURN(MPI_Get_address(&angl.theta,    &data_size[size_b++]));
#endif
#if defined(MEAM) || defined(ANG)
  CHECK_RETURN(MPI_Get_address(&angl.slot,   &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&angl.shift,  &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&angl.step,   &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&angl.g,      &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&angl.dg,     &data_size[size_b++]));
#endif // MEAM

  // clang-format on

  if (size_a != size_b) {
    printf("Sizes of item array and address array differ in %s:%d", __FILE__,
           __LINE__);
    MPI_Finalize();
    exit(EXIT_SUCCESS);
  }

  // calculate displacements from addresses
  for (int i = 1; i < size_a; ++i)
    data_size[i] -= data_size[0];
  data_size[0] = 0;

  CHECK_RETURN(MPI_Type_create_struct(size_a, data_len, data_size, data_type,
                                      &g_mpi.MPI_ANGL));
  CHECK_RETURN(MPI_Type_commit(&g_mpi.MPI_ANGL));
#endif  // THREEBODY

  // clang-format off

  // MPI_ATOM
  size_a = 0;
  size_b = 0;
  atom_t atom;

  data_len[size_a] = 1;     data_type[size_a++] = MPI_INT;          // type
  data_len[size_a] = 1;     data_type[size_a++] = MPI_INT;          // num_neigh
  data_len[size_a] = 1;     data_type[size_a++] = g_mpi.MPI_VECTOR; // pos
  data_len[size_a] = 1;     data_type[size_a++] = g_mpi.MPI_VECTOR; // force
  data_len[size_a] = 1;     data_type[size_a++] = MPI_DOUBLE;       // absforce
  data_len[size_a] = 1;     data_type[size_a++] = MPI_INT;          // conf
#if defined(CONTRIB)
  data_len[size_a] = 1;     data_type[size_a++] = MPI_INT;          // contrib
#endif // CONTRIB
#if defined(EAM) || defined(ADP) || defined(MEAM)
  data_len[size_a] = 1;     data_type[size_a++] = MPI_DOUBLE;       // rho
  data_len[size_a] = 1;     data_type[size_a++] = MPI_DOUBLE;       // gradF
#endif // EAM || ADP || MEAM
#if defined(ADP)
  data_len[size_a] = 1;     data_type[size_a++] = g_mpi.MPI_VECTOR; // mu
  data_len[size_a] = 1;     data_type[size_a++] = g_mpi.MPI_STENS;  // lambda
  data_len[size_a] = 1;     data_type[size_a++] = MPI_DOUBLE;       // nu
#endif // ADP
#if defined(DIPOLE)
  data_len[size_a] = 1;     data_type[size_a++] = g_mpi.MPI_VECTOR; // E_stat
  data_len[size_a] = 1;     data_type[size_a++] = g_mpi.MPI_VECTOR; // p_sr
  data_len[size_a] = 1;     data_type[size_a++] = g_mpi.MPI_VECTOR; // E_ind
  data_len[size_a] = 1;     data_type[size_a++] = g_mpi.MPI_VECTOR; // p_ind
  data_len[size_a] = 1;     data_type[size_a++] = g_mpi.MPI_VECTOR; // E_old
  data_len[size_a] = 1;     data_type[size_a++] = g_mpi.MPI_VECTOR; // E_tot
#endif // DIPOLE
#if defined(THREEBODY)
  data_len[size_a] = 1;     data_type[size_a++] = MPI_INT;          // num_angl
#if defined(MEAM)
  data_len[size_a] = 1;     data_type[size_a++] = MPI_DOUBLE;       // rho_eam
#endif // MEAM
#endif // THREEBODY

  // DO NOT BROADCAST NEIGHBORS !!! DYNAMIC ALLOCATION
  // DO NOT BROADCAST ANGLES !!! DYNAMIC ALLOCATION

  CHECK_RETURN(MPI_Get_address(&atom.type,       &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&atom.num_neigh,  &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&atom.pos,        &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&atom.force,      &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&atom.absforce,   &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&atom.conf,       &data_size[size_b++]));
#if defined(CONTRIB)
  CHECK_RETURN(MPI_Get_address(&atom.contrib,    &data_size[size_b++]));
#endif // CONTRIB
#if defined(EAM) || defined(ADP) || defined(MEAM)
  CHECK_RETURN(MPI_Get_address(&atom.rho,        &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&atom.gradF,      &data_size[size_b++]));
#endif // EAM || ADP
#if defined(ADP)
  CHECK_RETURN(MPI_Get_address(&atom.mu,         &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&atom.lambda,     &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&atom.nu,         &data_size[size_b++]));
#endif // ADP
#if defined(DIPOLE)
  CHECK_RETURN(MPI_Get_address(&atom.E_stat,     &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&atom.p_sr,       &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&atom.E_ind,      &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&atom.p_ind,      &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&atom.E_old,      &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&atom.E_tot,      &data_size[size_b++]));
#endif // DIPOLE
#if defined(THREEBODY)
  CHECK_RETURN(MPI_Get_address(&atom.num_angles, &data_size[size_b++]));
#if defined(MEAM)
  CHECK_RETURN(MPI_Get_address(&atom.rho_eam,    &data_size[size_b++]));
#endif // MEAM
#endif // THREEBODY

  // clang-format on

  if (size_a != size_b) {
    printf("Sizes of item array and address array differ in %s:%d", __FILE__,
           __LINE__);
    MPI_Finalize();
    exit(EXIT_SUCCESS);
  }

  // calculate displacements from addresses
  for (int i = 1; i < size_a; ++i)
    data_size[i] -= data_size[0];
  data_size[0] = 0;

  CHECK_RETURN(MPI_Type_create_struct(size_a, data_len, data_size, data_type,
                                      &g_mpi.MPI_ATOM));
  CHECK_RETURN(MPI_Type_commit(&g_mpi.MPI_ATOM));

  return MPI_SUCCESS;
}

/****************************************************************
    broadcast_basic_data
****************************************************************/

int broadcast_basic_data()
{
  // Distribute fundamental parameters
  CHECK_RETURN(MPI_Bcast(&g_calc.mdim, 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(&g_param.ntypes, 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(&g_config.natoms, 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(&g_config.nconf, 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(&g_calc.paircol, 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(&g_param.opt, 1, MPI_INT, 0, MPI_COMM_WORLD));

  // allocate and broadcast config metadata
  if (g_mpi.myid > 0) {
    g_config.inconf = (int*)Malloc(g_config.nconf * sizeof(int));
    g_config.cnfstart = (int*)Malloc(g_config.nconf * sizeof(int));
    g_config.force_0 = (double*)Malloc(g_calc.mdim * sizeof(double));
    g_config.conf_weight = (double*)Malloc(g_config.nconf * sizeof(double));
  }
  CHECK_RETURN(
      MPI_Bcast(g_config.inconf, g_config.nconf, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(
      MPI_Bcast(g_config.cnfstart, g_config.nconf, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(
      MPI_Bcast(g_config.force_0, g_calc.mdim, MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(g_config.conf_weight, g_config.nconf, MPI_DOUBLE, 0,
                         MPI_COMM_WORLD));

  // Broadcast weights...
  CHECK_RETURN(MPI_Bcast(&g_param.eweight, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(&g_param.sweight, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));

  // Broadcast the potential...
  CHECK_RETURN(MPI_Bcast(&g_pot.format_type, 1, MPI_INT, 0, MPI_COMM_WORLD));

#if defined(COULOMB)
  CHECK_RETURN(MPI_Bcast(&g_config.dp_cut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
#endif  // COULOMB
#if defined(DIPOLE)
  CHECK_RETURN(MPI_Bcast(&g_config.dp_tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(&g_config.dp_mix, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
#endif  // DIPOLE

  return MPI_SUCCESS;
}

/****************************************************************
    broadcast_calcpot_table
****************************************************************/

int broadcast_calcpot_table()
{
  CHECK_RETURN(MPI_Bcast(&g_pot.calc_pot.len, 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(&g_pot.calc_pot.ncols, 1, MPI_INT, 0, MPI_COMM_WORLD));

  int ncols = g_pot.calc_pot.ncols;
  int calclen = g_pot.calc_pot.len;

  if (g_mpi.myid > 0) {
    g_pot.calc_pot.begin = (double*)Malloc(ncols * sizeof(double));
    g_pot.calc_pot.end = (double*)Malloc(ncols * sizeof(double));
    g_pot.calc_pot.step = (double*)Malloc(ncols * sizeof(double));
    g_pot.calc_pot.invstep = (double*)Malloc(ncols * sizeof(double));
    g_pot.calc_pot.first = (int*)Malloc(ncols * sizeof(int));
    g_pot.calc_pot.last = (int*)Malloc(ncols * sizeof(int));
    g_pot.calc_pot.table = (double*)Malloc(calclen * sizeof(double));
    g_pot.calc_pot.xcoord = (double*)Malloc(calclen * sizeof(double));
    g_pot.calc_pot.d2tab = (double*)Malloc(calclen * sizeof(double));
  }

  CHECK_RETURN(
      MPI_Bcast(g_pot.calc_pot.begin, ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(
      MPI_Bcast(g_pot.calc_pot.end, ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(
      MPI_Bcast(g_pot.calc_pot.step, ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(
      MPI_Bcast(g_pot.calc_pot.invstep, ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(
      MPI_Bcast(g_pot.calc_pot.first, ncols, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(
      MPI_Bcast(g_pot.calc_pot.last, ncols, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(
      MPI_Bcast(g_pot.calc_pot.table, calclen, MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(
      MPI_Bcast(g_pot.calc_pot.d2tab, calclen, MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(
      MPI_Bcast(g_pot.calc_pot.xcoord, calclen, MPI_DOUBLE, 0, MPI_COMM_WORLD));

  return MPI_SUCCESS;
}

/****************************************************************
    broadcast_apot_table
****************************************************************/

int broadcast_apot_table()
{
#if defined(APOT)
  CHECK_RETURN(MPI_Bcast(&g_param.enable_cp, 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(&g_pot.opt_pot.len, 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(
      MPI_Bcast(&g_pot.apot_table.number, 1, MPI_INT, 0, MPI_COMM_WORLD));
#if defined(COULOMB)
  CHECK_RETURN(
      MPI_Bcast(&g_pot.apot_table.total_ne_par, 1, MPI_INT, 0, MPI_COMM_WORLD));
#endif  // COULOMB

  if (g_param.enable_cp) {
    if (g_mpi.myid > 0) {
      g_config.na_type = (int**)Malloc((g_config.nconf + 1) * sizeof(int*));
      for (int i = 0; i < (g_config.nconf + 1); i++)
        g_config.na_type[i] = (int*)Malloc(g_param.ntypes * sizeof(int));
    }
    for (int i = 0; i < (g_config.nconf + 1); i++)
      CHECK_RETURN(MPI_Bcast(g_config.na_type[i], g_param.ntypes, MPI_INT, 0,
                             MPI_COMM_WORLD));
  }

  if (g_mpi.myid > 0) {
    g_pot.calc_list = (double*)Malloc(g_pot.opt_pot.len * sizeof(double));
    g_pot.apot_table.n_par =
        (int*)Malloc(g_pot.apot_table.number * sizeof(int));
    g_pot.apot_table.begin =
        (double*)Malloc(g_pot.apot_table.number * sizeof(double));
    g_pot.apot_table.end =
        (double*)Malloc(g_pot.apot_table.number * sizeof(double));
    g_pot.apot_table.idxpot =
        (int*)Malloc(g_pot.apot_table.number * sizeof(int));
#if defined(COULOMB)
    g_pot.apot_table.ratio = (double*)Malloc(g_param.ntypes * sizeof(double));
#endif  // COULOMB
    g_pot.smooth_pot = (int*)Malloc(g_pot.apot_table.number * sizeof(int));
    g_pot.invar_pot = (int*)Malloc(g_pot.apot_table.number * sizeof(int));
    g_config.rcut =
        (double*)Malloc(g_param.ntypes * g_param.ntypes * sizeof(double));
    g_config.rmin =
        (double*)Malloc(g_param.ntypes * g_param.ntypes * sizeof(double));
    g_pot.apot_table.fvalue = (fvalue_pointer*)Malloc(g_pot.apot_table.number *
                                                      sizeof(fvalue_pointer));
    g_pot.opt_pot.table = (double*)Malloc(g_pot.opt_pot.len * sizeof(double));
    g_pot.opt_pot.first = (int*)Malloc(g_pot.apot_table.number * sizeof(int));
  }

  CHECK_RETURN(MPI_Bcast(g_pot.smooth_pot, g_pot.apot_table.number, MPI_INT, 0,
                         MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(g_pot.invar_pot, g_pot.apot_table.number, MPI_INT, 0,
                         MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(g_pot.calc_list, g_pot.opt_pot.len, MPI_DOUBLE, 0,
                         MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(g_pot.apot_table.n_par, g_pot.apot_table.number,
                         MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(g_config.rcut, g_param.ntypes * g_param.ntypes,
                         MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(g_config.rmin, g_param.ntypes * g_param.ntypes,
                         MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(g_pot.apot_table.fvalue, g_pot.apot_table.number,
                         MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(g_pot.apot_table.end, g_pot.apot_table.number,
                         MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(g_pot.apot_table.begin, g_pot.apot_table.number,
                         MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(g_pot.apot_table.idxpot, g_pot.apot_table.number,
                         MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(&g_pot.cp_start, 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(&g_pot.have_globals, 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(&g_pot.global_idx, 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(
      MPI_Bcast(&g_pot.apot_table.globals, 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(g_pot.opt_pot.first, g_pot.apot_table.number, MPI_INT,
                         0, MPI_COMM_WORLD));
#if defined(COULOMB)
  CHECK_RETURN(MPI_Bcast(&g_pot.apot_table.last_charge, 1, MPI_DOUBLE, 0,
                         MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Bcast(g_pot.apot_table.ratio, g_param.ntypes, MPI_DOUBLE, 0,
                         MPI_COMM_WORLD));
#endif  // COULOMB

  if (g_pot.have_globals) {
    if (g_mpi.myid > 0) {
      g_pot.apot_table.n_glob =
          (int*)Malloc(g_pot.apot_table.globals * sizeof(int));
      g_pot.apot_table.global_idx =
          (int***)Malloc(g_pot.apot_table.globals * sizeof(int**));
    }
    CHECK_RETURN(MPI_Bcast(g_pot.apot_table.n_glob, g_pot.apot_table.globals,
                           MPI_INT, 0, MPI_COMM_WORLD));
    if (g_mpi.myid > 0) {
      for (int i = 0; i < g_pot.apot_table.globals; i++)
        g_pot.apot_table.global_idx[i] =
            (int**)Malloc(g_pot.apot_table.n_glob[i] * sizeof(int*));
      for (int i = 0; i < g_pot.apot_table.globals; i++)
        for (int j = 0; j < g_pot.apot_table.n_glob[i]; j++)
          g_pot.apot_table.global_idx[i][j] = (int*)Malloc(2 * sizeof(int));
    }
    for (int i = 0; i < g_pot.apot_table.globals; i++)
      for (int j = 0; j < g_pot.apot_table.n_glob[i]; j++)
        CHECK_RETURN(MPI_Bcast(g_pot.apot_table.global_idx[i][j], 2, MPI_INT, 0,
                               MPI_COMM_WORLD));
  }
#endif  // APOT

  return MPI_SUCCESS;
}

/****************************************************************
    broadcast_configurations
****************************************************************/

int broadcast_configurations()
{
  // Each node: nconf/num_cpus configurations.
  // Last nconf%num_cpus nodes: 1 additional config

  if (g_mpi.myid == 0) {
    int each = (g_config.nconf / g_mpi.num_cpus);
    int odd = (g_config.nconf % g_mpi.num_cpus) - g_mpi.num_cpus;

    g_mpi.atom_len = (int*)Malloc(g_mpi.num_cpus * sizeof(int));
    g_mpi.atom_dist = (int*)Malloc(g_mpi.num_cpus * sizeof(int));
    g_mpi.conf_len = (int*)Malloc(g_mpi.num_cpus * sizeof(int));
    g_mpi.conf_dist = (int*)Malloc(g_mpi.num_cpus * sizeof(int));
    for (int i = 0; i < g_mpi.num_cpus; i++)
      g_mpi.conf_dist[i] = i * each + (((i + odd) > 0) ? (i + odd) : 0);
    for (int i = 0; i < g_mpi.num_cpus - 1; i++)
      g_mpi.conf_len[i] = g_mpi.conf_dist[i + 1] - g_mpi.conf_dist[i];
    g_mpi.conf_len[g_mpi.num_cpus - 1] =
        g_config.nconf - g_mpi.conf_dist[g_mpi.num_cpus - 1];
    for (int i = 0; i < g_mpi.num_cpus; i++)
      g_mpi.atom_dist[i] =
          g_config.cnfstart[i * each + (((i + odd) > 0) ? (i + odd) : 0)];
    for (int i = 0; i < g_mpi.num_cpus - 1; i++)
      g_mpi.atom_len[i] = g_mpi.atom_dist[i + 1] - g_mpi.atom_dist[i];
    g_mpi.atom_len[g_mpi.num_cpus - 1] =
        g_config.natoms - g_mpi.atom_dist[g_mpi.num_cpus - 1];
  }

  CHECK_RETURN(MPI_Scatter(g_mpi.atom_len, 1, MPI_INT, &g_mpi.myatoms, 1,
                           MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Scatter(g_mpi.atom_dist, 1, MPI_INT, &g_mpi.firstatom, 1,
                           MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Scatter(g_mpi.conf_len, 1, MPI_INT, &g_mpi.myconf, 1,
                           MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Scatter(g_mpi.conf_dist, 1, MPI_INT, &g_mpi.firstconf, 1,
                           MPI_INT, 0, MPI_COMM_WORLD));

  g_config.conf_vol = (double*)Malloc(g_mpi.myconf * sizeof(double));
  g_config.conf_uf = (int*)Malloc(g_mpi.myconf * sizeof(int));

  CHECK_RETURN(MPI_Scatterv(g_config.volume, g_mpi.conf_len, g_mpi.conf_dist,
                            MPI_DOUBLE, g_config.conf_vol, g_mpi.myconf,
                            MPI_DOUBLE, 0, MPI_COMM_WORLD));
  CHECK_RETURN(MPI_Scatterv(g_config.useforce, g_mpi.conf_len, g_mpi.conf_dist,
                            MPI_INT, g_config.conf_uf, g_mpi.myconf, MPI_INT, 0,
                            MPI_COMM_WORLD));

#if defined(STRESS)
  g_config.conf_us = (int*)Malloc(g_mpi.myconf * sizeof(double));

  CHECK_RETURN(MPI_Scatterv(g_config.usestress, g_mpi.conf_len, g_mpi.conf_dist,
                            MPI_INT, g_config.conf_us, g_mpi.myconf, MPI_INT, 0,
                            MPI_COMM_WORLD));
#endif  // STRESS

  return MPI_SUCCESS;
}

/****************************************************************
    scatter atom table
****************************************************************/

int broadcast_atoms()
{
  atom_t atom;

  g_config.conf_atoms = (atom_t*)Malloc(g_mpi.myatoms * sizeof(atom_t));

  for (int i = 0; i < g_config.natoms; i++) {
    if (g_mpi.myid == 0)
      atom = g_config.atoms[i];
    CHECK_RETURN(MPI_Bcast(&atom, 1, g_mpi.MPI_ATOM, 0, MPI_COMM_WORLD));
    if (i >= g_mpi.firstatom && i < (g_mpi.firstatom + g_mpi.myatoms))
      g_config.conf_atoms[i - g_mpi.firstatom] = atom;
  }

  return MPI_SUCCESS;
}

/****************************************************************
    scatter dynamic neighbor table
****************************************************************/

int broadcast_neighbors()
{
  int num_neighs = 0;
  neigh_t neigh;
  atom_t* atom = NULL;

  memset(&neigh, 0, sizeof(neigh));

  for (int i = 0; i < g_config.natoms; ++i) {
    atom = g_config.conf_atoms + i - g_mpi.firstatom;
    if (g_mpi.myid == 0)
      num_neighs = g_config.atoms[i].num_neigh;
    CHECK_RETURN(MPI_Bcast(&num_neighs, 1, MPI_INT, 0, MPI_COMM_WORLD));
    if (num_neighs > 0 && i >= g_mpi.firstatom &&
        i < (g_mpi.firstatom + g_mpi.myatoms)) {
      atom->neigh = (neigh_t*)Malloc(num_neighs * sizeof(neigh_t));
      for (int j = 0; j < num_neighs; ++j)
        memset(atom->neigh + j, 0, sizeof(neigh_t));
    }
    for (int j = 0; j < num_neighs; ++j) {
      if (g_mpi.myid == 0)
        neigh = g_config.atoms[i].neigh[j];
      CHECK_RETURN(MPI_Bcast(&neigh, 1, g_mpi.MPI_NEIGH, 0, MPI_COMM_WORLD));
      if (i >= g_mpi.firstatom && i < (g_mpi.firstatom + g_mpi.myatoms))
        atom->neigh[j] = neigh;
    }
  }

  return MPI_SUCCESS;
}

/***************************************************************************
    scatter dynamic angle table
**************************************************************************/

int broadcast_angles()
{
#if defined(THREEBODY)
  int num_angles = 0;
  angle_t angle;
  atom_t* atom = NULL;

  memset(&angle, 0, sizeof(angle));

  for (int i = 0; i < g_config.natoms; ++i) {
    atom = g_config.conf_atoms + i - g_mpi.firstatom;
    if (g_mpi.myid == 0)
      num_angles = g_config.atoms[i].num_angles;
    CHECK_RETURN(MPI_Bcast(&num_angles, 1, MPI_INT, 0, MPI_COMM_WORLD));
#if defined(ANG)
    if (num_angles == 0)
	    continue;
#endif
    if (i >= g_mpi.firstatom && i < (g_mpi.firstatom + g_mpi.myatoms)) {
      atom->angle_part = (angle_t*)Malloc(num_angles * sizeof(angle_t));
      for (int j = 0; j < num_angles; ++j)
        memset(atom->angle_part + j, 0, sizeof(angle_t));
    }
    for (int j = 0; j < num_angles; ++j) {
      if (g_mpi.myid == 0)
        angle = g_config.atoms[i].angle_part[j];
      CHECK_RETURN(MPI_Bcast(&angle, 1, g_mpi.MPI_ANGL, 0, MPI_COMM_WORLD));
      if (i >= g_mpi.firstatom && i < (g_mpi.firstatom + g_mpi.myatoms))
        atom->angle_part[j] = angle;
    }
  }
#endif  // THREEBODY

  return MPI_SUCCESS;
}

#endif  // MPI
