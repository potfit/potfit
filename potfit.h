/****************************************************************
 *
 * potfit.h: potfit header file
 *
 ****************************************************************
 *
 * Copyright 2002-2010 Peter Brommer, Franz G"ahler, Daniel Schopf
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

#define NRANSI

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifdef MPI
#include <mpi.h>
#define REAL MPI_DOUBLE
#endif /* MPI */

#include "random.h"

#ifdef APOT
#define APOT_STEPS 200		/* number of sampling points for analytic pot */
#define APOT_PUNISH 10e6	/* general value for apot punishments */
#endif

#if defined EAM || defined ADP
#define DUMMY_WEIGHT 100.
#endif /* EAM || ADP */

#define FORCE_EPS .1

#if defined PAIR
#define SLOTS 1
#elif defined EAM
#define SLOTS 2
#elif defined ADP
#define SLOTS 4
#endif /* PAIR */

/****************************************************************
 *
 *  type definitions
 *
 ****************************************************************/

typedef double real;

typedef enum ParamType {
  PARAM_STR, PARAM_STRPTR,
  PARAM_INT, PARAM_DOUBLE
} PARAMTYPE;

typedef struct {
  real  x;
  real  y;
  real  z;
} vector;

/* This is the order of vasp for stresses */
typedef struct {
  real  xx;
  real  yy;
  real  zz;
  real  xy;
  real  yz;
  real  zx;
} sym_tens;

typedef struct {
  int   typ;
  int   nr;
  real  r;
  vector dist;			/* distance divided by r */
  int   slot[SLOTS];
  real  shift[SLOTS];
  real  step[SLOTS];
  int   col[SLOTS];		/* coloumn of interaction for this neighbor */
#ifdef ADP
  vector rdist;			/* real distance */
  sym_tens sqrdist;		/* real squared distance */
  real  u_val, u_grad;		/* value and gradient of u(r) */
  real  w_val, w_grad;		/* value and gradient of w(r) */
#endif				/* ADP */
} neigh_t;

typedef struct {
  int   typ;
  int   n_neigh;
  vector pos;
  vector force;
  real  absforce;
  int   conf;			/* Which configuration... */
#if defined EAM || defined ADP
  real  rho;			/* embedding electron density */
  real  gradF;			/* gradient of embedding fn. */
#endif				/* EAM || ADP */
#ifdef ADP
  vector mu;
  sym_tens lambda;
  real  nu;
#endif				/* ADP */
  neigh_t *neigh;		/* dynamic array for neighbors */
} atom_t;

typedef struct {
  int   len;			/* total length of the table */
  int   idxlen;			/* number of changeable potential values */
  int   ncols;			/* number of columns */
  real *begin;			/* first value in the table */
  real *end;			/* last value in the table */
  real *step;			/* table increment */
  real *invstep;		/* inverse of increment */
  int  *first;			/* index of first entry */
  int  *last;			/* index of last entry */
  real *xcoord;			/* the x-coordinates of sampling points */
  real *table;			/* the actual data */
  real *d2tab;			/* second derivatives of table data for spline int */
  int  *idx;			/* indirect indexing */
} pot_table_t;

#ifdef APOT
/* function pointer for analytic potential evaluation */
typedef void (*fvalue_pointer) (real, real *, real *);

typedef struct {
  /* potentials */
  int   number;			/* number of analytic potentials */
  int   invar_pots;		/* number of invariant analytic potentials */
  int  *idxpot;			/* indirect index for potentials */
  char **names;			/* name of analytic potentials */
  real *begin;			/* starting position of potential */
  real *end;			/* end position of potential = cutoff radius */
  int  *n_par;			/* number of parameters for analytic potential */

  /* parameters */
  int   total_par;		/* total number of parameters for all potentials */
  int  *idxparam;		/* indirect index for potential parameters */
  int **invar_par;		/* array of invariant parameters */
  char ***param_name;		/* name of parameters */
  real **pmin;			/* minimum values for parameters */
  real **values;		/* parameter values for analytic potentials */
  real **pmax;			/* maximum values for parameters */

  /* global parameters */
  int   globals;		/* number of global parameters */
  int  *n_glob;			/* number of global parameter usage */
  int ***global_idx;		/* index of global parameters */

#ifdef PAIR
  real *chempot;		/* chemical potentials */
#endif				/* PAIR */

  fvalue_pointer *fvalue;	/* function pointers for analytic potentials */
} apot_table_t;
#endif /* APOT */

#define MAX(a,b)   ((a) > (b) ? (a) : (b))
#define MIN(a,b)   ((a) < (b) ? (a) : (b))
#define SPROD(a,b) (((a).x * (b).x) + ((a).y * (b).y) + ((a).z * (b).z))
#define SQR(a)     ((a)*(a))
#define SWAP(A,B,C) (C)=(A);(A)=(B);(B)=(C);

/****************************************************************
 *
 *  global variables
 *
 ****************************************************************/

/* MAIN is defined only once in the main module */
#ifdef MAIN
#define EXTERN			/* define Variables in main */
#define INIT(data) =data	/* initialize data only in main */
#else
#define EXTERN extern		/* declare them extern otherwise */
#define INIT(data)		/* skip initialization otherwise */
#endif /* MAIN */

/* system variables */
EXTERN int myid INIT(0);	/* Who am I? (0 if serial) */
EXTERN int num_cpus INIT(1);	/* How many cpus are there */
#ifdef MPI
EXTERN MPI_Datatype MPI_ATOM;
EXTERN MPI_Datatype MPI_NEIGH;
EXTERN MPI_Datatype MPI_TRANSMIT_NEIGHBOR;
EXTERN MPI_Datatype MPI_STENS;
EXTERN MPI_Datatype MPI_VEKTOR;
#endif /* MPI */

/* general settings (from parameter file) */
EXTERN char config[255] INIT("\0");	/* file with atom configuration */
EXTERN char distfile[255] INIT("\0");	/* file for distributions */
EXTERN char endpot[255] INIT("\0");	/* file for end potential */
EXTERN char flagfile[255] INIT("\0");	/* break if file exists */
EXTERN char imdpot[255] INIT("\0");	/* file for IMD potential */
EXTERN char maxchfile[255] INIT("\0");	/* file with maximal changes */
EXTERN char output_prefix[255] INIT("\0");	/* prefix for all output files */
EXTERN char plotfile[255] INIT("\0");	/* file for plotting */
EXTERN char plotpointfile[255] INIT("\0");	/* write points for plotting */
EXTERN char startpot[255] INIT("\0");	/* file with start potential */
EXTERN char tempfile[255] INIT("\0");	/* backup potential file */
EXTERN int imdpotsteps INIT(1000);	/* resolution of IMD potential */
EXTERN int ntypes INIT(1);	/* number of atom types */
EXTERN int opt INIT(0);		/* optimization flag */
EXTERN int seed INIT(4);	/* seed for RNG */
EXTERN int usemaxch INIT(0);	/* use maximal changes file */
EXTERN int write_output_files INIT(0);
EXTERN int write_pair INIT(0);
EXTERN int writeimd INIT(0);
#ifdef SIMANN
EXTERN real anneal_temp INIT(1.);
#endif
EXTERN real evo_threshold INIT(1.e-6);
EXTERN real eweight INIT(100.);
EXTERN real extend INIT(2.);	/* how far should one extend imd pot */
EXTERN real sweight INIT(10.);
#ifdef APOT
EXTERN int compnodes INIT(0);	/* how many additional composition nodes */
EXTERN int enable_cp INIT(0);	/* switch chemical potential on/off */
EXTERN real apot_punish_value INIT(0.);
EXTERN real plotmin INIT(0.);	/* minimum for plotfile */
#endif /* APOT */

/* configurations */
EXTERN atom_t *atoms INIT(NULL);	/* atoms array */
EXTERN atom_t *conf_atoms INIT(NULL);	/* Atoms in configuration */
EXTERN char **elements INIT(NULL);	/* element names from vasp2force */
EXTERN int **na_type INIT(NULL);	/* number of atoms per type */
EXTERN int *cnfstart INIT(NULL);	/* Nr. of first atom in config */
EXTERN int *conf_uf INIT(NULL);
EXTERN int *conf_us INIT(NULL);
EXTERN int *inconf INIT(NULL);	/* Nr. of atoms in each config */
EXTERN int *useforce INIT(NULL);	/* Should we use force/stress */
EXTERN int *usestress INIT(NULL);	/* Should we use force/stress */
EXTERN int have_elements INIT(0);	/* do we have the elements ? */
EXTERN int maxneigh INIT(0);	/* maximum number of neighbors */
EXTERN int natoms INIT(0);	/* number of atoms */
EXTERN int nconf INIT(0);	/* number of configurations */
EXTERN real *coheng INIT(NULL);	/* Cohesive energy for each config */
EXTERN real *conf_vol INIT(NULL);
EXTERN real *conf_weight INIT(NULL);	/* weight of configuration */
EXTERN real *force_0 INIT(NULL);	/* the forces we aim at */
EXTERN real *rcut INIT(NULL);
EXTERN real *rmin INIT(NULL);
EXTERN real *volumen INIT(NULL);	/* Volume of cell */
EXTERN real rcutmax INIT(0.);	/* maximum of all cutoff values */
EXTERN sym_tens *conf_stress INIT(NULL);
EXTERN sym_tens *stress INIT(NULL);	/* Stresses in each config */
EXTERN vector box_x, box_y, box_z;
EXTERN vector tbox_x, tbox_y, tbox_z;

/* potential variables */
EXTERN char interaction[10] INIT("\0");
EXTERN int *gradient INIT(NULL);	/* Gradient of potential fns.  */
EXTERN int *invar_pot INIT(NULL);
EXTERN int format INIT(-1);	/* format of potential table */
EXTERN int have_grad INIT(0);	/* Is gradient specified?  */
EXTERN int have_invar INIT(0);	/* Are invariant pots specified?  */
#ifdef APOT
EXTERN int *smooth_pot INIT(NULL);
EXTERN int cp_start INIT(0);	/* cp in opt_pot.table */
EXTERN int do_smooth INIT(0);	/* smooth cutoff option enabled? */
EXTERN int global_idx INIT(0);	/* index for global parameters in opt_pot table */
EXTERN int global_pot INIT(0);	/* number of "potential" for global parameters */
EXTERN int have_globals INIT(0);	/* do we have global parameters? */
EXTERN real *calc_list INIT(NULL);	/* list of current potential in the calc table */
EXTERN real *compnodelist INIT(NULL);	/* list of the composition nodes */
#endif /* APOT */

/* potential tables */
EXTERN pot_table_t opt_pot;	/* potential in the internal */
				/* representation used for  */
				/* minimisation */
EXTERN pot_table_t calc_pot;	/* the potential table used */
				/* for force calculations */
#ifdef APOT
EXTERN apot_table_t apot_table;	/* potential in analytic form */
#endif /* APOT */

/* optimization variables */
EXTERN int fcalls INIT(0);
EXTERN int mdim INIT(0);
EXTERN int ndim INIT(0);
EXTERN int ndimtot INIT(0);
EXTERN int paircol INIT(0);	/* How manc columns for pair potential */
EXTERN real d_eps INIT(0.);

/* general variables */
EXTERN int firstatom INIT(0);
EXTERN int firstconf INIT(0);
EXTERN int myatoms INIT(0.);
EXTERN int myconf INIT(0.);
EXTERN real *rms INIT(NULL);

/* pointers for force-vector */
EXTERN int energy_p INIT(0);	/* pointer to energies */
EXTERN int stress_p INIT(0);	/* pointer to stresses */
#if defined EAM || defined ADP
EXTERN int dummy_p INIT(0);	/* pointer to dummy constraints */
EXTERN int limit_p INIT(0);	/* pointer to limiting constraints */
#endif /* EAM || ADP */
#ifdef APOT
EXTERN int punish_par_p INIT(0);	/* pointer to parameter punishment contraints */
EXTERN int punish_pot_p INIT(0);	/* pointer to potential punishment constraints */
#endif /* APOT */

/* memory management */
EXTERN char **pointer_names INIT(NULL);
EXTERN int num_pointers INIT(0);
EXTERN void **all_pointers INIT(NULL);

/* variables needed for atom distribution with mpi */
#ifdef MPI
EXTERN int *atom_dist INIT(NULL);
EXTERN int *atom_len INIT(NULL);
EXTERN int *conf_dist INIT(NULL);
EXTERN int *conf_len INIT(NULL);
#endif /* MPI */

/* misc. stuff - has to belong somewhere */
EXTERN int *idx INIT(NULL);
EXTERN int init_done INIT(0);
EXTERN int plot INIT(0);	/* plot output flag */
EXTERN real *lambda INIT(NULL);	/* embedding energy slope... */
EXTERN real *maxchange INIT(NULL);	/* Maximal permissible change */
EXTERN dsfmt_t dsfmt;		/* random number generator */
EXTERN char *component[6];	/* componentes of vectors and tensors */

/****************************************************************
 *
 *  global function pointers
 *
 ****************************************************************/

EXTERN real (*calc_forces) (real *, real *, int);
EXTERN real (*splint) (pot_table_t *, real *, int, real);
EXTERN real (*splint_grad) (pot_table_t *, real *, int, real);
EXTERN real (*splint_comb) (pot_table_t *, real *, int, real, real *);
#ifdef APOT
EXTERN void (*write_pot_table) (apot_table_t *, char *);
#else
EXTERN void (*write_pot_table) (pot_table_t *, char *);
#endif /* APOT */
#ifdef PARABEL
EXTERN real (*parab) (pot_table_t *, real *, int, real);
EXTERN real (*parab_comb) (pot_table_t *, real *, int, real, real *);
EXTERN real (*parab_grad) (pot_table_t *, real *, int, real);
#endif /* PARABEL */

/****************************************************************
 *
 *  function prototypes
 *
 ****************************************************************/

/* general functions [potfit.c] */
void  error(char *);
void  warning(char *);

/* reading parameter file [param.c] */
int   getparam(char *, void *, PARAMTYPE, int, int);
void  read_parameters(int, char **);
void  read_paramfile(FILE *);

/* reading potential file [potential.c] */
void  read_pot_table(pot_table_t *, char *);
#ifdef APOT
void  read_pot_table0(pot_table_t *, apot_table_t *, char *, FILE *);
#endif /* APOT */
void  read_pot_table3(pot_table_t *, int, int, int *, char *, FILE *);
void  read_pot_table4(pot_table_t *, int, int, int *, char *, FILE *);

/* calculating potential tables [potential.c] */
void  init_calc_table(pot_table_t *, pot_table_t *);
void  update_calc_table(real *, real *, int);

/* parabolic interpolation [potential.c] */
#ifdef PARABEL
real  parab_comb_ed(pot_table_t *, real *, int, real, real *);
real  parab_grad_ed(pot_table_t *, real *, int, real);
real  parab_ed(pot_table_t *, real *, int, real);
real  parab_comb_ne(pot_table_t *, real *, int, real, real *);
real  parab_grad_ne(pot_table_t *, real *, int, real);
real  parab_ne(pot_table_t *, real *, int, real);
#endif /* PARABEL */

/* writing potentials to files [potential.c] */
#ifdef APOT
void  write_pot_table0(apot_table_t *, char *);
#endif /* APOT */
void  write_pot_table3(pot_table_t *, char *);
void  write_pot_table4(pot_table_t *, char *);
void  write_pot_table_imd(pot_table_t *, char *);
void  write_plotpot_pair(pot_table_t *, char *);
void  write_altplot_pair(pot_table_t *, char *);
#ifdef PDIST
void  write_pairdist(pot_table_t *, char *);
#endif /* PDIST */

/* read atomic configuration file [config.c] */
real  make_box(void);
void  read_config(char *);
#ifdef APOT
void  update_slots(void);	/* new slots */
#endif /* APOT */

/* force routines for different potential models [force_xxx.c] */
#ifdef PAIR
real  calc_forces_pair(real *, real *, int);
#elif defined EAM
real  calc_forces_eam(real *, real *, int);
#elif defined ADP
real  calc_forces_adp(real *, real *, int);
#endif /* PAIR */

/* simulated annealing [simann.c] */
#ifdef SIMANN
#ifdef APOT
void  randomize_parameter(int, real *, real *);
#else
void  makebump(real *, real, real, int);
#endif /* APOT */
void  anneal(real *);
#endif /* SIMANN */

/* powell least squares [powell_lsq.c] */
void  powell_lsq(real *);
int   gamma_init(real **, real **, real *, real *);
int   gamma_update(real **, real, real, real *, real *, real *, int, int, int,
  real);
void  lineqsys_init(real **, real **, real *, real *, int, int);
void  lineqsys_update(real **, real **, real *, real *, int, int, int);
void  copy_matrix(real **, real **, int, int);
void  copy_vector(real *, real *, int);
void  matdotvec(real **, real *, real *, int, int);
real  normalize_vector(real *, int);

/* differential evolution [diff_evo.c] */
void  init_population(real **, real *, real *);
#ifdef APOT
void  opposite_check(real **, real *, int);
#endif /* APOT */
void  diff_evo(real *);

/* spline interpolation [splines.c] */
void  spline_ed(real, real *, int, real, real, real *);
real  splint_ed(pot_table_t *, real *, int, real);
real  splint_grad_ed(pot_table_t *, real *, int, real);
real  splint_comb_ed(pot_table_t *, real *, int, real, real *);
real  splint_dir(pot_table_t *, real *, int, real, real);
real  splint_comb_dir(pot_table_t *, real *, int, real, real, real *);
real  splint_grad_dir(pot_table_t *, real *, int, real, real);
void  spline_ne(real *, real *, int, real, real, real *);
real  splint_ne(pot_table_t *, real *, int, real);
real  splint_comb_ne(pot_table_t *, real *, int, real, real *);
real  splint_grad_ne(pot_table_t *, real *, int, real);

/* rescaling functions for EAM [rescale.c] */
#ifdef EAM
real  rescale(pot_table_t *, real, int);
void  embed_shift(pot_table_t *);
#endif /* EAM */

/* MPI parallelization [mpi_utils.c] */
#ifdef MPI
void  init_mpi(int, char **);
void  shutdown_mpi(void);
void  broadcast_params(void);
void  debug_mpi(int);
void  broadcast_neighbors(void);
void  potsync(void);
#endif /* MPI */

#ifdef APOT
/* analytic functions [functions.c] */
int   apot_parameters(char *);
int   apot_assign_functions(apot_table_t *);
int   apot_check_params(real *);
real  apot_punish(real *, real *);
real  apot_grad(real, real *, void (*function) (real, real *, real *));
real  cutoff(real, real, real);
#ifdef DEBUG
void  debug_apot();
#endif /* DEBUG */

#ifdef PAIR
/* chemical potential [chempot.c] */
int   swap_chem_pot(int, int);
int   sort_chem_pot_2d(void);
real  chemical_potential_1d(int *, real *);
real  chemical_potential_2d(int *, real *);
real  chemical_potential_3d(int *, real *, int);
real  chemical_potential(int, int *, real *);
void  init_chemical_potential(int);
#endif /* PAIR */

/* actual functions for different potentials [functions.c] */

void  lj_value(real, real *, real *);
void  eopp_value(real, real *, real *);
void  morse_value(real, real *, real *);
void  ms_value(real, real *, real *);
void  softshell_value(real, real *, real *);
void  eopp_exp_value(real, real *, real *);
void  meopp_value(real, real *, real *);
void  power_decay_value(real, real *, real *);
void  exp_decay_value(real, real *, real *);
void  pohlong_value(real, real *, real *);
void  parabola_value(real, real *, real *);
void  csw_value(real, real *, real *);
void  universal_value(real, real *, real *);
void  const_value(real, real *, real *);
void  sqrt_value(real, real *, real *);
void  mexp_decay_value(real, real *, real *);
void  strmm_value(real, real *, real *);
void  double_morse_value(real, real *, real *);
void  double_exp_value(real, real *, real *);
void  poly_5_value(real, real *, real *);
void  cbb_value(real, real *, real *);
void  exp_plus_value(real, real *, real *);
void  mishin_value(real, real *, real *);
void  gen_lj_value(real, real *, real *);
void  gljm_value(real, real *, real *);

/* template for new potential function called newpot */

/* newpot potential */
void  newpot_value(real, real *, real *);
/* end of template */

#endif /* APOT */
