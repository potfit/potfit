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

#ifdef OMP
#include <omp.h>
#endif

#ifdef MPI
#include <mpi.h>
#define REAL MPI_DOUBLE
#endif

#include "random.h"

#if defined EAM || defined ADP
#define DUMMY_WEIGHT 100.
#endif

#define FORCE_EPS .1

#if defined PAIR || defined DIPOLE
#define SLOTS 1
#elif defined EAM
#define SLOTS 2
#elif defined ADP
#define SLOTS 4
#endif

/******************************************************************************
*
*  type definitions
*
******************************************************************************/

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

// This is the order of vasp for stresses
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
#endif
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
#endif
#ifdef ADP
  vector mu;
  sym_tens lambda;
  real  nu;
#endif
#ifdef DIPOLE
  vector E_stat;                /* static field-contribution */
  vector p_stat;                /* short-range dipole moment */
#endif
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
#ifdef DIPOLE
  real *table_dipole;           /* static data for tail of coulomb potential */
#endif
  real *d2tab;			/* second derivatives of table data for spline int */
  int  *idx;			/* indirect indexing */
} pot_table_t;

#ifdef APOT
/* function pointer for analytic potential evaluation */
typedef void (*fvalue_pointer) (real, real *, real *);

typedef struct {
//   potentials
  int   number;			/* number of analytic potentials */
  int  *idxpot;			/* indirect index for potentials */
  char **names;			/* name of analytic potentials */
  real *begin;			/* starting position of potential */
  real *end;			/* end position of potential = cutoff radius */
  int  *n_par;			/* number of parameters for analytic potential */

//   parameters
  int   total_par;		/* total number of parameters for all potentials */
  int  *idxparam;		/* indirect index for potential parameters */
  int **invar_par;		/* array of invariant parameters */
  char ***param_name;		/* name of parameters */
  real **pmin;			/* minimum values for parameters */
  real **values;		/* parameter values for analytic potentials */
  real **pmax;			/* maximum values for parameters */

//   global parameters
  int   globals;		/* number of global parameters */
  int  *n_glob;			/* number of global parameter usage */
  int ***global_idx;		/* index of global parameters */

#ifdef PAIR
  real *chempot;		/* chemical potentials */
#endif

#ifdef DIPOLE
  real *charge;			/* charges */
  real *dp_alpha;		/* polarisability */
  real *dp_b;			/* parameter for short-range-dipole-moment */
  real *dp_c;			/* parameter for short-range-dipole-moment */
#endif

  fvalue_pointer *fvalue;	/* function pointers for analytic potentials */
} apot_table_t;
#endif

#define MAX(a,b)   ((a) > (b) ? (a) : (b))
#define MIN(a,b)   ((a) < (b) ? (a) : (b))
#define SQR(a)     ((a)*(a))
#define SPROD(a,b) (((a).x * (b).x) + ((a).y * (b).y) + ((a).z * (b).z))

/******************************************************************************
*
*  global variables
*
******************************************************************************/

/* MAIN is defined only once in the main module */
#ifdef MAIN
#define EXTERN			/* define Variables in main */
#define INIT(data) =data	/* initialize data only in main */
#else
#define EXTERN extern		/* declare them extern otherwise */
#define INIT(data)		/* skip initialization otherwise */
#endif

// system variables
EXTERN int myid INIT(0);	/* Who am I? (0 if serial) */
EXTERN int num_cpus INIT(1);	/* How many cpus are there */
#ifdef MPI
EXTERN MPI_Datatype MPI_ATOM;
EXTERN MPI_Datatype MPI_NEIGH;
EXTERN MPI_Datatype MPI_TRANSMIT_NEIGHBOR;
EXTERN MPI_Datatype MPI_STENS;
EXTERN MPI_Datatype MPI_VEKTOR;
#endif

// general settings (from parameter file)
EXTERN char config[255];	/* file with atom configuration */
EXTERN char distfile[255];	/* file for distributions */
EXTERN char endpot[255];	/* file for end potential */
EXTERN char flagfile[255] INIT("\0");	/* break if file exists */
EXTERN char imdpot[255];	/* file for IMD potential */
EXTERN char maxchfile[255];	/* file with maximal changes */
EXTERN char output_prefix[255] INIT("\0");	/* prefix for all output files */
EXTERN char plotfile[255];	/* file for plotting */
EXTERN char plotpointfile[255] INIT("\0");	/* write points for plotting */
EXTERN char startpot[255];	/* file with start potential */
EXTERN char tempfile[255] INIT("\0");	/* backup potential file */
EXTERN int imdpotsteps;		/* resolution of IMD potential */
EXTERN int ntypes INIT(1);	/* number of atom types */
EXTERN int opt INIT(0);		/* optimization flag */
EXTERN int seed INIT(123456);	/* seed for RNG */
EXTERN int usemaxch INIT(0);	/* use maximal changes file */
EXTERN int write_output_files INIT(0);
EXTERN int write_pair INIT(0);
EXTERN int writeimd INIT(0);
EXTERN real anneal_temp INIT(1.);
EXTERN real eweight INIT(100.);
EXTERN real extend INIT(2.);	/* how far should one extend imd pot */
EXTERN real sweight INIT(10.);
#ifdef APOT
EXTERN int compnodes INIT(0);	/* how many additional composition nodes */
EXTERN int enable_cp INIT(0);	/* switch chemical potential on/off */
EXTERN real plotmin INIT(0.);	/* minimum for plotfile */
#endif
#ifdef EVO
EXTERN real evo_width INIT(1.);
#endif

// configurations
EXTERN atom_t *atoms INIT(NULL);	/* atoms array */
EXTERN atom_t *conf_atoms INIT(NULL);	/* Atoms in configuration */
EXTERN char **elements INIT(NULL);	/* element names from vasp2force */
EXTERN int **na_typ INIT(NULL);	/* number of atoms per type */
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

// potential variables
EXTERN char interaction[10] INIT("\0");
EXTERN int *gradient INIT(NULL);	/* Gradient of potential fns.  */
EXTERN int *invar_pot INIT(NULL);
EXTERN int format;		/* format of potential table */
EXTERN int have_grad INIT(0);	/* Is gradient specified?  */
EXTERN int have_invar INIT(0);	/* Are invariant pots specified?  */
#ifdef APOT
EXTERN int *pot_index INIT(NULL);	/* index to access i*n+j from i*(i+1)/2 */
EXTERN int *smooth_pot INIT(NULL);
EXTERN int cp_start INIT(0);	/* cp in opt_pot.table */
EXTERN int do_smooth INIT(0);	/* smooth cutoff option enabled? */
EXTERN int global_idx INIT(0);	/* index for global parameters in opt_pot table */
EXTERN int global_pot INIT(0);	/* number of "potential" for global parameters */
EXTERN int have_globals INIT(0);	/* do we have global parameters? */
EXTERN real *calc_list INIT(NULL);	/* list of current potential in the calc table */
EXTERN real *compnodelist INIT(NULL);	/* list of the composition nodes */
#endif

// potential tables
EXTERN pot_table_t opt_pot;	/* potential in the internal */
					 /* representation used for  */
					 /* minimisation */
EXTERN pot_table_t calc_pot;	/* the potential table used */
					 /* for force calculations */
#ifdef APOT
EXTERN apot_table_t apot_table;	/* potential in analytic form */
#endif

// optimization variables
EXTERN int fcalls INIT(0);
EXTERN int mdim INIT(0);
EXTERN int ndim INIT(0);
EXTERN int ndimtot INIT(0);
EXTERN int paircol INIT(0);	/* How manc columns for pair pot. */
EXTERN real d_eps INIT(0);

// general variables
EXTERN int firstatom INIT(0);
EXTERN int firstconf INIT(0);
EXTERN int myatoms INIT(0.);
EXTERN int myconf INIT(0.);
EXTERN real *rms INIT(NULL);

// pointers for force-vector
EXTERN int energy_p INIT(0);	/* pointer to energies */
EXTERN int stress_p INIT(0);	/* pointer to stresses */
#if defined EAM || defined ADP
EXTERN int dummy_p INIT(0);	/* pointer to dummy constraints */
EXTERN int limit_p INIT(0);	/* pointer to limiting constraints */
#endif
#ifdef APOT
EXTERN int punish_par_p INIT(0);	/* pointer to parameter punishment contraints */
EXTERN int punish_pot_p INIT(0);	/* pointer to potential punishment constraints */
#endif

// memory management
EXTERN char **pointer_names INIT(NULL);
EXTERN int num_pointers INIT(0);
EXTERN void **all_pointers INIT(NULL);

// variables needed for atom distribution with mpi
#ifdef MPI
EXTERN int *atom_dist;
EXTERN int *atom_len;
EXTERN int *conf_dist;
EXTERN int *conf_len;
#endif

// misc. stuff - has to belong somewhere
EXTERN int *idx INIT(NULL);
EXTERN int init_done INIT(0);
EXTERN int plot INIT(0);	/* plot output flag */
EXTERN real *lambda INIT(NULL);	/* embedding energy slope... */
EXTERN real *maxchange INIT(NULL);	/* Maximal permissible change */
EXTERN dsfmt_t dsfmt;		/* random number generator */
#ifdef STRESS
EXTERN char *stress_comp[6];
#endif

// variables needed for option dipole
#ifdef DIPOLE
EXTERN real dp_kappa INIT(0.10);	/* parameter kappa */
EXTERN real dp_tol INIT(1.e-9); /* dipole iteration precision */
EXTERN real dp_cut INIT(0);      /* cutoff-radius for long-range interactions */
EXTERN real dp_eps INIT(14.40);	/* this is e^2/(4*pi*epsilon_0) in eV A */
#endif

/******************************************************************************
*
*  global function pointers
*
******************************************************************************/

EXTERN real (*calc_forces) (real *, real *, int);
EXTERN real (*splint) (pot_table_t *, real *, int, real);
EXTERN real (*splint_grad) (pot_table_t *, real *, int, real);
EXTERN real (*splint_comb) (pot_table_t *, real *, int, real, real *);
#ifdef APOT
EXTERN void (*write_pot_table) (apot_table_t *, char *);
#else
EXTERN void (*write_pot_table) (pot_table_t *, char *);
#endif
#ifdef PARABEL
EXTERN real (*parab) (pot_table_t *, real *, int, real);
EXTERN real (*parab_comb) (pot_table_t *, real *, int, real, real *);
EXTERN real (*parab_grad) (pot_table_t *, real *, int, real);
#endif /* PARABEL */

/******************************************************************************
*
*  function prototypes
*
******************************************************************************/

void  error(char *);
void  warning(char *);
int   getparam(char *, void *, PARAMTYPE, int, int);
void  read_parameters(int, char **);
void  read_paramfile(FILE *);
#ifdef APOT
void  read_apot_table(pot_table_t *pt, apot_table_t *apt, char *filename,
		      FILE *infile);
void  write_apot_table(apot_table_t *, char *);
#endif
void  read_pot_table(pot_table_t *, char *);
void  read_pot_table3(pot_table_t *pt, int size, int ncols, int *nvals,
		      char *filename, FILE *infile);
void  read_pot_table4(pot_table_t *pt, int size, int ncols, int *nvals,
		      char *filename, FILE *infile);
void  init_calc_table(pot_table_t *optt, pot_table_t *calct);
void  update_calc_table(real *xi_opt, real *xi_calc, int);
void  write_pot_table3(pot_table_t *, char *);
void  write_pot_table4(pot_table_t *, char *);
void  write_pot_table_imd(pot_table_t *, char *);
void  write_plotpot_pair(pot_table_t *, char *);
void  write_altplot_pair(pot_table_t *, char *);
real  normdist(void);
real  grad2(pot_table_t *, real *, int, real);
real  grad3(pot_table_t *, real *, int, real);
real  pot2(pot_table_t *, int, real);
real  pot3(pot_table_t *, int, real);
vector vec_prod(vector, vector);
real  make_box(void);
void  read_config(char *);
void  read_config2(char *);
#ifdef PAIR
real  calc_forces_pair(real *, real *, int);
#elif defined EAM
real  calc_forces_eam(real *, real *, int);
#elif defined ADP
real  calc_forces_adp(real *, real *, int);
#elif defined DIPOLE
real  calc_forces_dipole(real *, real *, int);
#endif
#ifdef APOT
void  randomize_parameter(int, real *, real *);
#else
void  makebump(real *, real, real, int);
#endif
void  anneal(real *xi);
void  powell_lsq(real *xi);
#if defined EVO
real *calc_vect(real *x);
void  init_population(real **pop, real *xi, int size, real scale);
void  diff_evo(real *xi);
#endif
void  spline_ed(real xstep, real y[], int n, real yp1, real ypn, real y2[]);
real  splint_ed(pot_table_t *pt, real *xi, int col, real r);
real  splint_grad_ed(pot_table_t *pt, real *xi, int col, real r);
real  splint_comb_ed(pot_table_t *pt, real *xi, int col, real r, real *grad);
real  splint_dir(pot_table_t *pt, real *xi, int k, real b, real step);
real  splint_comb_dir(pot_table_t *pt, real *xi, int k, real b, real step,
		      real *grad);
real  splint_grad_dir(pot_table_t *pt, real *xi, int k, real b, real step);
void  spline_ne(real x[], real y[], int n, real yp1, real ypn, real y2[]);
real  splint_ne(pot_table_t *pt, real *xi, int col, real r);
real  splint_comb_ne(pot_table_t *pt, real *xi, int col, real r, real *grad);
real  splint_grad_ne(pot_table_t *pt, real *xi, int col, real r);

#ifdef PARABEL
real  parab_comb_ed(pot_table_t *pt, real *xi, int col, real r, real *grad);
real  parab_grad_ed(pot_table_t *pt, real *xi, int col, real r);
real  parab_ed(pot_table_t *pt, real *xi, int col, real r);
real  parab_comb_ne(pot_table_t *pt, real *xi, int col, real r, real *grad);
real  parab_grad_ne(pot_table_t *pt, real *xi, int col, real r);
real  parab_ne(pot_table_t *pt, real *xi, int col, real r);
#endif
#ifdef EAM
real  rescale(pot_table_t *pt, real upper, int flag);
void  embed_shift(pot_table_t *pt);
#endif
#ifdef MPI
void  init_mpi(int argc, char **argv);
void  shutdown_mpi(void);
void  broadcast_params(void);
void  debug_mpi(int i);
void  broadcast_neighbors();
void  potsync();
#endif
#ifdef PDIST
void  write_pairdist(pot_table_t *pt, char *filename);
#endif

/******************************************************************************
*
*  additional stuff for analytic potentials
*
******************************************************************************/

#ifdef APOT

#define APOT_STEPS 1000		/* number of sampling points for analytic pot */
#define APOT_PUNISH 10e6	/* general value for apot punishments */

real apot_punish_value INIT(0.);

/* functions.c */
int   apot_parameters(char *);
int   apot_assign_functions(apot_table_t *);
int   apot_check_params(real *);
real  apot_punish(real *, real *);
real  apot_grad(real, real *, void (*function) (real, real *, real *));

/* potential.c */
void  update_slots();		/* new slots for smooth cutoff */

#ifdef PAIR
/* chempot.c */
int   swap_chem_pot(int, int);
int   sort_chem_pot_2d(void);
real  chemical_potential_1d(int *, real *);
real  chemical_potential_2d(int *, real *);
real  chemical_potential_3d(int *, real *, int);
real  chemical_potential(int, int *, real *);
void  init_chemical_potential(int);
#endif

/* smooth.c */
real  cutoff(real, real, real);

#ifdef DEBUG
void  debug_apot();
#endif /* DEBUG */

/* actual functions for different potentials */

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

#ifdef DIPOLE
/* functions for dipole calculations  */
void ms_shift(real, real*, real*);
real shortrange_value(real, real *, real *, real *);
void coulomb_value(real, real *, real *);
void coulomb_shift(real, real*);
#endif

/* template for new potential function called newpot */

/* newpot potential */
void  newpot_value(real, real *, real *);
/* end of template */

#endif /* APOT */
