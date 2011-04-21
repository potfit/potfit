/****************************************************************
 *
 * potfit.h: potfit header file
 *
 ****************************************************************
 *
 * Copyright 2002-2011 Peter Brommer, Franz G"ahler, Daniel Schopf
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

#define POTFIT_H

#include <stdarg.h>
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
#define APOT_STEPS 300		/* number of sampling points for analytic pot */
#define APOT_PUNISH 10e6	/* general value for apot punishments */
#endif /* APOT */

#if defined EAM || defined ADP
#define DUMMY_WEIGHT 100.
#endif /* EAM || ADP */

#define FORCE_EPS .1

#if defined PAIR || defined COULOMB
#define SLOTS 1			/* pair potential = 0 */
#elif defined EAM		/* transfer function = 1 */
#define SLOTS 2			/* dipole term = 2 */
#elif defined ADP		/* quadrupole term = 3 */
#define SLOTS 4
#endif /* PAIR || COULOMB */

/****************************************************************
 *
 *  type definitions
 *
 ****************************************************************/

typedef double real;

typedef enum Param_T { PARAM_STR, PARAM_INT, PARAM_DOUBLE } param_t;

typedef enum Interaction_T { I_PAIR, I_EAM, I_ADP, I_ELSTAT } Interaction_T;

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
#ifdef COULOMB
  real  r2;			/* r^2 */
  real  fnval_el;		/* stores tail of electrostatic potential */
  real  grad_el;		/* stores tail of first derivative of electrostatic potential */
  real  ggrad_el;		/* stores tail of second derivative of electrostatic potential */
#endif
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
  vector E_stat;		/* static field-contribution */
  vector p_sr;			/* short-range dipole moment */
  vector E_ind;			/* induced field-contribution */
  vector p_ind;			/* induced dipole moment */
  vector E_old;			/* stored old induced field */
  vector E_tot;			/* temporary induced field */
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
#ifdef COULOMB
  int   total_ne_par;		/* total number of non-electrostatic parameters */
#endif
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
#endif

#ifdef COULOMB
  real *ratio;			/* stoichiometric ratio */
  real *charge;			/* charges */
  real  last_charge;		/* last charge determined on the basis of charge neutrality */
#endif
#ifdef DIPOLE
  real *dp_alpha;		/* polarisability */
  real *dp_b;			/* parameter for short-range-dipole-moment */
  real *dp_c;			/* parameter for short-range-dipole-moment */
  real *dp_kappa;        	/* parameter kappa */
  int   sum_t;			/* dipole-convergency output */
  int   sw_kappa;               /* switch for kappa-optimization */
#endif

  fvalue_pointer *fvalue;	/* function pointers for analytic potentials */
} apot_table_t;

typedef struct {
  char **name;			/* identifier of the potential */
  int  *n_par;			/* number of parameters */
  fvalue_pointer *fvalue;	/* function pointer */
} function_table_t;

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
#define INIT(data) = data	/* initialize data only in main */
#else
#define EXTERN extern		/* declare them extern otherwise */
#define INIT(data)		/* skip initialization otherwise */
#endif /* MAIN */

/* define interaction type */
#ifdef PAIR
EXTERN Interaction_T interaction INIT(I_PAIR);
#elif defined EAM
EXTERN Interaction_T interaction INIT(I_EAM);
#elif defined ADP
EXTERN Interaction_T interaction INIT(I_ADP);
#elif defined COULOMB
EXTERN Interaction_T interaction INIT(I_ELSTAT);
#endif /* interaction type */

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
EXTERN int ntypes INIT(-1);	/* number of atom types */
EXTERN int opt INIT(0);		/* optimization flag */
EXTERN int seed INIT(4);	/* seed for RNG */
EXTERN int usemaxch INIT(0);	/* use maximal changes file */
EXTERN int write_output_files INIT(0);
EXTERN int write_pair INIT(0);
EXTERN int writeimd INIT(0);
#ifdef EVO
EXTERN real evo_threshold INIT(1.e-6);
#else /* EVO */
EXTERN char anneal_temp[20] INIT("\0");
#endif /* EVO */
EXTERN real eweight INIT(-1.);
EXTERN real sweight INIT(-1.);
EXTERN real extend INIT(2.);	/* how far should one extend imd pot */
#ifdef APOT
EXTERN int compnodes INIT(0);	/* how many additional composition nodes */
EXTERN int enable_cp INIT(0);	/* switch chemical potential on/off */
EXTERN real apot_punish_value INIT(0.);
EXTERN real plotmin INIT(0.);	/* minimum for plotfile */
#endif /* APOT */

/* configurations */
EXTERN atom_t *atoms;		/* atoms array */
EXTERN atom_t *conf_atoms;	/* Atoms in configuration */
EXTERN char **elements;		/* element names from vasp2force */
EXTERN int **na_type;		/* number of atoms per type */
EXTERN int *cnfstart;		/* Nr. of first atom in config */
EXTERN int *conf_uf;
EXTERN int *conf_us;
EXTERN int *inconf;		/* Nr. of atoms in each config */
EXTERN int *useforce;		/* Should we use force/stress */
EXTERN int *usestress;		/* Should we use force/stress */
EXTERN int have_elements INIT(0);	/* do we have the elements ? */
EXTERN int maxneigh INIT(0);	/* maximum number of neighbors */
EXTERN int natoms INIT(0);	/* number of atoms */
EXTERN int nconf INIT(0);	/* number of configurations */
EXTERN real *coheng;		/* Cohesive energy for each config */
EXTERN real *conf_vol;
EXTERN real *conf_weight;	/* weight of configuration */
EXTERN real *force_0;		/* the forces we aim at */
EXTERN real *rcut;
EXTERN real *rmin;
EXTERN real *volumen;		/* Volume of cell */
EXTERN real rcutmax INIT(0.);	/* maximum of all cutoff values */
EXTERN sym_tens *conf_stress;
EXTERN sym_tens *stress;	/* Stresses in each config */
EXTERN vector box_x, box_y, box_z;
EXTERN vector tbox_x, tbox_y, tbox_z;

/* potential variables */
EXTERN char interaction_name[10] INIT("\0");
EXTERN int *gradient;		/* Gradient of potential fns.  */
EXTERN int *invar_pot;
EXTERN int format INIT(-1);	/* format of potential table */
EXTERN int have_grad INIT(0);	/* Is gradient specified?  */
EXTERN int have_invar INIT(0);	/* Are invariant pots specified?  */
#ifdef APOT
EXTERN int *smooth_pot;
EXTERN int cp_start INIT(0);	/* cp in opt_pot.table */
EXTERN int global_idx INIT(0);	/* index for global parameters in opt_pot table */
EXTERN int global_pot INIT(0);	/* number of "potential" for global parameters */
EXTERN int have_globals INIT(0);	/* do we have global parameters? */
EXTERN real *calc_list;		/* list of current potential in the calc table */
EXTERN real *compnodelist;	/* list of the composition nodes */
#endif /* APOT */

/* potential tables */
EXTERN pot_table_t opt_pot;	/* potential in the internal */
				/* representation used for  */
				/* minimisation */
EXTERN pot_table_t calc_pot;	/* the potential table used */
				/* for force calculations */
#ifdef APOT
EXTERN apot_table_t apot_table;	/* potential in analytic form */
EXTERN int n_functions INIT(0);	/* number of analytic function prototypes */
EXTERN function_table_t function_table;	/* table with all functions */
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
EXTERN int myatoms INIT(0);
EXTERN int myconf INIT(0);

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
EXTERN char **pointer_names;
EXTERN int num_pointers INIT(0);
EXTERN void **all_pointers;

/* variables needed for atom distribution with mpi */
#ifdef MPI
EXTERN int *atom_dist;
EXTERN int *atom_len;
EXTERN int *conf_dist;
EXTERN int *conf_len;
#endif /* MPI */

/* misc. stuff - has to belong somewhere */
EXTERN int *idx;
EXTERN int init_done INIT(0);
EXTERN int plot INIT(0);	/* plot output flag */
EXTERN real *lambda;		/* embedding energy slope... */
EXTERN real *maxchange;		/* Maximal permissible change */
EXTERN dsfmt_t dsfmt;		/* random number generator */
EXTERN char *component[6];	/* componentes of vectors and tensors */

/* variables needed for electrostatic options */
#ifdef COULOMB
EXTERN real dp_eps INIT(14.40);	/* this is e^2/(4*pi*epsilon_0) in eV A */
EXTERN real dp_cut INIT(10);	/* cutoff-radius for long-range interactions */
#endif /* COULOMB */
#ifdef DIPOLE
EXTERN real dp_tol INIT(1.e-7);	/* dipole iteration precision */
EXTERN real dp_mix INIT(0.2);	/* mixing parameter (other than that one in IMD) */
#endif /* DIPOLE */

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
void  error(char *, ...);
void  warning(char *, ...);

/* reading parameter file [param.c] */
int   getparam(char *, void *, param_t, int, int);
void  check_parameters_complete(char *);
void  read_parameters(int, char **);
void  read_paramfile(FILE *);

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
#elif defined COULOMB
real  calc_forces_elstat(real *, real *, int);
#endif /* interaction type */

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
void  broadcast_neighbors(void);
void  potsync(void);
#endif /* MPI */
