/****************************************************************
 *
 * potfit.h: potfit header file
 *
 ****************************************************************
 *
 * Copyright 2002-2012
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.itap.physik.uni-stuttgart.de/
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

#ifdef MPI
#include <mpi.h>
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

/****************************************************************
 *
 *  SLOTS: number of different distance tables used in the force calculations
 *
 *  pair potential 	= 	0 		PAIR pair distance
 *  transfer function 	= 	1  		EAM transfer function
 *  dipole term 	= 	2 		ADP dipole term
 *  quadrupole term 	= 	3 		ADP quadrupole term
 *
 ****************************************************************/

#define SLOTS 1

#if defined EAM
#undef SLOTS
#define SLOTS 2
#elif defined ADP
#undef SLOTS
#define SLOTS 4
#endif

/****************************************************************
 *
 *  type definitions
 *
 ****************************************************************/

typedef enum Param_T { PARAM_STR, PARAM_INT, PARAM_DOUBLE } param_t;

typedef enum Interaction_T { I_PAIR, I_EAM, I_ADP, I_ELSTAT } Interaction_T;

typedef struct {
  double x;
  double y;
  double z;
} vector;

/* This is the order of VASP for stresses */
typedef struct {
  double xx;
  double yy;
  double zz;
  double xy;
  double yz;
  double zx;
} sym_tens;

typedef struct {
  int   typ;
  int   nr;
  double r;
  vector dist;			/* distance divided by r */
  int   slot[SLOTS];
  double shift[SLOTS];
  double step[SLOTS];
  int   col[SLOTS];		/* coloumn of interaction for this neighbor */
#ifdef ADP
  vector rdist;			/* real distance */
  sym_tens sqrdist;		/* real squared distance */
  double u_val, u_grad;		/* value and gradient of u(r) */
  double w_val, w_grad;		/* value and gradient of w(r) */
#endif
#ifdef COULOMB
  double r2;			/* r^2 */
  double fnval_el;		/* stores tail of electrostatic potential */
  double grad_el;		/* stores tail of first derivative of electrostatic potential */
  double ggrad_el;		/* stores tail of second derivative of electrostatic potential */
#endif
} neigh_t;

typedef struct {
  int   typ;
  int   n_neigh;
  vector pos;
  vector force;
  double absforce;
  int   conf;			/* Which configuration... */
#if defined EAM || defined ADP
  double rho;			/* embedding electron density */
  double gradF;			/* gradient of embedding fn. */
#endif
#ifdef ADP
  vector mu;
  sym_tens lambda;
  double nu;
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
  double *begin;		/* first value in the table */
  double *end;			/* last value in the table */
  double *step;			/* table increment */
  double *invstep;		/* inverse of increment */
  int  *first;			/* index of first entry */
  int  *last;			/* index of last entry */
  double *xcoord;		/* the x-coordinates of sampling points */
  double *table;		/* the actual data */
  double *d2tab;		/* second derivatives of table data for spline int */
  int  *idx;			/* indirect indexing */
} pot_table_t;

#ifdef APOT
/* function pointer for analytic potential evaluation */
typedef void (*fvalue_pointer) (double, double *, double *);

typedef struct {
  /* potentials */
  int   number;			/* number of analytic potentials */
  int   invar_pots;		/* number of invariant analytic potentials */
  int  *idxpot;			/* indirect index for potentials */
  char **names;			/* name of analytic potentials */
  double *begin;		/* starting position of potential */
  double *end;			/* end position of potential = cutoff radius */
  int  *n_par;			/* number of parameters for analytic potential */

  /* parameters */
  int   total_par;		/* total number of parameters for all potentials */
#ifdef COULOMB
  int   total_ne_par;		/* total number of non-electrostatic parameters */
#endif
  int  *idxparam;		/* indirect index for potential parameters */
  int **invar_par;		/* array of invariant parameters */
  char ***param_name;		/* name of parameters */
  double **pmin;		/* minimum values for parameters */
  double **values;		/* parameter values for analytic potentials */
  double **pmax;		/* maximum values for parameters */

  /* global parameters */
  int   globals;		/* number of global parameters */
  int  *n_glob;			/* number of global parameter usage */
  int ***global_idx;		/* index of global parameters */

#ifdef PAIR
  double *chempot;		/* chemical potentials */
#endif

#ifdef COULOMB
  double *ratio;		/* stoichiometric ratio */
  double *charge;		/* charges */
  double last_charge;		/* last charge determined on the basis of charge neutrality */
  double *dp_kappa;		/* parameter kappa */
  int   sw_kappa;		/* switch for kappa-optimization */
#endif
#ifdef DIPOLE
  double *dp_alpha;		/* polarisability */
  double *dp_b;			/* parameter for short-range-dipole-moment */
  double *dp_c;			/* parameter for short-range-dipole-moment */
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
EXTERN int write_lammps INIT(0);	/* write output also in LAMMPS format */
#ifdef EVO
EXTERN double evo_threshold INIT(1.e-6);
#else /* EVO */
EXTERN char anneal_temp[20] INIT("\0");
#endif /* EVO */
EXTERN double eweight INIT(-1.);
EXTERN double sweight INIT(-1.);
EXTERN double extend INIT(2.);	/* how far should one extend imd pot */
#ifdef APOT
EXTERN int compnodes INIT(0);	/* how many additional composition nodes */
EXTERN int enable_cp INIT(0);	/* switch chemical potential on/off */
EXTERN double apot_punish_value INIT(0.);
EXTERN double plotmin INIT(0.);	/* minimum for plotfile */
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
EXTERN double *coheng;		/* Cohesive energy for each config */
EXTERN double *conf_vol;
EXTERN double *conf_weight;	/* weight of configuration */
EXTERN double *force_0;		/* the forces we aim at */
EXTERN double *rcut;
EXTERN double *rmin;
EXTERN double *volumen;		/* Volume of cell */
EXTERN double rcutmin INIT(999.);	/* minimum of all cutoff values */
EXTERN double rcutmax INIT(0.);	/* maximum of all cutoff values */
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
EXTERN double *calc_list;	/* list of current potential in the calc table */
EXTERN double *compnodelist;	/* list of the composition nodes */
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
EXTERN double d_eps INIT(0.);

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
EXTERN double *lambda;		/* embedding energy slope... */
EXTERN double *maxchange;	/* Maximal permissible change */
EXTERN dsfmt_t dsfmt;		/* random number generator */
EXTERN char *component[6];	/* componentes of vectors and tensors */

/* variables needed for electrostatic options */
#ifdef COULOMB
EXTERN double dp_eps INIT(14.40);	/* this is e^2/(4*pi*epsilon_0) in eV A */
EXTERN double dp_cut INIT(10);	/* cutoff-radius for long-range interactions */
#endif /* COULOMB */
#ifdef DIPOLE
EXTERN double dp_tol INIT(1.e-7);	/* dipole iteration precision */
EXTERN double dp_mix INIT(0.2);	/* mixing parameter (other than that one in IMD) */
#endif /* DIPOLE */

/****************************************************************
 *
 *  global function pointers
 *
 ****************************************************************/

EXTERN double (*calc_forces) (double *, double *, int);
EXTERN double (*splint) (pot_table_t *, double *, int, double);
EXTERN double (*splint_grad) (pot_table_t *, double *, int, double);
EXTERN double (*splint_comb) (pot_table_t *, double *, int, double, double *);
#ifdef APOT
EXTERN void (*write_pot_table) (apot_table_t *, char *);
#else
EXTERN void (*write_pot_table) (pot_table_t *, char *);
#endif /* APOT */
#ifdef PARABEL
EXTERN double (*parab) (pot_table_t *, double *, int, double);
EXTERN double (*parab_comb) (pot_table_t *, double *, int, double, double *);
EXTERN double (*parab_grad) (pot_table_t *, double *, int, double);
#endif /* PARABEL */

/****************************************************************
 *
 *  function prototypes
 *
 ****************************************************************/

/* general functions [potfit.c] */
void  error(int, char *, ...);
void  warning(int, char *, ...);

/* reading parameter file [param.c] */
int   getparam(char *, void *, param_t, int, int);
void  check_parameters_complete(char *);
void  read_parameters(int, char **);
void  read_paramfile(FILE *);

/* read atomic configuration file [config.c] */
double make_box(void);
void  read_config(char *);
#ifdef APOT
void  update_slots(void);	/* new slots */
#endif /* APOT */

/* force routines for different potential models [force_xxx.c] */
#ifdef PAIR
double calc_forces_pair(double *, double *, int);
#elif defined EAM
double calc_forces_eam(double *, double *, int);
#elif defined ADP
double calc_forces_adp(double *, double *, int);
#elif defined COULOMB
double calc_forces_elstat(double *, double *, int);
#endif /* interaction type */

/* rescaling functions for EAM [rescale.c] */
#ifdef EAM
double rescale(pot_table_t *, double, int);
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
