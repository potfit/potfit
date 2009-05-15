/****************************************************************
* 
*  potfit.h: potfit header file
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
* $Revision: 1.66 $
* $Date: 2009/05/15 16:28:50 $
*****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
//#include "utils.h"
#ifdef OMP
#include <omp.h>
#endif
#ifdef MPI
#include <mpi.h>
#endif
#ifdef MPI
#define REAL MPI_DOUBLE
#endif /* MPI */
#define NRANSI
#define MAXNEIGH 400
#ifdef EAM

#define DUMMY_WEIGHT 100.
#endif
#define ENG_WEIGHT 10.
#define STRESS_WEIGHT 10.
#define FORCE_EPS .1

/******************************************************************************
*
*  type definitions
*
******************************************************************************/

typedef double real;
typedef struct {
  real  x;
  real  y;
  real  z;
} vektor;
typedef struct {
  real  xx;
  real  yy;
  real  zz;
  real  xy;
  real  yz;
  real  zx;
} stens;

typedef struct {
  int   typ;
  int   nr;
  real  r;
  vektor dist;
  int   slot[2];
  real  shift[2];
  real  step[2];
} neigh_t;

typedef struct {
  int   typ;
  int   n_neigh;
  vektor pos;
  vektor force;
  real  absforce;
  neigh_t neigh[MAXNEIGH];
  int   conf;			/* Which configuration... */
#ifdef EAM
  real  rho;			/* embedding electron density */
  real  gradF;			/* gradient of embedding fn. */
#endif
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
  int   number;			/* number of analytic potentials */
  int   total_par;		/* total number of parameters for all potentials */
  int  *n_par;			/* number of parameters for analytic potential */
  char **names;			/* name of analytic potentials */
  char ***param_name;		/* name of parameters */
  real **values;		/* parameter values for analytic potentials */
  int **invar_par;		/* parameter values for analytic potentials */
  real *chempot;		/* chemical potentials */
  real *begin;			/* starting position of potential */
  real *end;			/* end position of potential = cutoff radius */
  real **pmin;			/* minimum values for parameters */
  real **pmax;			/* maximum values for parameters */
  int  *idxpot;			/* indirect index for potentials */
  int  *idxparam;		/* indirect index for potential parameters */
  fvalue_pointer *fvalue;	/* function pointers for analytic potentials */
} apot_table_t;
#endif

#define MAX(a,b)   ((a) > (b) ? (a) : (b))
#define MIN(a,b)   ((a) < (b) ? (a) : (b))
#define SQR(a)     ((a)*(a))
#define SPROD(a,b) (((a).x * (b).x) + ((a).y * (b).y) + ((a).z * (b).z))
static real sqrreal;
#define SQRREAL(x) ((sqrreal=(x)) == 0.0 ? 0.0 : sqrreal*sqrreal)

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
EXTERN int myid INIT(0);	/* Who am I? (0 if serial) */
EXTERN int num_cpus INIT(1);	/* How many cpus are there */
#ifdef MPI
EXTERN MPI_Datatype MPI_VEKTOR;
EXTERN MPI_Datatype MPI_STENS;
//EXTERN MPI_Datatype MPI_POTTABLE;
EXTERN MPI_Datatype MPI_ATOM;
EXTERN MPI_Datatype MPI_NEIGH;
#endif
EXTERN int *gradient INIT(NULL);	/* Gradient of potential fns.  */
EXTERN int have_grad INIT(0);	/* Is gradient specified?  */
EXTERN int *invar_pot INIT(NULL);
EXTERN int have_invar INIT(0);	/* Are invariant pots specified?  */
EXTERN int do_smooth INIT(0);	/* smooth cutoff option enabled? */
EXTERN int *smooth_pot INIT(NULL);
EXTERN int write_pair INIT(0);
EXTERN atom_t *conf_atoms INIT(NULL);	/* Atoms in configuration */
EXTERN real *conf_eng INIT(NULL);
EXTERN real *conf_vol INIT(NULL);
EXTERN stens *conf_stress INIT(NULL);
EXTERN int *atom_len;
EXTERN int *atom_dist;
EXTERN int *conf_len;
EXTERN int *conf_dist;
EXTERN int *conf_uf INIT(NULL);
EXTERN int *conf_us INIT(NULL);
EXTERN real sweight INIT(STRESS_WEIGHT);
EXTERN real eweight INIT(ENG_WEIGHT);
EXTERN int myconf INIT(0.);
EXTERN int myatoms INIT(0.);
EXTERN int firstconf INIT(0);
EXTERN int firstatom INIT(0);
EXTERN real *rms INIT(NULL);
EXTERN real pi INIT(0.);
EXTERN real extend INIT(2.);	/* how far should one extend imd pot */
EXTERN int writeimd INIT(0);
EXTERN real anneal_temp INIT(1.);
EXTERN int seed INIT(123456);	/* seed for RNG */
EXTERN int fcalls INIT(0);
EXTERN int ndim INIT(0);
EXTERN int ndimtot INIT(0);
EXTERN int mdim INIT(0);
EXTERN int usemaxch INIT(0);	/* use maximal changes file */
EXTERN int ntypes INIT(1);	/* number of atom types */
EXTERN int natoms INIT(0);	/* number of atoms */
EXTERN char **elements INIT(NULL);	/* element names from vasp2force */
EXTERN int have_elements INIT(0);	/* do we have the elements ? */
EXTERN int **na_typ INIT(NULL);	/* number of atoms per type */
EXTERN int nconf INIT(0);	/* number of configurations */
EXTERN int paircol INIT(0);	/* How manc columns for pair pot. */
EXTERN real *maxchange INIT(NULL);	/* Maximal permissible change */
EXTERN atom_t *atoms INIT(NULL);	/* atoms array */
EXTERN real *force_0 INIT(NULL);	/* the forces we aim at */
EXTERN real *coheng INIT(NULL);	/* Cohesive energy for each config */
EXTERN real *volumen INIT(NULL);	/* Volume of cell */
EXTERN stens *stress INIT(NULL);	/* Stresses in each config */
EXTERN int *inconf INIT(NULL);	/* Nr. of atoms in each config */
EXTERN int *cnfstart INIT(NULL);	/* Nr. of first atom in config */
EXTERN int *useforce INIT(NULL);	/* Should we use force/stress */
EXTERN int *usestress INIT(NULL);	/* Should we use force/stress */
EXTERN char startpot[255];	/* file with start potential */
EXTERN char maxchfile[255];	/* file with maximal changes */
EXTERN char endpot[255];	/* file for end potential */
EXTERN char imdpot[255];	/* file for IMD potential */
EXTERN char config[255];	/* file with atom configuration */
EXTERN char plotfile[255];	/* file for plotting */
EXTERN char distfile[255];	/* file for distributions */
EXTERN int write_output_files INIT(0);
EXTERN char output_prefix[255] INIT("");	/* prefix for all output files */
EXTERN char **config_name INIT(NULL);
EXTERN int config_name_max INIT(0);

EXTERN char flagfile[255] INIT("potfit.break");
					 /* break if file exists */
EXTERN char tempfile[255] INIT("\0");	/* backup potential file */
EXTERN char plotpointfile[255] INIT("\0");
					/* write points for plotting */
EXTERN int imdpotsteps;		/* resolution of IMD potential */
EXTERN pot_table_t opt_pot;	/* potential in the internal */
					 /* representation used for  */
					 /* minimisation */
EXTERN pot_table_t calc_pot;	/* the potential table used */
					 /* for force calculations */
#ifdef APOT
EXTERN apot_table_t apot_table;	/* potential in analytic form */
EXTERN int ***pot_list INIT(NULL);	/* list for pairs in potential */
EXTERN int *pot_list_length INIT(NULL);	/* length of pot_list */
EXTERN real plotmin INIT(0.);	/* minimum for plotfile */
EXTERN real *calc_list INIT(NULL);	/* list of current potential in the calc table */
EXTERN int disable_cp INIT(0);	/* switch chemical potential on/off */
EXTERN int compnodes INIT(0);	/* how many additional composition nodes */
EXTERN real *compnodelist INIT(NULL);	/* list of the composition nodes */
EXTERN int cp_start INIT(0);	/* cp in opt_pot.table */
EXTERN int *pot_index INIT(NULL);	/* index to access i*n+j from i*(i+1)/2 */
#endif
EXTERN int format;		/* format of potential table */
EXTERN int opt INIT(0);		/* optimization flag */
EXTERN int plot INIT(0);	/* plot output flag */
EXTERN vektor box_x, box_y, box_z;
EXTERN vektor tbox_x, tbox_y, tbox_z;
EXTERN real rcutmax INIT(0.);	/* maximum of all cutoff values */
EXTERN real *rcut INIT(NULL);
EXTERN real *rmin INIT(NULL);
EXTERN real *lambda INIT(NULL);	/* embedding energy slope... */
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
EXTERN int *idx INIT(NULL);
EXTERN int init_done INIT(0);
/* EXTERN real   *dummy_phi INIT(NULL);     /\* Dummy Constraints for PairPot *\/ */
/* EXTERN real   dummy_rho INIT(1.);        /\* Dummy Constraint for rho *\/ */
/* EXTERN real   dummy_r  INIT(2.5);        /\* Distance of Dummy Constraints *\/ */


/******************************************************************************
*
*  function prototypes
*
******************************************************************************/

void  error(char *);
void  warning(char *);
void  read_parameters(int, char **);
void  read_paramfile(FILE *);
#ifdef APOT
void  read_pot_table(pot_table_t *, apot_table_t *, char *, int);
void  read_apot_table(pot_table_t *pt, apot_table_t *apt, char *filename,
		      FILE *infile);
void  write_apot_table(apot_table_t *, char *);
#else
void  read_pot_table(pot_table_t *, char *, int);
#endif
void  read_pot_table3(pot_table_t *pt, int size, int ncols, int *nvals,
		      char *filename, FILE *infile);
void  read_pot_table4(pot_table_t *pt, int size, int ncols, int *nvals,
		      char *filename, FILE *infile);
void  read_pot_table5(pot_table_t *pt, int size, int ncols, int *nvals,
		      char *filename, FILE *infile);
void  init_calc_table(pot_table_t *optt, pot_table_t *calct);
void  update_calc_table(real *xi_opt, real *xi_calc, int);
void  write_pot_table3(pot_table_t *, char *);
void  write_pot_table4(pot_table_t *, char *);
void  write_pot_table5(pot_table_t *, char *);
void  write_pot_table_imd(pot_table_t *, char *);
void  write_plotpot_pair(pot_table_t *, char *);
void  write_altplot_pair(pot_table_t *, char *);
real  normdist(void);
real  grad2(pot_table_t *, real *, int, real);
real  grad3(pot_table_t *, real *, int, real);
real  pot2(pot_table_t *, int, real);
real  pot3(pot_table_t *, int, real);
void  read_config(char *);
void  read_config2(char *);
real  calc_forces_pair(real *, real *, int);
void  powell_lsq(real *xi);
void  anneal(real *xi);
void  spline_ed(real xstep, real y[], int n, real yp1, real ypn, real y2[]);
real  splint_ed(pot_table_t *pt, real *xi, int col, real r);
real  splint_grad_ed(pot_table_t *pt, real *xi, int col, real r);
real  splint_comb_ed(pot_table_t *pt, real *xi, int col, real r, real *grad);
real  splint_dir(pot_table_t *pt, real *xi, int col, int k, real b,
		 real step);
real  splint_comb_dir(pot_table_t *pt, real *xi, int col, int k, real b,
		      real step, real *grad);
real  splint_grad_dir(pot_table_t *pt, real *xi, int col, int k, real b,
		      real step);
/* real splint_dir_ed(pot_table_t *pt, real *xi, int col, int k, real b); */
/* real splint_grad_dir_ed(pot_table_t *pt, real *xi, int col, int k, real b); */
/* real splint_comb_dir_ed(pot_table_t *pt, real *xi, int col, int k, real b, real *grad); */
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
void  init_mpi(int *argc_pointer, char **argv);
void  shutdown_mpi(void);
void  broadcast_params(void);
void  dbb(int i);
void  potsync();
#endif
#ifdef PDIST
void  write_pairdist(pot_table_t *pt, char *filename);
#endif

#ifdef APOT

#define APOT_STEPS 200		/* number of sampling points for analytic pot */
#define CUTOFF_MARGIN 1.4	/* number to multiply cutoff radius with */
#define APOT_PUNISH 10e10	/* punishment for out of bounds */

int   apot_parameters(char *);
int   apot_assign_functions(apot_table_t *);
int   apot_validate(int, real);
void  new_slots(int, int);	/* new slots for smooth cutoff */
real  chemical_potential(int, int *, real *);
void  init_chemical_potential(int);

#ifdef DEBUG

void  debug_apot();

#endif

/* global counting variables for misc. purposes */

//int count_1 INIT(0);
//int count_2 INIT(0);

/* end of counting variables */

real  smooth(void (*function) (real, real *, real *), real, real *, real *);

#ifdef MPI
void  potsync_apot();
#endif

/* actual functions for different potentials */

void  lj_value(real, real *, real *);
void  eopp_value(real, real *, real *);
void  morse_value(real, real *, real *);
void  softshell_value(real, real *, real *);
void  eopp_exp_value(real, real *, real *);
void  meopp_value(real, real *, real *);
void  power_decay_value(real, real *, real *);
void  pohlong_value(real, real *, real *);
void  parabola_value(real, real *, real *);
void  csw_value(real, real *, real *);
void  universal_value(real, real *, real *);

/* template for new potential function called newpot */

/* newpot potential */
void  newpot_value(real, real *, real *);

/* end of template */

#endif
