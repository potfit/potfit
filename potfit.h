/****************************************************************
* 
*  potfit.h: potfit header file
*
*****************************************************************/


/****************************************************************
* $Revision: 1.39 $
* $Date: 2004/11/17 17:28:49 $
*****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
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
#define MAXNEIGH 170
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
typedef struct { real x; real y; real z; } vektor;
typedef struct { real xx; real yy; real zz; real xy; real yz; real zx; } stens;

typedef struct {
  int    typ;
  int    nr;
  real   r;
  vektor dist;
  int    slot[2];
  real   shift[2];
  real   step[2];
} neigh_t;

typedef struct {
  int    typ;
  int    n_neigh;
  vektor pos;
  vektor force;
  real   absforce;
  neigh_t neigh[MAXNEIGH];
  int    conf; 			/* Which configurarion... */
#ifdef EAM
  real   rho;			/* embedding electron density */
  real   gradF;			/* gradient of embedding fn. */
#endif
} atom_t;

typedef struct {
  int  len;         /* total length of the table */
  int  idxlen;      /* number of changeable potential values */
  int  ncols;       /* number of columns */
  real *begin;      /* first value in the table */
  real *end;        /* last value in the table */
  real *step;       /* table increment */
  real *invstep;    /* inverse of increment */
  int  *first;      /* index of first entry */
  int  *last;       /* index of last entry */
  real *xcoord;     /* the x-coordinates of sampling points */
  real *table;      /* the actual data */
  real *d2tab;      /* second derivatives of table data for spline int */
  int  *idx;        /* indirect indexing */
} pot_table_t;

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
#define EXTERN             /* define Variables in main */
#define INIT(data) =data   /* initialize data only in main */
#else
#define EXTERN extern      /* declare them extern otherwise */
#define INIT(data)         /* skip initialization otherwise */
#endif
EXTERN int myid INIT(0);                  /* Who am I? (0 if serial) */
EXTERN int num_cpus INIT(1);              /* How many cpus are there */
#ifdef MPI
EXTERN MPI_Datatype MPI_VEKTOR;
EXTERN MPI_Datatype MPI_STENS;
//EXTERN MPI_Datatype MPI_POTTABLE;
EXTERN MPI_Datatype MPI_ATOM;
EXTERN MPI_Datatype MPI_NEIGH;
#endif
EXTERN atom_t *conf_atoms INIT(NULL); /* Atoms in configuration */
EXTERN real *conf_eng INIT(NULL);
EXTERN real *conf_vol INIT(NULL);
EXTERN stens *conf_stress INIT(NULL);
EXTERN int *atom_len;
EXTERN int *atom_dist;
EXTERN int *conf_len;
EXTERN int *conf_dist;
EXTERN real sweight INIT(STRESS_WEIGHT);
EXTERN real eweight INIT(ENG_WEIGHT);
EXTERN int myconf INIT(0.);
EXTERN int myatoms INIT(0.);
EXTERN int firstconf INIT(0);
EXTERN int firstatom INIT(0);
EXTERN real   pi       INIT(0.);
EXTERN real   anneal_temp INIT(1.);
EXTERN int    seed     INIT(123456);     /* seed for RNG */
EXTERN int    fcalls   INIT(0);
EXTERN int    ndim     INIT(0);
EXTERN int    ndimtot  INIT(0);
EXTERN int    mdim     INIT(0);
EXTERN int    ntypes   INIT(1);          /* number of atom types */
EXTERN int    natoms   INIT(0);          /* number of atoms */
EXTERN int    nconf    INIT(0);	         /* number of configurations */
EXTERN int    eam      INIT(0);	         /* EAM usage flag */
EXTERN int    paircol  INIT(0);          /* How manc columns for pair pot.*/
EXTERN atom_t *atoms   INIT(NULL);       /* atoms array */
EXTERN real   *force_0 INIT(NULL);       /* the forces we aim at */
EXTERN real   *coheng  INIT(NULL);       /* Cohesive energy for each config */
EXTERN real   *volumen INIT(NULL);       /* Volume of cell*/
EXTERN stens  *stress  INIT(NULL);       /* Stresses in each config */
EXTERN int    *inconf  INIT(NULL);       /* Nr. of atoms in each config */
EXTERN int    *cnfstart INIT(NULL);       /* Nr. of first atom in config */
EXTERN char startpot[255];               /* file with start potential */
EXTERN char endpot[255];                 /* file for end potential */
EXTERN char imdpot[255];                 /* file for IMD potential */
EXTERN char config[255];                 /* file with atom configuration */
EXTERN char plotfile[255];               /* file for plotting */
EXTERN char distfile[255];	         /* file for distributions*/
EXTERN char flagfile[255] INIT("potfit.break");
				         /* break if file exists */
EXTERN  char tempfile[255] INIT("\0");   /* backup potential file */
EXTERN char plotpointfile[255] INIT("\0");
                                         /* write points for plotting */
EXTERN int  imdpotsteps;                 /* resolution of IMD potential */
EXTERN pot_table_t pair_pot;             /* the potential table */
EXTERN int format;		         /* format of potential table */
EXTERN int  opt INIT(0);                 /* optimization flag */
EXTERN int plot INIT(0);                 /* plot output flag */
EXTERN vektor box_x,  box_y,  box_z;
EXTERN vektor tbox_x, tbox_y, tbox_z;
EXTERN real *rcut INIT(NULL);
EXTERN real *rmin INIT(NULL);
EXTERN real *lambda INIT(NULL);	/* embedding energy slope... */
EXTERN real (*calc_forces)(real*,real*,int);
EXTERN real (*splint)(pot_table_t*,real*,int,real);
EXTERN real (*splint_grad)(pot_table_t*,real*,int,real);
EXTERN real (*splint_comb)(pot_table_t*,real*,int,real,real*);
EXTERN void (*write_pot_table)(pot_table_t*,char*);
#ifdef PARABEL
EXTERN real (*parab)(pot_table_t*,real*,int,real);
EXTERN real (*parab_comb)(pot_table_t*,real*,int,real,real*);
EXTERN real (*parab_grad)(pot_table_t*,real*,int,real);
#endif /* PARABEL */
EXTERN int  *idx INIT(NULL);
/* EXTERN real   *dummy_phi INIT(NULL);     /\* Dummy Constraints for PairPot *\/ */
/* EXTERN real   dummy_rho INIT(1.);        /\* Dummy Constraint for rho *\/ */
/* EXTERN real   dummy_r  INIT(2.5);        /\* Distance of Dummy Constraints *\/ */


/******************************************************************************
*
*  function prototypes
*
******************************************************************************/

void error(char*);
void warning(char*);
void read_parameters(int, char**);
void read_paramfile(FILE*);
void read_pot_table(pot_table_t*, char*, int);
void read_pot_table3(pot_table_t *pt, int size, int ncols, int *nvals, 
		     char *filename, FILE *infile);
void read_pot_table4(pot_table_t *pt, int size, int ncols, int *nvals, 
		     char *filename, FILE *infile);
void write_pot_table3(pot_table_t*, char*);
void write_pot_table4(pot_table_t*, char*);
void write_pot_table_imd(pot_table_t*, char*);
void write_plotpot_pair(pot_table_t*, char*);
void write_altplot_pair(pot_table_t*, char*);
real normdist(void);
real grad2(pot_table_t*, real*, int, real);
real grad3(pot_table_t*, real*, int, real);
real pot2 (pot_table_t*, int, real);
real pot3 (pot_table_t*, int, real);
void read_config(char*);
void read_config2(char*);
real calc_forces_pair(real*,real*,int);
void powell_lsq(real *xi);
void anneal(real *xi);
void spline_ed(real xstep, real y[], int n, real yp1, real ypn, real y2[]);
real splint_ed(pot_table_t *pt, real *xi,int col, real r);
real splint_grad_ed(pot_table_t *pt, real *xi, int col, real r);
real splint_comb_ed(pot_table_t *pt, real *xi, int col, real r, real *grad);
real splint_dir(pot_table_t *pt, real *xi, int col, int k, real b, real step);
real splint_comb_dir(pot_table_t *pt, real *xi, int col, int k, real b, 
		     real step, real *grad);
real splint_grad_dir(pot_table_t *pt, real *xi, int col, int k, real b, 
		     real step);
/* real splint_dir_ed(pot_table_t *pt, real *xi, int col, int k, real b); */
/* real splint_grad_dir_ed(pot_table_t *pt, real *xi, int col, int k, real b); */
/* real splint_comb_dir_ed(pot_table_t *pt, real *xi, int col, int k, real b, real *grad); */
void spline_ne(real x[], real y[], int n, real yp1, real ypn, real y2[]);
real splint_ne(pot_table_t *pt, real *xi, int col, real r);
real splint_comb_ne(pot_table_t *pt, real *xi, int col, real r, real *grad);
real splint_grad_ne(pot_table_t *pt, real *xi, int col, real r);


#ifdef PARABEL
real parab_comb_ed(pot_table_t *pt,  real *xi, int col, real r, real *grad);
real parab_grad_ed(pot_table_t *pt,  real *xi, int col, real r);
real parab_ed(pot_table_t *pt,  real *xi, int col, real r);
real parab_comb_ne(pot_table_t *pt,  real *xi, int col, real r, real *grad);
real parab_grad_ne(pot_table_t *pt,  real *xi, int col, real r);
real parab_ne(pot_table_t *pt,  real *xi, int col, real r);
#endif
#ifdef EAM
real rescale(pot_table_t *pt, real upper,int flag);
void embed_shift(pot_table_t *pt);
#endif
#ifdef MPI
void init_mpi(int *argc_pointer, char **argv);
void shutdown_mpi(void);
void broadcast_params(void);
void dbb(int i);
void potsync();
#endif
#ifdef PDIST
void write_pairdist(pot_table_t *pt, char *filename);
#endif
