
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>



/******************************************************************************
*
*  type definitions
*
******************************************************************************/

typedef double real;
typedef struct { real x; real y; real z; } vektor;

typedef struct {
  int    typ;
  real   r;
  vektor dist;
} neigh_t;

typedef struct {
  int    typ;
  int    n_neigh;
  vektor pos;
  neigh_t *neigh;
} atom_t;

typedef struct {
  real *begin;      /* first value in the table */
  real *end;        /* last value in the table */
  real *step;       /* table increment */
  real *invstep;    /* inverse of increment */
  int  *first;      /* index of first entry */
  int  *last;       /* index of last entry */
  int  len;         /* total length of the table */
  int  ncols;       /* number of columns */
  real *table;      /* the actual data */
} pot_table_t;


#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
static real sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


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
EXTERN int    fcalls   INIT(0);
EXTERN int    ndim;
EXTERN int    mdim;
EXTERN int    ntypes   INIT(1);          /* number of atom types */
EXTERN int    natoms   INIT(0);          /* number of atoms */
EXTERN atom_t *atoms   INIT(NULL);       /* atoms array */
EXTERN real   *force_0 INIT(NULL);       /* the forces we aim at */
EXTERN char startpot[255];               /* file with start potential */
EXTERN char endpot[255];                 /* file for end potential */
EXTERN char imdpot[255];                 /* file for IMD potential */
EXTERN char config[255];                 /* file with atom configuration */
EXTERN int  imdpotsteps;                 /* resolution of IMD potential */
EXTERN pot_table_t pair_pot;             /* the potential table */

/******************************************************************************
*
*  function prototypes
*
******************************************************************************/

void error(char*);
void read_parameters(int, char**);
void read_paramfile(FILE*);
void read_pot_table(pot_table_t*, char*, int);
void write_pot_table(pot_table_t*, char*);
void write_pot_table_imd(pot_table_t*, char*);
real grad2(pot_table_t*, real*, int, real);
real grad3(pot_table_t*, real*, int, real);
real pot2 (pot_table_t*, int, real);
real pot3 (pot_table_t*, int, real);
void read_config(char*);
real calc_forces_pair(real[],real[]);
void powell_lsq(real *xi, real (*fcalc)(real[],real[]));







