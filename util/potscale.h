/******************************************************************
 *
 *  potscale.h: potscale header file
*
*****************************************************************/

/****************************************************************
* $Revision: 1.1 $
* $Date: 2004/12/03 17:31:41 $
*****************************************************************/

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

/******************************************************************************
*
*  global variables
*
******************************************************************************/


#ifdef MAIN 
#define EXTERN             /* define Variables in main */
#define INIT(data) =data   /* initialize data only in main */
#else
#define EXTERN extern      /* declare them extern otherwise */
#define INIT(data)         /* skip initialization otherwise */
#endif

EXTERN real   extend INIT(2.);	/* how far should one extend imd pot */
EXTERN int    ndim     INIT(0);
EXTERN int    ndimtot  INIT(0);
EXTERN int    ntypes   INIT(1);          /* number of atom types */

EXTERN int    paircol  INIT(0);          /* How manc columns for pair pot.*/
EXTERN char startpot[255];               /* file with start potential */
EXTERN char imdpot[255];                 /* file for IMD potential */
EXTERN char plotfile[255];               /* file for plotting */
EXTERN int  imdpotsteps;                 /* resolution of IMD potential */
EXTERN pot_table_t pair_pot;             /* the potential table */
EXTERN int format;		         /* format of potential table */

EXTERN real *lambda INIT(NULL);	/* embedding energy slope... */
EXTERN real (*splint)(pot_table_t*,real*,int,real);
EXTERN real (*splint_grad)(pot_table_t*,real*,int,real);
EXTERN real (*splint_comb)(pot_table_t*,real*,int,real,real*);

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
void write_pot_table_imd(pot_table_t*, char*);
void write_plotpot_pair(pot_table_t*, char*);
void write_altplot_pair(pot_table_t*, char*);
void spline_ed(real xstep, real y[], int n, real yp1, real ypn, real y2[]);
real splint_ed(pot_table_t *pt, real *xi,int col, real r);
real splint_grad_ed(pot_table_t *pt, real *xi, int col, real r);
real splint_comb_ed(pot_table_t *pt, real *xi, int col, real r, real *grad);
real splint_dir(pot_table_t *pt, real *xi, int col, int k, real b, real step);
real splint_comb_dir(pot_table_t *pt, real *xi, int col, int k, real b, 
		     real step, real *grad);
real splint_grad_dir(pot_table_t *pt, real *xi, int col, int k, real b, 
		     real step);
void spline_ne(real x[], real y[], int n, real yp1, real ypn, real y2[]);
real splint_ne(pot_table_t *pt, real *xi, int col, real r);
real splint_comb_ne(pot_table_t *pt, real *xi, int col, real r, real *grad);
real splint_grad_ne(pot_table_t *pt, real *xi, int col, real r);

