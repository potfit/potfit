/****************************************************************
* 
*  powell_lsq.h: Header file used by powell routines.
* 
*****************************************************************/

/****************************************************************
* $Revision: 1.6 $
* $Date: 2003/03/19 09:05:42 $
*****************************************************************/



#ifdef POWELL
#define EXPOW             /* define Variables in powell */
#define INPOW(data) =data   /* initialize data only in powell */
#else
#define EXPOW extern      /* declare them extern otherwise */
#define INPOW(data)         /* skip initialization otherwise */
#endif

real brent_r(real ax, real bx, real cx, real fbx, real tol, 
	real *xmin, real *xmin2, real *fxmin, real *fxmin2);

void copy_matrix(real **a, real **b, int n, int m);

void copy_vector(real *a, real *b, int n);

real f1dim_r(real x);

int gamma_init(real **gamma, real **d, real *xi, real *force_xi, int n, int m);

int gamma_update(real **gamma, real a, real b, real *fa, real *fb,
		int j, int m);
		
void lineqsys_init(real **gamma, real **lineqsys, real *deltaforce, 
		real *p, int n, int m);

void lineqsys_update(real **gamma, real **lineqsys, real *force_xi,
		real *p, int i, int n, int m);
		
real linmin_r(real p[], real xi[], real fxi1, int n, int m,
	real *x1, real *x2, real *fret1, real *fret2);

void lubksb_r(real **a, int n, int *indx, real b[]);

void ludcmp_r(real **a, int n, int *indx, real *d);

void matdotvec(real **a, real *x, real *y, int n, int m);

void mnbrak(real *ax, real *bx, real *cx, real *fa, real *fb, 
	real *fc, real (*func)(real));

void mprove_r(real **a, real **alud, int n, int indx[], real b[], 
		real x[]);

real normalize_vector(real *v, int n);
