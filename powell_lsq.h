double brent_d(double ax, double bx, double cx, double fbx, double tol, 
	double *xmin, double *xmin2, double *fxmin, double *fxmin2);

void copy_matrix(double **a, double **b, int n, int m);

void copy_vector(double *a, double *b, int n);

double f1dim_d(double x);

double force_calc(double *xi, double *force, int n, int m);

/*
void fxi_init(double *xi, double *fsoll, double *force_xi, int n, int m);
*/

void gamma_init(double **gamma, double **d, double *xi, double *force_xi, int n,
		int m);

int gamma_update(double **gamma, double a, double b, double *fa, double *fb,
		int j, int n, int m);
		
void lineqsys_init(double **gamma, double **lineqsys, double *deltaforce, 
		double *p, int n, int m);

void lineqsys_update(double **gamma, double **lineqsys, double *force_xi,
		double *p, int i, int n, int m);
		
double linmin_d(double p[], double xi[], double fxi1, int n, int m,
	double *x1, double *x2, double *fret1, double *fret2, 
	double (*func)(double [], double [], int, int));

void lubksb_d(double **a, int n, int *indx, double b[]);

void ludcmp_d(double **a, int n, int *indx, double *d);

void matdotvec(double **a, double *x, double *y, int n, int m);

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, 
	double *fc, double (*func)(double));

void mprove_d(double **a, double **alud, int n, int indx[], double b[], 
		double x[]);

double normalize_vector(double *v, int n);
