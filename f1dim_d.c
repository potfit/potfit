/**** rewritten for double precision and zero-offset vectors and matrices ****
***** by Peter Brommer, ITAP, 2002-10-10                                  ****
***** adapted to Powell requirements (return vector ...) 2002-10-11       ***/


#define NRANSI
#include "nrutil.h"

extern int ncom, mcom;
extern double *pcom,*xicom;
extern double (*nrfunc)(double [], double [], int, int);

double f1dim_d(double x)
{
	int j;
	double f,*xt;
	double *res;

	xt=dvector(0,ncom-1);
	res=dvector(0,mcom-1);

	for (j=0;j<ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt,res,ncom, mcom);
	free_dvector(xt,0,ncom-1);
	free_dvector(res,0,mcom-1);
	return f;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */
