/**** rewritten for rdouble precision and zero-offset vectors and matrices ****
***** by Peter Brommer, ITAP, 2002-10-10                                  ****
***** adapted to Powell requirements (return vector ...) 2002-10-11       ***/
/**** adapted to real variables (ITAP standard) by PB, ITAP, 2002-10-24   ***/


#define NRANSI
#include "nrutil_r.h"

extern int ncom, mcom;
extern real *pcom,*xicom;
extern real (*nrfunc)(real [], real [], int, int);

real f1dim_r(real x)
{
	int j;
	real f,*xt;
	real *res;

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
