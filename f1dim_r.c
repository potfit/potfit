/****************************************************************
* 
*  f1dim_r.c: One-dimensional replacement function. 
*
*****************************************************************/
/****************************************************************
* $Revision: 1.5 $
* $Date: 2003/03/19 09:05:34 $
*****************************************************************/
/**** rewritten for rdouble precision and zero-offset vectors and matrices ****
***** by Peter Brommer, ITAP, 2002-10-10                                  ****
***** adapted to Powell requirements (return vector ...) 2002-10-11       ***/
/**** adapted to real variables (ITAP standard) by PB, ITAP, 2002-10-24   ***/

#include "potfit.h"
#include "nrutil_r.h"

extern real *xicom, *delcom;

real f1dim_r(real x)
{
	int j;
	static real *res = NULL, *xt = NULL;

	if (xt  == NULL) xt  = dvector(0,ndimtot-1);
	if (res == NULL) res = dvector(0,mdim-1);

	for (j=0; j<ndimtot; j++) xt[j]=xicom[j]+x*delcom[j];
	return (*calc_forces)(xt,res);

}

