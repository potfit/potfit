/****************************************************************
* 
* force.c: Routines used for calculating forces/energies in various 
*     interpolation schemes. 
*
*****************************************************************/
/****************************************************************
* $Revision: 1.17 $
* $Date: 2003/04/10 12:43:35 $
*****************************************************************/

#include "potfit.h"

/*****************************************************************************
*
*  compute forces using pair potentials with polynomial interpolation
*
******************************************************************************/

/*********************OBSOLETE*************************/

real calc_forces_pair_poly(real *xi, real *forces)
/******************************************************/
{
  int     i, j, k, typ1, typ2, col;
  atom_t  *atom;
  neigh_t *neigh;
  real    grad, sum=0.0;
  
  for (i=0; i<natoms; i++) {

    atom = atoms + i;
    typ1 = atom->typ;
    k    = 3*i;
    forces[k  ] = -force_0[k  ];
    forces[k+1] = -force_0[k+1];
    forces[k+2] = -force_0[k+2];

    for (j=0; j<atom->n_neigh; j++) {

      neigh = atom->neigh+j;
      typ2  = neigh->typ;
      col   = (typ1 <= typ2) ? typ1 * ntypes + typ2 - ((typ1 * (typ1 + 1))/2) 
			     : typ2 * ntypes + typ1 - ((typ2 * (typ2 + 1))/2);

      if (neigh->r < pair_pot.end[col]) {
        grad = grad3( &pair_pot, xi, col, neigh->r);
        forces[k  ] += neigh->dist.x * grad;
        forces[k+1] += neigh->dist.y * grad;
        forces[k+2] += neigh->dist.z * grad;
      }
    }
    /* Returned force is difference between calculated and input force */
    sum += SQR(forces[k]) + SQR(forces[k+1]) + SQR(forces[k+2]);
  }
  fcalls++;			/* Increase function call counter */
  return sum;
}

/*****************************************************************************
*
*  compute forces using pair potentials with spline interpolation
*
******************************************************************************/

real calc_forces_pair(real *xi, real *forces)
{
  real sum=0.0;
  int first,col;
  real grad,y0,y1,x0,x1;
  /* init second derivatives for splines */
  for (col=0; col<pair_pot.ncols; col++){
    first=pair_pot.first[col];
    x0=pair_pot.begin[col];
    x1=x0+pair_pot.step[col];
    y0=xi[first];
    y1=xi[first+1];
    if (y0*y1>0)
      grad=(y0*log(y0/y1)) /  (x0*log(x0/x1));
    else
      grad=1e30;
    if (!((grad>-1e10) && (grad<1e10))) grad=1e30;
    if (  (grad>-1e-20)&& (grad<1e-20)) grad=0.;
    spline_ed(pair_pot.step[col], xi+first,
	      pair_pot.last[col]-first+1,
	      grad, 0.0, pair_pot.d2tab+first);
  }
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    int h,k,i,j,typ1,typ2,col,config;
    real fnval, grad;
    atom_t *atom;
    neigh_t *neigh;
#ifdef _OPENMP
#pragma omp for reduction(+:sum)
#endif 
    for (h=0; h<nconf; h++) {
      config = 3*natoms + h;
      forces[config]=0.;
      for (i=0; i<inconf[h]; i++) {
	atom = atoms + i + cnfstart[h];
	typ1 = atom->typ;
	k    = 3*(cnfstart[h]+i);
	forces[k  ] = -force_0[k  ];
	forces[k+1] = -force_0[k+1];
	forces[k+2] = -force_0[k+2];
      
	for (j=0; j<atom->n_neigh; j++) {
	
	  neigh = atom->neigh+j;
	  typ2  = neigh->typ;
	  col = (typ1 <= typ2) ? typ1 * ntypes + typ2 - ((typ1 * (typ1 + 1))/2)
	    : typ2 * ntypes + typ1 - ((typ2 * (typ2 + 1))/2);

	  if (neigh->r < pair_pot.end[col] + pair_pot.step[col]) {
	    /* not a real force: cohesive energy */
	    /* grad is calculated in the same step */
	    fnval = splint_comb_ed(&pair_pot,xi,col,neigh->r,&grad);
	    forces[config] += fnval;
	    forces[k  ] += neigh->dist.x * grad;
	    forces[k+1] += neigh->dist.y * grad;
	    forces[k+2] += neigh->dist.z * grad;
	  }
	}
	/* Returned force is difference between calculated and input force */
	sum += SQR(forces[k]) + SQR(forces[k+1]) + SQR(forces[k+2]);
      }
      forces[config]/=(real) inconf[h]*2.0; /* double counting... */
      forces[config]-=force_0[config];
      sum += SQR(forces[config]);

    }
  }
  fcalls++;			/* Increase function call counter */
  return sum;
}

