/****************************************************************
* 
* force.c: Routines used for calculating forces/energies in various 
*     interpolation schemes. 
*
*****************************************************************/
/****************************************************************
* $Revision: 1.15 $
* $Date: 2003/03/19 10:12:16 $
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
  int     i, j, k, typ1, typ2, col,first;
  int config;
  atom_t  *atom;
  neigh_t *neigh;
  real    pot,grad, y0, y1, x0, x1, sum=0.0;
  
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
  for (i=0; i<nconf;i++) forces[3*natoms+i] = 0;
  for (i=0; i<natoms; i++) {

    atom = atoms + i;
    typ1 = atom->typ;
    config=3*natoms+atom->conf;
    k    = 3*i;
    forces[k  ] = -force_0[k  ];
    forces[k+1] = -force_0[k+1];
    forces[k+2] = -force_0[k+2];

    for (j=0; j<atom->n_neigh; j++) {
      
      neigh = atom->neigh+j;
      typ2  = neigh->typ;
      col   = (typ1 <= typ2) ? typ1 * ntypes + typ2 - ((typ1 * (typ1 + 1))/2) 
			     : typ2 * ntypes + typ1 - ((typ2 * (typ2 + 1))/2);

      if (neigh->r < pair_pot.end[col] + pair_pot.step[col]) {
	/* not a real force: cohesive energy */
        forces[config] += splint_ed(&pair_pot,xi,col,neigh->r);
	grad = splint_grad_ed(&pair_pot, xi, col, neigh->r);
        forces[k  ] += neigh->dist.x * grad;
        forces[k+1] += neigh->dist.y * grad;
        forces[k+2] += neigh->dist.z * grad;
      }
    }
    /* Returned force is difference between calculated and input force */
    sum += SQR(forces[k]) + SQR(forces[k+1]) + SQR(forces[k+2]);
  }
  for (i=0;i<nconf;i++) {
    forces[3*natoms+i]/=(real) inconf[i]*2.0; /* double counting... */
    forces[3*natoms+i]-=force_0[3*natoms+i];
    sum += SQR(forces[3*natoms+i]);
}
  fcalls++;			/* Increase function call counter */
  return sum;
}

