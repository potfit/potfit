
#include "potfit.h"

/*****************************************************************************
*
*  compute forces using pair potentials with polynomial interpolation
*
******************************************************************************/

real calc_forces_pair_poly(real *xi, real *forces)
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
  int     i, j, k, typ1, typ2, col;
  atom_t  *atom;
  neigh_t *neigh;
  real    grad, sum=0.0;
  
  /* init second derivatives for splines */

  for (col=0; col<pair_pot.ncols; col++) 
      spline_ed(pair_pot.step[col],xi+pair_pot.first[col],
		pair_pot.last[col]-pair_pot.first[col]+1,
		1e30,0,pair_pot.d2tab+pair_pot.first[col]);

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

      if (neigh->r < pair_pot.end[col] + pair_pot.step[col]) {
        grad = splint_grad_ed(&pair_pot, xi, col, neigh->r);
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

