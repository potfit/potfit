
#include "potfit.h"

/*****************************************************************************
*
*  compute forces using pair potentials 
*
******************************************************************************/

real calc_forces_pair(real *xi, real *forces)
{
  int     i, j, k, typ1, typ2, col;
  atom_t  *atom;
  neigh_t *neigh;
  real    grad, sum=0.0;

  for (i=0; i<natoms; i++) {

    atom = atoms + i;
    typ1 = atom->typ;
    k    = 3*i;
    forces[k  ] = 0.0;
    forces[k+1] = 0.0;
    forces[k+2] = 0.0;

    for (j=0; j<atom->n_neigh; j++) {

      neigh = atom->neigh+j;
      typ2  = neigh->typ;
      col   = (typ1 <= typ2) ? typ1 * ntypes + typ2 : typ2 * ntypes + typ1;

      if (neigh->r < pair_pot.end[col] + pair_pot.step[col]) {
        grad = grad3( &pair_pot, xi, col, neigh->r);
        forces[k  ] += neigh->dist.x * grad;
        forces[k+1] += neigh->dist.y * grad;
        forces[k+2] += neigh->dist.z * grad;
      }
    }
  }
  fcalls++;			/* Increase function call counter */
  for (i=0; i<3*natoms; i++) sum += SQR(forces[i]-force_0[i]);
  return sum;
}
