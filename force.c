/****************************************************************
* 
* force.c: Routines used for calculating forces/energies in various 
*     interpolation schemes. 
*
*****************************************************************/
/****************************************************************
* $Revision: 1.18 $
* $Date: 2003/04/17 13:59:26 $
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
  int paircol=(ntypes*(ntypes+1))/2;
  real grad,y0,y1,x0,x1;
  /* init second derivatives for splines */
  for (col=0; col<paircol; col++){  /* just pair potentials */
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
#ifdef EAM
  if (eam) {
    for (col=paircol; col<paircol+ntypes; col++) { /* rho */
      first=pair_pot.first[col];
      spline_ed(pair_pot.step[col], xi+first, 
		pair_pot.last[col]-first+1,
		1e30,0.0,pair_pot.d2tab+first);
    }
    for  (col=paircol+ntypes; col<paircol+2*ntypes; col++) { /* F */
      first=pair_pot.first[col];
      spline_ed(pair_pot.step[col], xi+first, 
		pair_pot.last[col]-first+1,
		1e30,0.,pair_pot.d2tab+first);
    }
  }

#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    int h,k,i,j,typ1,typ2,col,config;
#ifdef EAM
    int col2;
    real gradF;
#endif
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
#ifdef EAM
        atom->rho=0.0;
#endif
	for (j=0; j<atom->n_neigh; j++) {
	
	  neigh = atom->neigh+j;
	  typ2  = neigh->typ;
	  col = (typ1 <= typ2) ? typ1 * ntypes + typ2 - ((typ1 * (typ1 + 1))/2)
	    : typ2 * ntypes + typ1 - ((typ2 * (typ2 + 1))/2);
#ifdef EAM
	  col2 = paircol+typ2;
#endif EAM
	  if (neigh->r < pair_pot.end[col]) {
	    /* not a real force: cohesive energy */
	    /* grad is calculated in the same step */
	    fnval = splint_comb_ed(&pair_pot,xi,col,neigh->r,&grad);
	    forces[config] += fnval;
	    forces[k  ] += neigh->dist.x * grad;
	    forces[k+1] += neigh->dist.y * grad;
	    forces[k+2] += neigh->dist.z * grad;
	  }
#ifdef EAM
	  if (eam && neigh->r < pair_pot.end[col2]) {
	    atom->rho += splint_ed(&pair_pot,xi,col2,neigh->r);
	  }	
#endif
	}
#ifndef EAM /*we can't do that here right now*/
	/* Returned force is difference between calculated and input force */
	sum += SQR(forces[k]) + SQR(forces[k+1]) + SQR(forces[k+2]);
#else
	if (!eam) sum += SQR(forces[k]) + SQR(forces[k+1]) + SQR(forces[k+2]);
#endif
      }	/* first loop over atoms */
#ifdef EAM
      if (eam) {
	for (i=0; i<inconf[h]; i++) {
	  atom = atoms + i + cnfstart[h];
	  typ1 = atom->typ;
	  k    = 3*(cnfstart[h]+i);
	  col  = paircol+ntypes+typ1; /* column of F */
	  if (atom->rho > pair_pot.end[col])
	    atom->rho = pair_pot.end[col]; /* rho too big... bad thing */
	  if (atom->rho < pair_pot.begin[col]) {
	    forces[config] += 1e10; /* rho too small... really bad thing*/
	    atom->rho = pair_pot.begin[col];
	  }
	  fnval=splint_comb_ed(&pair_pot,xi,col,atom->rho,&gradF);
	  forces[config]+= 2 * fnval;
	  for (j=0; j<atom->n_neigh; j++) {
	    neigh = atom->neigh+j;
	    typ2  = neigh->typ;
	    col2  = paircol+typ2;
	    if (neigh->r < pair_pot.end[col2]) {
	      grad=splint_grad_ed(&pair_pot,xi,col2,neigh->r)*gradF;
	      forces[k  ] += neigh->dist.x * grad;
	      forces[k+1] += neigh->dist.y * grad;
	      forces[k+2] += neigh->dist.z * grad;
	    }
	  } /* loop over neighbours */
	sum += SQR(forces[k]) + SQR(forces[k+1]) + SQR(forces[k+2]);
	}	/* loop over atoms */
      }	/* if eam is used */
#endif
      forces[config]/=(real) inconf[h]*2.0; /* double counting... */
      forces[config]-=force_0[config];
      sum += SQR(forces[config]);
    } /* loop over configurations */
  } /* parallel region */
#ifdef EAM
  if (eam) {
    forces[mdim-2]=
      splint_ed(&pair_pot,xi,paircol+DUMMY_COL_RHO,DUMMY_R_RHO)
      - force_0[mdim-2];
    forces[mdim-1]=splint_ed(&pair_pot,xi,DUMMY_COL_PHI,DUMMY_R_PHI)
      - force_0[mdim-1];
    sum+= SQR(forces[mdim-2]);
    sum+= SQR(forces[mdim-1]);
  }
#endif 
  fcalls++;			/* Increase function call counter */
  return sum;
}

