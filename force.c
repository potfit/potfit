/****************************************************************
* 
* force.c: Routines used for calculating forces/energies in various 
*     interpolation schemes. 
*
*****************************************************************/
/****************************************************************
* $Revision: 1.20 $
* $Date: 2003/05/16 12:16:59 $
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
  int first,col1;
  int paircol=(ntypes*(ntypes+1))/2;
  real grad0,y0,y1,x0,x1;
  /* init second derivatives for splines */
  for (col1=0; col1<paircol; col1++){  /* just pair potentials */
    first=pair_pot.first[col1];
    x0=pair_pot.begin[col1];
    x1=x0+pair_pot.step[col1];
    y0=xi[first];
    y1=xi[first+1];
    if (y0*y1>0)
      grad0=(y0*log(y0/y1)) /  (x0*log(x0/x1));
    else
      grad0=1e30;
    if (!((grad0>-1e10) && (grad0<1e10))) grad0=1e30;
    if (  (grad0>-1e-20)&& (grad0<1e-20)) grad0=0.;
    spline_ed(pair_pot.step[col1], xi+first,
	      pair_pot.last[col1]-first+1,
	      grad0, 0.0, pair_pot.d2tab+first);
  }
#ifdef EAM
  if (eam) {
    for (col1=paircol; col1<paircol+ntypes; col1++) { /* rho */
      first=pair_pot.first[col1];
      spline_ed(pair_pot.step[col1], xi+first, 
		pair_pot.last[col1]-first+1,
		1e30,0.0,pair_pot.d2tab+first);
    }
#ifndef PARABEL
    for  (col1=paircol+ntypes; col1<paircol+2*ntypes; col1++) { /* F */
      first=pair_pot.first[col1];
      spline_ed(pair_pot.step[col1], xi+first, 
		pair_pot.last[col1]-first+1,
		0.,0.,pair_pot.d2tab+first);
    }
#endif
  }

#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
#ifdef EAM
    vektor d;
    int col2,l;
    real gradF,gradF2,grad2,r,eamforce;
    atom_t *atom2;
#endif
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
      forces[config+nconf]=force_0[config+nconf];

      for (i=0; i<inconf[h]; i++) {
	k    = 3*(cnfstart[h]+i);
	forces[k  ] = -force_0[k  ];
	forces[k+1] = -force_0[k+1];
	forces[k+2] = -force_0[k+2];
        atoms[cnfstart[h]+i].rho=0.0;
      }
      for (i=0; i<inconf[h]; i++) {
	atom = atoms + i + cnfstart[h];
	typ1 = atom->typ;
	k    = 3*(cnfstart[h]+i);
#ifdef EAM

#endif
	for (j=0; j<atom->n_neigh; j++) {
	
	  neigh = atom->neigh+j;
	  if (neigh->nr > i+cnfstart[h]) {
	    typ2  = neigh->typ;
	    col = (typ1 <= typ2) ? typ1 * ntypes + typ2 - ((typ1 * (typ1 + 1))/2)
	      : typ2 * ntypes + typ1 - ((typ2 * (typ2 + 1))/2);
#ifdef EAM

#endif EAM
	    if (neigh->r < pair_pot.end[col]) {
	      /* not a real force: cohesive energy */
	      /* grad is calculated in the same step */
	      fnval = splint_comb_ed(&pair_pot,xi,col,neigh->r,&grad);
	      forces[config] += fnval;
	      forces[k  ] += neigh->dist.x * grad;
	      forces[k+1] += neigh->dist.y * grad;
	      forces[k+2] += neigh->dist.z * grad;
	      l    = 3*neigh->nr;
	      forces[l  ] -= neigh->dist.x * grad; 
	      forces[l+1] -= neigh->dist.y * grad; 
	      forces[l+2] -= neigh->dist.z * grad; 
	    }
	  
#ifdef EAM
	    if (eam) {
	      col2 = paircol+typ2;
	      if (typ2==typ1) {
		if (neigh->r < pair_pot.end[col2]) {
		  fnval = splint_ed(&pair_pot,xi,col2,neigh->r);
		  atom->rho += fnval;
		  atoms[neigh->nr].rho += fnval;
		}
	      } else {
		col = paircol+typ1;
		if (neigh->r < pair_pot.end[col2])
		  atom->rho += splint_ed(&pair_pot,xi,col2,neigh->r);
		if (neigh->r < pair_pot.end[col])
		  atoms[neigh->nr].rho += splint_ed(&pair_pot,xi,col,neigh->r);
	      }	
	    }
	    
#endif
	  }
	}
#ifndef EAM /*we can't do that here right now*/
	/* Returned force is difference between calculated and input force */
	sum += SQR(forces[k]) + SQR(forces[k+1]) + SQR(forces[k+2]);
#else
	if (!eam) {
	  sum += SQR(forces[k]) + SQR(forces[k+1]) + SQR(forces[k+2]);
	} else {
	  col2=paircol+ntypes+typ1;
	  if (atom->rho > pair_pot.end[col2]) {
#ifdef LIMIT
	    forces[config+nconf]+=1000*SQR(atom->rho - pair_pot.end[col2]);
#endif
	    atom->rho = pair_pot.end[col2];
	  }
/* rho too big... bad thing */
	  if (atom->rho < pair_pot.begin[col2]) {
	    /*  forces[config] += 1e10;  rho too small... really bad thing*/
#ifdef LIMIT
	    forces[config+nconf]+=1000*SQR(pair_pot.begin[col2]-atom->rho);
#endif
	    atom->rho = pair_pot.begin[col2];
	  }
	} 
#endif
      }	/* first loop over atoms */
#ifdef EAM
      if (eam) {
	for (i=0; i<inconf[h]; i++) {
	  atom = atoms + i + cnfstart[h];
	  typ1 = atom->typ;
	  k    = 3*(cnfstart[h]+i);
	  col  = paircol+ntypes+typ1; /* column of F */
#ifdef PARABEL
	  fnval=parab_comb_ed(&pair_pot,xi,col,atom->rho,&gradF);
#else
	  fnval=splint_comb_ed(&pair_pot,xi,col,atom->rho,&gradF);
#endif
	  forces[config]+= fnval;
	  for (j=0; j<atom->n_neigh; j++) {
	    neigh = atom->neigh+j;
	    if (neigh->nr > i+cnfstart[h]) {
	      typ2  = neigh->typ;
	      col2  = paircol+typ2;
	      r = neigh->r;
	      if ((r < pair_pot.end[col2]) || (r < pair_pot.end[col-ntypes])) {
		grad = (r<pair_pot.end[col2]) ? 
		  splint_grad_ed(&pair_pot,xi,col2,r) : 0.;
#ifdef PARABEL
		gradF2=parab_grad_ed(&pair_pot,xi,col2+ntypes,
				      atoms[neigh->nr].rho);
#else
		gradF2=splint_grad_ed(&pair_pot,xi,col2+ntypes,
				      atoms[neigh->nr].rho);
#endif
		if (typ2 == typ1) 
		  grad2=grad;
		else
		  grad2 = (r < pair_pot.end[col-ntypes]) ? 
		    splint_grad_ed(&pair_pot,xi,col-ntypes,r) : 0.;
		eamforce =  (grad * gradF + grad2 * gradF2) ;
		forces[k  ] += neigh->dist.x * eamforce; 
		forces[k+1] += neigh->dist.y * eamforce; 
		forces[k+2] += neigh->dist.z * eamforce; 
		l    = 3*neigh->nr;
		forces[l  ] -= neigh->dist.x * eamforce; 
		forces[l+1] -= neigh->dist.y * eamforce; 
		forces[l+2] -= neigh->dist.z * eamforce;
	      }
	    }
	  } /* loop over neighbours */
	sum += SQR(forces[k]) + SQR(forces[k+1]) + SQR(forces[k+2]);
	}	/* loop over atoms */
      }	/* if eam is used */
#endif
      forces[config]/=(real) inconf[h]; /* double counting... */
      forces[config]-=force_0[config];
      sum += SQR(forces[config]) + SQR(forces[config+nconf]);
    } /* loop over configurations */
  } /* parallel region */
#ifdef EAM
  if (eam) {
    forces[mdim-2]= DUMMY_WEIGHT * 
      splint_ed(&pair_pot,xi,paircol+DUMMY_COL_RHO,DUMMY_R_RHO)
      - force_0[mdim-2];
    forces[mdim-1]= DUMMY_WEIGHT *
      splint_ed(&pair_pot,xi,DUMMY_COL_PHI,DUMMY_R_PHI)
      - force_0[mdim-1];
    sum+= SQR(forces[mdim-2]);
    sum+= SQR(forces[mdim-1]);
  }
#endif 
  fcalls++;			/* Increase function call counter */
  return sum;
}

