/****************************************************************
* 
* force.c: Routines used for calculating forces/energies in various 
*     interpolation schemes. 
*
*****************************************************************/
/****************************************************************
* $Revision: 1.33 $
* $Date: 2004/11/17 16:00:53 $
*****************************************************************/

#include "potfit.h"

/*****************************************************************************
*
*  compute forces using pair potentials with spline interpolation
*
*  returns sum of squares of differences between calculated and reference
*     values
*
*  arguments: *xi - pointer to potential
*             *forces - pointer to forces calculated from potential
*             flag - used for special tasks
*
* When using the mpi-parallelized version of potfit, all processes but the
* root process jump into this function immediately after initialization and
* stay in here for an infinite loop, to exit only when a certain flag value
* is passed from process 0. When a set of forces needs to be calculated, 
* the root process enters the function with a flag value of 0, broadcasts
* the current potential table xi and the flag value to the other processes,
* thus initiating a force calculation. Whereas the root process returns with
* the result, the other processes stay in the loop. If the root process is
* called with flag value 1, all processes exit the function without 
* calculating the forces.
* If anything changes about the potential beyond the values of the parameters,
* e.g. the location of the sampling points, these changes have to be broadcast
* from rank 0 process to the higher ranked processes. This is done when the
* root process is called with flag value 2. Then a potsync function call is
* initiated by all processes to get the new potential from root.
* 
* xi is the array storing the potential parameters (usually it is the 
*     pair_pot.table - part of the struct pair_pot, but it can also be
*     modified from the current potential.
*
* forces is the array storing the deviations from the reference data, not
*     only for forces, but also for energies, stresses or dummy constraints
*     (if applicable). 
*
* flag is an integer controlling the behaviour of calc_forces_pair. 
*    flag == 1 will cause all processes to exit calc_forces_pair after 
*             calculation of forces.
*    flag == 2 will cause all processes to perform a potsync (i.e. broadcast
*             any changed potential parameters from process 0 to the others)
*             before calculation of forces
*    all other values will cause a set of forces to be calculated. The root 
*             process will return with the sum of squares of the forces,
*             while all other processes remain in the function, waiting for 
*             the next communication initiating another force calculation
*             loop
*
******************************************************************************/

real calc_forces_pair(real *xi, real *forces, int flag)
{
  real tmpsum,sum=0.;
  int first,col1,g;
  real grad0,y0,y1,x0,x1;

  /* This is the start of an infinite loop */
  while (1) {
    tmpsum=0.; 			/* sum of squares of local process */

#ifdef MPI
    /* exchange potential and flag value */
    MPI_Bcast(xi,ndimtot,REAL,0,MPI_COMM_WORLD);
    MPI_Bcast(&flag,1,MPI_INT,0,MPI_COMM_WORLD);
#ifdef EAM
    /* if flag==2 then the potential parameters have changed -> sync */
    if (flag==2) potsync();	
#endif /* EAM */
#endif /* MPI */

    /* init second derivatives for splines */
    for (col1=0; col1<paircol; col1++){  /* just pair potentials */
      first=pair_pot.first[col1];
      x0=pair_pot.begin[col1];
      x1=pair_pot.xcoord[pair_pot.first[col1]+1];
      y0=xi[first];
      y1=xi[first+1];
      /* use power law for inclination at left border */
      if (y0*y1>0)
	grad0=(y0*log(y0/y1)) /  (x0*log(x0/x1));
      else
	grad0=1e30; 		/* natural spline: curvature 0 */
      if (!((grad0>-1e10) && (grad0<1e10))) grad0=1e30;
      if (  (grad0>-1e-20)&& (grad0<1e-20)) grad0=0.;
      if (format==3) 
	spline_ed(pair_pot.step[col1], xi+first,
		  pair_pot.last[col1]-first+1,
		  grad0, 0.0, pair_pot.d2tab+first);
      else 			/* format == 4 ! */
	spline_ne(pair_pot.xcoord+first,xi+first,
		  pair_pot.last[col1]-first+1,
		  grad0, 0.0, pair_pot.d2tab+first);
    }

#ifdef EAM
    for (col1=paircol; col1<paircol+ntypes; col1++) { /* rho */
      first=pair_pot.first[col1];
      if (format==3)
	spline_ed(pair_pot.step[col1], xi+first, 
		  pair_pot.last[col1]-first+1,
		  1e30,0.0,pair_pot.d2tab+first);
      else                   /* format == 4 ! */
	spline_ne(pair_pot.xcoord+first,xi+first,
		  pair_pot.last[col1]-first+1,
		  1e30,0.0,pair_pot.d2tab+first);		  
    }
#ifndef PARABEL 		
/* if we have parabolic interpolation, we don't need that */
    for  (col1=paircol+ntypes; col1<paircol+2*ntypes; col1++) { /* F */
      first=pair_pot.first[col1];
      /* Steigung am linken Rand ist 0 */
      if (format==3)
	spline_ed(pair_pot.step[col1], xi+first, 
		  pair_pot.last[col1]-first+1,
		  0.,1e30,pair_pot.d2tab+first);
      else                   /* format == 4 ! */
	spline_ne(pair_pot.xcoord+first,xi+first,
		  pair_pot.last[col1]-first+1,
		  0.,1e30,pair_pot.d2tab+first);
    }
#endif /* PARABEL */
#endif /* EAM */

#ifndef MPI
    myconf=nconf;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif

    /* region containing loop over configurations, 
       also OMP-parallelized region */
    {

#ifdef EAM
      vektor d;
      int col2;
      real gradF,gradF2,grad2,r,eamforce;
      atom_t *atom2;
#endif
      vektor tmp_force;
      int h,k,i,l,j,typ1,typ2,col,config,stresses;
      real fnval, grad;
      atom_t *atom;
      neigh_t *neigh;

#ifdef _OPENMP
#pragma omp for reduction(+:tmpsum)
#endif
      /* loop over configurations */
      for (h=firstconf; h<firstconf+myconf; h++) {
	config = 3*natoms + h; 	/* slot for energies */
	stresses = 3*natoms + nconf + 6*h; /* slot for stresses */

	/* reset forces */
	forces[config]=0.;
	for (i=stresses; i<stresses+6; i++) 
	  forces[i]=0.;

#ifdef LIMIT
	/* set dummy constraints */
	forces[config+7*nconf]=-force_0[config+7*nconf];
#endif
	/* first loop over atoms: reset forces, densities */
	for (i=0; i<inconf[h]; i++) {
	  k    = 3*(cnfstart[h]+i);
	  forces[k  ] = -force_0[k  ];
	  forces[k+1] = -force_0[k+1];
	  forces[k+2] = -force_0[k+2];
#ifdef EAM
	  conf_atoms[cnfstart[h]-firstatom+i].rho=0.0;
#endif
	} /* end first loop */

        /* 2nd loop: calculate pair forces and energies, atomic densities. */
	for (i=0; i<inconf[h]; i++) {
	  atom = conf_atoms + i + cnfstart[h]-firstatom;
	  typ1 = atom->typ;
	  k    = 3*(cnfstart[h]+i);

	  /* loop over neighbours */
	  for (j=0; j<atom->n_neigh; j++) {
	    neigh = atom->neigh+j;
	    /* only use neigbours with higher numbers, 
	       others are calculated by actio=reactio */
	    if (neigh->nr > i+cnfstart[h]) {
	      typ2  = neigh->typ;
	      /* find correct column */
	      col = (typ1 <= typ2) ? 
		typ1 * ntypes + typ2 - ((typ1 * (typ1 + 1))/2)
		: typ2 * ntypes + typ1 - ((typ2 * (typ2 + 1))/2);
	      if (neigh->r < pair_pot.end[col]) {
		/* fn value and grad are calculated in the same step */
		fnval = splint_comb_dir(&pair_pot,xi,col,
					neigh->slot[0],neigh->shift[0],
					neigh->step[0],&grad);
		forces[config] += fnval; /* not real force: cohesive energy */
		tmp_force.x  = neigh->dist.x * grad;
		tmp_force.y  = neigh->dist.y * grad;
		tmp_force.z  = neigh->dist.z * grad;
		forces[k  ] += tmp_force.x;
		forces[k+1] += tmp_force.y;
		forces[k+2] += tmp_force.z;
		l    = 3*neigh->nr; /* actio = reactio */
		forces[l  ] -= tmp_force.x;
		forces[l+1] -= tmp_force.y;
		forces[l+2] -= tmp_force.z;
#ifdef STRESS
		/* also calculate pair stresses */
		tmp_force.x        *=neigh->r;
		tmp_force.y        *=neigh->r;
		tmp_force.z        *=neigh->r;
		forces[stresses]   -= neigh->dist.x * tmp_force.x;
		forces[stresses+1] -= neigh->dist.y * tmp_force.y;
		forces[stresses+2] -= neigh->dist.z * tmp_force.z;
		forces[stresses+3] -= neigh->dist.x * tmp_force.y;
		forces[stresses+4] -= neigh->dist.y * tmp_force.z;
		forces[stresses+5] -= neigh->dist.z * tmp_force.x;
#endif /* STRESS */
	      }
#ifdef EAM 
	      /* calculate atomic densities */
	      col2 = paircol+typ2;
	      if (typ2==typ1) { /* then transfer(a->b)==transfer(b->a) */
		if (neigh->r < pair_pot.end[col2]) {
		  fnval = splint_dir(&pair_pot,xi,col2,
				     neigh->slot[1],neigh->shift[1],
				     neigh->step[1]);
		  atom->rho += fnval;
		  conf_atoms[neigh->nr - firstatom ].rho += fnval;
		}
	      } else { 		/* transfer(a->b)!=transfer(b->a) */
		col = paircol+typ1;
		if (neigh->r < pair_pot.end[col2]) {
		  atom->rho += splint_dir(&pair_pot,xi,col2,
					  neigh->slot[1],neigh->shift[1],
					  neigh->step[1]);
		}
		/* cannot use slot/shift to access splines */
		if (neigh->r < pair_pot.end[col])
		  conf_atoms[neigh->nr - firstatom].rho += 
		    splint(&pair_pot,xi,col,neigh->r);
	      }	
	      
#endif /* EAM */
	    } /*  neighbours with bigger atom nr */
	  } /* loop over neighbours */
#ifndef EAM /*then we can calculate contribution of forces right away*/

#ifdef FWEIGHT 		   
	  /* Weigh by absolute value of force */
	  forces[k]    /= FORCE_EPS + atom->absforce;
	  forces[k+1]  /= FORCE_EPS + atom->absforce;
	  forces[k+2]  /= FORCE_EPS + atom->absforce;
#endif /* FWEIGHT */

	  /* Returned force is difference between calculated and input force */
	  tmpsum += SQR(forces[k]) + SQR(forces[k+1]) + SQR(forces[k+2]);

#else  /* EAM */

	  col2=paircol+ntypes+typ1; /* column of F */

	  if (atom->rho > pair_pot.end[col2]) {
#ifdef LIMIT   /* then punish target function -> bad potential */
	    forces[config+7*nconf]+=1000*SQR(atom->rho - pair_pot.end[col2]);
#endif /* LIMIT */
#ifndef PARABEL /* then we use the final value, with PARABEL: extrapolate */
	    atom->rho = pair_pot.end[col2];
#endif /* PARABEL */
	  }

	  if (atom->rho < pair_pot.begin[col2]) {
#ifdef LIMIT  /* then punish target function -> bad potential */
	    forces[config+7*nconf]+=1000*SQR(pair_pot.begin[col2]-atom->rho);
#endif /* LIMIT */
#ifndef PARABEL /* then we use the final value, with PARABEL: extrapolate */
	    atom->rho = pair_pot.begin[col2];
#endif /* PARABEL */
	  }

	/* embedding energy, embedding gradient */
	  /* contribution to cohesive energy is difference between 
	     F(n) and F(0) */
#ifdef PARABEL
	  fnval=parab_comb(&pair_pot,xi,col2,atom->rho,atom->gradF)
	    - parab(&pair_pot,xi,col2,0.);
#else
	  fnval=splint_comb(&pair_pot,xi,col2,atom->rho,atom->gradF)
	    - (pair_pot.begin[col2]<=0 ? 
	       splint(&pair_pot,xi,col2,0.) :
	       xi[pair_pot.first[col2]]);
#endif

	  forces[config]+= fnval;

	}




#endif /* EAM */
	}	/* second loop over atoms */

#ifdef EAM /* if we don't have EAM, we're done */
	/* 3rd loop over atom: EAM force */
	for (i=0; i<inconf[h]; i++) {
	  atom = conf_atoms + i + cnfstart[h] - firstatom;
	  typ1 = atom->typ;
	  k    = 3*(cnfstart[h]+i);
	  col  = paircol+ntypes+typ1; /* column of F */


	  for (j=0; j<atom->n_neigh; j++) { /* loop over neighbours */
	    neigh = atom->neigh+j;
            /* only neigbours higher than current atom are of interest */
	    if (neigh->nr > i+cnfstart[h]) { 
	      typ2  = neigh->typ;
	      col2  = paircol+typ2;
	      r = neigh->r;

	      /* are we within reach? */
	      if ((r < pair_pot.end[col2]) || (r < pair_pot.end[col-ntypes])) {
		grad = (r<pair_pot.end[col2]) ? 
		  splint_grad_dir(&pair_pot,xi,col2,neigh->slot[1],
				     neigh->shift[1],neigh->step[1]) 
		  : 0.;

		if (typ2 == typ1) /* use actio = reactio */
		  grad2=grad;
		else
		  grad2 = (r < pair_pot.end[col-ntypes]) ? 
		    splint_grad(&pair_pot,xi,col-ntypes,r) : 0.;

		/* now we know everything - calculate forces */
		eamforce =  (grad * atom->gradF + grad2 * neigh->gradF) ;
		tmp_force.x  = neigh->dist.x * eamforce;
		tmp_force.y  = neigh->dist.y * eamforce;
		tmp_force.z  = neigh->dist.z * eamforce;
		forces[k  ] += tmp_force.x;
		forces[k+1] += tmp_force.y;
		forces[k+2] += tmp_force.z;
		l    = 3*neigh->nr;
		forces[l  ] -= tmp_force.x;
		forces[l+1] -= tmp_force.y;
		forces[l+2] -= tmp_force.z;
#ifdef STRESS
		/* and stresses */
		tmp_force.x        *=neigh->r;
		tmp_force.y        *=neigh->r;
		tmp_force.z        *=neigh->r;
		forces[stresses]   -= neigh->dist.x * tmp_force.x;
		forces[stresses+1] -= neigh->dist.y * tmp_force.y;
		forces[stresses+2] -= neigh->dist.z * tmp_force.z;
		forces[stresses+3] -= neigh->dist.x * tmp_force.y;
		forces[stresses+4] -= neigh->dist.y * tmp_force.z;
		forces[stresses+5] -= neigh->dist.z * tmp_force.x;
#endif /* STRESS */
	      }	/* within reach */
	    } /* higher neigbours */
	  } /* loop over neighbours */

#ifdef FWEIGHT
	  /* Weigh by absolute value of force */
	  forces[k]    /= FORCE_EPS + atom->absforce;
	  forces[k+1]  /= FORCE_EPS + atom->absforce;
	  forces[k+2]  /= FORCE_EPS + atom->absforce;
#endif /* FWEIGHT */

	  /* sum up forces  */
	  tmpsum += SQR(forces[k]) + SQR(forces[k+1]) + SQR(forces[k+2]);

	}	/* third loop over atoms */
#endif /* EAM */

	/* energy contributions */
	forces[config]*=eweight / (real) inconf[h]; 
	forces[config]-=force_0[config];
	tmpsum += SQR(forces[config]);

#ifdef STRESS
	/* stress contributions */
	for (i=stresses;i<stresses+6;i++) {
	  forces[i]*=sweight/conf_vol[h-firstconf];
	  forces[i]-=force_0[i];
	  tmpsum += SQR(forces[i]);

	}
#endif /* STRESS */
#ifdef LIMIT
	/* dummy constraints per configuration */
	tmpsum += SQR(forces[config+7*nconf]);

#endif
      } /* loop over configurations */
    } /* parallel region */

    /* dummy constraints (global) */
#ifdef EAM
    if (myid==0) {
      forces[mdim-(2*ntypes+1)]= DUMMY_WEIGHT * 
	splint(&pair_pot,xi,paircol+DUMMY_COL_RHO,dummy_r)
	- force_0[mdim-(2*ntypes+1)];
#ifdef LIMIT 			/* then we don't need that constraint */
      forces[mdim-(2*ntypes+1)]=0.;
#endif
      tmpsum+= SQR(forces[mdim-(2*ntypes+1)]);

      col1=0;

      for (g=0;g<ntypes;g++) {
	/* this dummy constraint is not used any more... */
	forces[mdim-2*ntypes+g]= 0.;

#ifdef PARABEL
	forces[mdim-ntypes+g]= DUMMY_WEIGHT * /* constraints on U(n) */
	  parab(&pair_pot,xi,paircol+ntypes+g,0.) 
	  - force_0[mdim-ntypes+g];
#else
	if (pair_pot.begin[paircol+ntypes+g]<=0.)
	  /* 0 in domain of U(n) */
	  forces[mdim-ntypes+g]= DUMMY_WEIGHT * /* constraints on U(n) */
	    splint(&pair_pot,xi,paircol+ntypes+g,0.) 
	    - force_0[mdim-ntypes+g];
	else
	  /* 0 not in domain of U(n) */
	  forces[mdim-ntypes+g]=DUMMY_WEIGHT * 
	    xi[pair_pot.first[paircol+ntypes+g]] 
	    - force_0[mdim-ntypes+g];
#endif
	tmpsum+= SQR(forces[mdim-2*ntypes+g]);
	tmpsum+= SQR(forces[mdim-ntypes+g]);

	col1+=ntypes-g;
      }	/* loop over types */
    } /* only root process */
#endif
    sum=tmpsum; 		/* global sum = local sum  */
#ifdef MPI
    /* reduce global sum */
    sum=0.; 
    MPI_Reduce(&tmpsum,&sum,1,REAL,MPI_SUM,0,MPI_COMM_WORLD);
    /* gather forces, energies, stresses */
    MPI_Gatherv(forces+firstatom*3,myatoms,MPI_VEKTOR, /* forces */
		forces,atom_len,atom_dist,MPI_VEKTOR,0,MPI_COMM_WORLD);
    MPI_Gatherv(forces+natoms*3+firstconf,myconf,REAL, /* energies */
		forces+natoms*3,conf_len,conf_dist,REAL,0,MPI_COMM_WORLD);
    /* stresses */
    MPI_Gatherv(forces+natoms*3+nconf+6*firstconf,myconf,MPI_STENS,
		forces+natoms*3+nconf,conf_len,conf_dist,MPI_STENS,
		0,MPI_COMM_WORLD);
#ifdef EAM
    /* punishment constraints */
#ifdef LIMIT
    MPI_Gatherv(forces+natoms*3+7*nconf+firstconf,myconf,REAL, 
		forces+natoms*3+7*nconf,conf_len,conf_dist,REAL,
		0,MPI_COMM_WORLD);
    /* no need to pick up dummy constraints - are already @ root */
#endif /* LIMIT */ 
#endif /* EAM */
#endif /* MPI */

    /* root process exits this function now */
    if (myid==0) {
      fcalls++;			/* Increase function call counter */
      return sum;
    }

    /* non root processes hang on, unless...  */
    if (flag==1) break;		/* Exception: flag 1 means clean up */
  }
  /* once a non-root process arrives here, all is done. */
  return -1.; 
}
