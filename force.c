/****************************************************************
* 
* force.c: Routines used for calculating forces/energies in various 
*     interpolation schemes. 
*
*****************************************************************/
/****************************************************************
* $Revision: 1.28 $
* $Date: 2004/03/17 13:40:29 $
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

real calc_forces_pair(real *xi, real *forces, int flag)
{
  real tmpsum,sum=0.;
  int first,col1,g;
  real grad0,y0,y1,x0,x1;
  /* Define variables globally if not OMP */
#ifndef _OPENMP
#ifdef EAM
  vektor d;
  int col2;
  real gradF,gradF2,grad2,r,eamforce;
  atom_t *atom2;
#endif
  real temp; 			/* dbg */
  vektor tmp_force;
  int h,k,i,l,j,typ1,typ2,col,config,stresses;
  real fnval, grad;
  atom_t *atom;
  neigh_t *neigh;
#endif

  while (1) {
    tmpsum=0.;
#ifdef MPI
    MPI_Bcast(xi,ndimtot,REAL,0,MPI_COMM_WORLD);
    MPI_Bcast(&flag,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
    switch (flag) {
#ifdef EAM
	case 2:			/* time to sync potential... */
#ifdef MPI
	  potsync();
#endif
	  break;
#endif	  
	case 1:
	  break; 		/* We're done - exit... */
	  
	case 0:			/* We actually calc force */
	  
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
		      1e30,1e30,pair_pot.d2tab+first);
	  }
#endif
#endif
#ifndef MPI
	  myconf=nconf;
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
	  {
	    /*  Define OMP thread variables locally */
#ifdef _OPENMP
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
#pragma omp for reduction(+:tmpsum)
#endif
	    for (h=firstconf; h<firstconf+myconf; h++) {
	      config = 3*natoms + h;
	      stresses = 3*natoms + nconf + 6*h; 
	      forces[config]=0.;
	      for (i=stresses; i<stresses+6; i++) 
		forces[i]=0.;
#ifdef LIMIT
	      forces[config+7*nconf]=-force_0[config+7*nconf];
#endif
	      for (i=0; i<inconf[h]; i++) {
		k    = 3*(cnfstart[h]+i);
		forces[k  ] = -force_0[k  ];
		forces[k+1] = -force_0[k+1];
		forces[k+2] = -force_0[k+2];
#ifdef EAM
#ifdef MPI
		conf_atoms[cnfstart[h]-firstatom+i].rho=0.0;
#else MPI
		atoms[cnfstart[h]+i].rho=0.0;
#endif MPI
#endif
	      }
	      for (i=0; i<inconf[h]; i++) {
#ifdef MPI
		atom = conf_atoms + i + cnfstart[h]-firstatom;
#else 
		atom = atoms + i + cnfstart[h];
#endif MPI
		typ1 = atom->typ;
		k    = 3*(cnfstart[h]+i);
#ifdef EAM
#endif
		for (j=0; j<atom->n_neigh; j++) {
		  neigh = atom->neigh+j;
		  if (neigh->nr > i+cnfstart[h]) {
		    typ2  = neigh->typ;
		    col = (typ1 <= typ2) ? 
		      typ1 * ntypes + typ2 - ((typ1 * (typ1 + 1))/2)
		      : typ2 * ntypes + typ1 - ((typ2 * (typ2 + 1))/2);
#ifdef EAM

#endif EAM
		    if (neigh->r < pair_pot.end[col]) {
		      /* not a real force: cohesive energy */
		      /* grad is calculated in the same step */
		      fnval = splint_comb_dir_ed(&pair_pot,xi,col,
						 neigh->slot[0],neigh->shift[0],
						 &grad);
		      forces[config] += fnval;
		      tmp_force.x  = neigh->dist.x * grad;
		      tmp_force.y  = neigh->dist.y * grad;
		      tmp_force.z  = neigh->dist.z * grad;
		      forces[k  ] += tmp_force.x;
		      forces[k+1] += tmp_force.y;
		      forces[k+2] += tmp_force.z;
		      l    = 3*neigh->nr;
		      forces[l  ] -= tmp_force.x;
		      forces[l+1] -= tmp_force.y;
		      forces[l+2] -= tmp_force.z;
#ifdef STRESS
		      tmp_force.x        *=neigh->r;
		      tmp_force.y        *=neigh->r;
		      tmp_force.z        *=neigh->r;
		      forces[stresses]   -= neigh->dist.x * tmp_force.x;
		      forces[stresses+1] -= neigh->dist.y * tmp_force.y;
		      forces[stresses+2] -= neigh->dist.z * tmp_force.z;
		      forces[stresses+3] -= neigh->dist.x * tmp_force.y;
		      forces[stresses+4] -= neigh->dist.y * tmp_force.z;
		      forces[stresses+5] -= neigh->dist.z * tmp_force.x;
#endif STRESS
		    }
#ifdef EAM
		    col2 = paircol+typ2;
		    if (typ2==typ1) {
		      if (neigh->r < pair_pot.end[col2]) {
			fnval = splint_dir_ed(&pair_pot,xi,col2,
					      neigh->slot[1],neigh->shift[1]);
			atom->rho += fnval;
#ifdef MPI
			conf_atoms[neigh->nr - firstatom ].rho += fnval;
#else MPI
			atoms[neigh->nr].rho += fnval;
#endif MPI
		      }
		    } else {
		      col = paircol+typ1;
		      if (neigh->r < pair_pot.end[col2]) {
			atom->rho += splint_dir_ed(&pair_pot,xi,col2,
						   neigh->slot[1],neigh->shift[1]);
		      }
		      if (neigh->r < pair_pot.end[col])
#ifdef MPI
			conf_atoms[neigh->nr - firstatom].rho += 
			  splint_ed(&pair_pot,xi,col,neigh->r);
#else MPI
		      atoms[neigh->nr].rho += 
			splint_ed(&pair_pot,xi,col,neigh->r);
#endif MPI
		    }	
	      
#endif
		  }
		}
#ifndef EAM /*we can't do that here right now*/
		/* Returned force is difference between calculated and input force */
		tmpsum += SQR(forces[k]) + SQR(forces[k+1]) + SQR(forces[k+2]);
#else
		col2=paircol+ntypes+typ1;
		if (atom->rho > pair_pot.end[col2]) {
#ifdef LIMIT
		  forces[config+7*nconf]+=1000*SQR(atom->rho - pair_pot.end[col2]);
#endif
#ifndef PARABEL
		  atom->rho = pair_pot.end[col2];
#endif
		}
/* rho too big... bad thing */
		if (atom->rho < pair_pot.begin[col2]) {
		  /*  forces[config] += 1e10;  rho too small... really bad thing*/
#ifdef LIMIT
		  forces[config+7*nconf]+=1000*SQR(pair_pot.begin[col2]-atom->rho);
#endif
#ifndef PARABEL
		  atom->rho = pair_pot.begin[col2];
#endif
		}
#endif
	      }	/* first loop over atoms */
#ifdef EAM
	      for (i=0; i<inconf[h]; i++) {
#ifdef MPI
		atom = conf_atoms + i + cnfstart[h] - firstatom;
#else MPI
		atom = atoms + i + cnfstart[h];
#endif MPI
		typ1 = atom->typ;
		k    = 3*(cnfstart[h]+i);
		col  = paircol+ntypes+typ1; /* column of F */
#ifdef PARABEL
		fnval=parab_comb_ed(&pair_pot,xi,col,atom->rho,&gradF)
		  - parab_ed(&pair_pot,xi,col,0.);
#else
		fnval=splint_comb_ed(&pair_pot,xi,col,atom->rho,&gradF)
		  - splint_ed(&pair_pot,xi,col,0.);
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
			splint_grad_dir_ed(&pair_pot,xi,col2,neigh->slot[1],
					   neigh->shift[1]) 
			: 0.;

#ifdef PARABEL
#ifdef MPI
		      gradF2=parab_grad_ed(&pair_pot,xi,col2+ntypes,
					   conf_atoms[neigh->nr-firstatom].rho);
#else
		      gradF2=parab_grad_ed(&pair_pot,xi,col2+ntypes,
					   atoms[neigh->nr].rho);
#endif MPI
#else
#ifdef MPI
		      gradF2=splint_grad_ed(&pair_pot,xi,col2+ntypes,
					    conf_atoms[neigh->nr-firstatom].rho);
#else
		      gradF2=splint_grad_ed(&pair_pot,xi,col2+ntypes,
					    atoms[neigh->nr].rho);
#endif MPI
#endif PARABEL
		      if (typ2 == typ1) 
			grad2=grad;
		      else
			grad2 = (r < pair_pot.end[col-ntypes]) ? 
			  splint_grad_ed(&pair_pot,xi,col-ntypes,r) : 0.;
		      eamforce =  (grad * gradF + grad2 * gradF2) ;
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
		      tmp_force.x        *=neigh->r;
		      tmp_force.y        *=neigh->r;
		      tmp_force.z        *=neigh->r;
		      forces[stresses]   -= neigh->dist.x * tmp_force.x;
		      forces[stresses+1] -= neigh->dist.y * tmp_force.y;
		      forces[stresses+2] -= neigh->dist.z * tmp_force.z;
		      forces[stresses+3] -= neigh->dist.x * tmp_force.y;
		      forces[stresses+4] -= neigh->dist.y * tmp_force.z;
		      forces[stresses+5] -= neigh->dist.z * tmp_force.x;
#endif STRESS
		    }
		  }
		} /* loop over neighbours */
		tmpsum += SQR(forces[k]) + SQR(forces[k+1]) + SQR(forces[k+2]);
	      }	/* loop over atoms */
#endif
	      forces[config]*=(real) ENG_WEIGHT / (real) inconf[h]; 
	      forces[config]-=force_0[config];
	      tmpsum += SQR(forces[config]);
#ifdef STRESS
	      for (i=stresses;i<stresses+6;i++) {
#ifdef MPI
		forces[i]*=((real) STRESS_WEIGHT)/conf_vol[h-firstconf];
#else
		forces[i]*=((real) STRESS_WEIGHT)/volumen[h];
#endif
		forces[i]-=force_0[i];
		tmpsum += SQR(forces[i]);
	      }
#endif
#ifdef LIMIT
	      tmpsum += SQR(forces[config+7*nconf]);
#endif
	    } /* loop over configurations */
	  } /* parallel region */
#ifdef EAM
	  if (myid==0) {
	    forces[mdim-(2*ntypes+1)]= DUMMY_WEIGHT * 
	      splint_ed(&pair_pot,xi,paircol+DUMMY_COL_RHO,dummy_r)
	      - force_0[mdim-(2*ntypes+1)];
#ifdef LIMIT
	    forces[mdim-(2*ntypes+1)]=0.;
#endif
	    tmpsum+= SQR(forces[mdim-(2*ntypes+1)]);
	    col1=0;
	    for (g=0;g<ntypes;g++) {
	      forces[mdim-2*ntypes+g]= DUMMY_WEIGHT * /* constraints on phi */
		splint_ed(&pair_pot,xi,col1,dummy_r)
		- force_0[mdim-2*ntypes+g];
#ifdef PARABEL
	      forces[mdim-ntypes+g]= DUMMY_WEIGHT * /* constraints on U(n) */
		parab_ed(&pair_pot,xi,paircol+ntypes+g,0.) 
		- force_0[mdim-ntypes+g];
#else
	      forces[mdim-ntypes+g]= DUMMY_WEIGHT * /* constraints on U(n) */
		splint_ed(&pair_pot,xi,paircol+ntypes+g,0.) 
		- force_0[mdim-ntypes+g];
#endif
	      tmpsum+= SQR(forces[mdim-2*ntypes+g]);
	      tmpsum+= SQR(forces[mdim-ntypes+g]);
	      col1+=ntypes-g;
	    }
	  }
#endif
	  sum=tmpsum;
#ifdef MPI
	  sum=0.;
	  MPI_Reduce(&tmpsum,&sum,1,REAL,MPI_SUM,0,MPI_COMM_WORLD);
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
#endif LIMIT
#endif EAM
#endif MPI
	  if (myid==0) {
//      printf("%1.16g\n",sum);
	    fcalls++;			/* Increase function call counter */
	    return sum;
	  }
	  break;
    }
    if (flag==1) break;
  }
  return -1.; 
}
