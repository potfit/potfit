/****************************************************************
 *
 * force_pair.c: Routines used for calculating pair forces/energies
 *
 ****************************************************************
 *
 * Copyright 2002-2014
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.sourceforge.net/
 *
 ****************************************************************
 *
 *   This file is part of potfit.
 *
 *   potfit is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   potfit is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

#include "../potfit.h"
#include "force_kim.h"

#ifdef PAIR

#include "../functions.h"
#include "../potential.h"
#include "../splines.h"
#include "../utils.h"


/* added */
#define DIM 3
/* added ends */


/****************************************************************
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
 * xi_opt is the array storing the potential parameters (usually it is the
 *     opt_pot.table - part of the struct opt_pot, but it can also be
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
 ****************************************************************/

double calc_forces(double *xi_opt, double *forces, int flag)
{
  int   first, col, i;
  double *xi = NULL;

  /* Some useful temp variables */
  double tmpsum = 0.0, sum = 0.0;

  atom_t *atom;
  int   h, j;
  int   n_i, n_j;
  int   self;
  int   uf;
#ifdef STRESS
  int   us, stresses;
#endif /* STRESS */

  /* pointer for neighbor table */
  neigh_t *neigh;

  /* pair variables */
  double phi_val, phi_grad;
  vector tmp_force;


/* added */
/***************************************************************************
*
* Declear(define) KIM local variables  
*
***************************************************************************/
	int status;	
	double kim_tmpsum = 0.0;   /* similar to tmp_sum, used by kim */
	double* kimenergy;
	double* kimforce; 			
	double* kimvirial;				 
	int numOfconf = nconf;		/* number of configurations in reference data */
/* added ends*/


  switch (format) {
      case 0:
	xi = calc_pot.table;
	break;
      case 3:			/* fall through */
      case 4:
	xi = xi_opt;		/* calc-table is opt-table */
	break;
      case 5:
	xi = calc_pot.table;	/* we need to update the calc-table */
  }


/* added */
/***************************************************************************
*
* Create KIM object
*
*	create them once for all 
***************************************************************************/
	if (! haveKIMObj) {
	
		char *modelname;
		modelname = "Pair_Lennard_Jones_Truncated_Nguyen_Ar__MO_398194508715_000";
		/*NOTE(change): need to read kim model name from potfit param file */
		
	  /* Allocate memory for KIM objects */
		pkimObj = (void**) malloc(numOfconf * sizeof(void *));
		if (NULL == pkimObj) {
		  KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
			exit(1);
		}
		
		/* Checking whether .kim files in test and model are compatible or not */
		for (i = 0; i < numOfconf; i++) {  				
			status = KIM_API_file_init(&pkimObj[i], "descriptor.kim", modelname);
  		if (KIM_STATUS_OK > status)
  		{
    		KIM_API_report_error(__LINE__, __FILE__, "KIM_API_file_init", status);
    		exit(1);
 	  	}
		}
		
		/* QUESTION: does each config has the same number of species? e.g. there are
		two types of species, but a specific config may one have one species). If
		not, then ntypes below need to be modified from config to config */

		/*  Initialize KIM objects */
		for (i = 0; i < numOfconf; i++) { 
			status = CreateKIMObj(pkimObj[i], inconf[i], ntypes, cnfstart[i]);
			if (KIM_STATUS_OK > status)	{
    		KIM_API_report_error(__LINE__, __FILE__,
															"KIM: initializing objects failed", status);
    		exit(1);
 	  	}
		}

		/* only need to initialize KIM objects once */
		haveKIMObj = 1;

		printf("Initializing KIM objects ... done\n");
		fflush(stdout);		
	}
/* added ends */


 /* unittest */
/* This is used to see that the neighbor list is correct. Note that c and cpp is
zerolist but fortran is onelist. So if fortran models are used(as is here), j 
starts from one. */ 
/*
	int request;
  int currentPart;
  int numOfPartNeigh;
  double* Rij_list;
  int* neighListOfCurrentPart;

for( j = 1; j <= inconf[1]; j++) {
 KIM_API_get_neigh(pkimObj[1], 1, j, &currentPart, &numOfPartNeigh,
                              &neighListOfCurrentPart, &Rij_list);
printf("which atom: %d\n",j);

i = 0; 
	printf("%d %f %f %f\n", neighListOfCurrentPart[i], 
													Rij_list[DIM*i],
													Rij_list[DIM*i+1],
											  	Rij_list[DIM*i+2]	);

printf("which atom: %d\n",j);
													
i = numOfPartNeigh-1; 
	printf("%d %f %f %f\n", neighListOfCurrentPart[i], 
													Rij_list[DIM*i],
													Rij_list[DIM*i+1],
											  	Rij_list[DIM*i+2]	);
}
exit(1);
*/
 /* unittest ends*/
 
 

/* added */
/***************************************************************************
*
* Publish KIM parameters
* 
* publish parameters for LJ potential: epsilon and sigma.
*
* Need to use xi_opt not xi (the second  and third argument) when calling
* PublishParam. xi_opt stores the original parameter for analytic
* potential, but xi here stores some parameter suitable for spline calculating.
* potfit uses splines to do force calculating, even though it is analytic
* potential. But KIM uses analytic potential directly. 
* 
* NOTE(change): the cutoff is also published here, since for the time being,
* we want touse the cutoff from potfit to see whether it can give the same 
*	resuts. In the future, it will come from KIM models directly, and do not 
*	need to be published. 
***************************************************************************/
 	for (i = 0; i < numOfconf; i++) {
		status = PublishParam(pkimObj[i], xi_opt, rcut[0]);	
		if (KIM_STATUS_OK > status) {
			KIM_API_report_error(__LINE__, __FILE__, 
														"KIM: publish parameters failed", status);
    	exit(1);
 	 	}
	}
/* added ends */


  /* This is the start of an infinite loop */
  while (1) {
    tmpsum = 0.0;		/* sum of squares of local process */


/* added */
		kim_tmpsum = 0.0;
/* added ends */


#if defined APOT && !defined MPI
    if (0 == format) {
      apot_check_params(xi_opt);
      update_calc_table(xi_opt, xi, 0);
    }
#endif /* APOT && !MPI */

#ifdef MPI
#ifndef APOT
    /* exchange potential and flag value */
    MPI_Bcast(xi, calc_pot.len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif /* APOT */
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (1 == flag)
      break;			/* Exception: flag 1 means clean up */

#ifdef APOT
    if (0 == myid)
      apot_check_params(xi_opt);
    MPI_Bcast(xi_opt, ndimtot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    update_calc_table(xi_opt, xi, 0);
#else /* APOT */
    /* if flag==2 then the potential parameters have changed -> sync */
    if (2 == flag)
      potsync();
#endif /* APOT */
#endif /* MPI */

    /* init second derivatives for splines */

    /* pair potentials */
    for (col = 0; col < paircol; col++) {
      first = calc_pot.first[col];
      if (0 == format || 3 == format)
	spline_ed(calc_pot.step[col], xi + first,
	  calc_pot.last[col] - first + 1, *(xi + first - 2), 0.0, calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
	  calc_pot.last[col] - first + 1, *(xi + first - 2), 0.0, calc_pot.d2tab + first);
    }

#ifndef MPI
    myconf = nconf;
#endif /* MPI */

    /* region containing loop over configurations */
    {

      /* loop over configurations */
      for (h = firstconf; h < firstconf + myconf; h++) {
	uf = conf_uf[h - firstconf];
#ifdef STRESS
	us = conf_us[h - firstconf];
#endif /* STRESS */
	/* reset energies and stresses */
	forces[energy_p + h] = 0.0;
#ifdef STRESS
	stresses = stress_p + 6 * h;
	for (i = 0; i < 6; i++)
	  forces[stresses + i] = 0.0;
#endif /* STRESS */

#ifdef APOT
	if (enable_cp)
	  forces[energy_p + h] += chemical_potential(ntypes, na_type[h], xi_opt + cp_start);
#endif /* APOT */

	/* first loop over atoms: reset forces, densities */
	for (i = 0; i < inconf[h]; i++) {
	  if (uf) {
	    n_i = 3 * (cnfstart[h] + i);
	    forces[n_i + 0] = -force_0[n_i + 0];
	    forces[n_i + 1] = -force_0[n_i + 1];
	    forces[n_i + 2] = -force_0[n_i + 2];
	  } else {
	    n_i = 3 * (cnfstart[h] + i);
	    forces[n_i + 0] = 0.0;
	    forces[n_i + 1] = 0.0;
	    forces[n_i + 2] = 0.0;
	  }
	}
	/* end first loop */

	/* 2nd loop: calculate pair forces and energies */
	for (i = 0; i < inconf[h]; i++) {
	  atom = conf_atoms + i + cnfstart[h] - firstatom;
	  n_i = 3 * (cnfstart[h] + i);
	  /* loop over neighbors */
	  for (j = 0; j < atom->num_neigh; j++) {		
	    neigh = atom->neigh + j;

	    /* In small cells, an atom might interact with itself */
	    self = (neigh->nr == i + cnfstart[h]) ? 1 : 0;

	    /* pair potential part */
	    if (neigh->r < calc_pot.end[neigh->col[0]]) {
	      /* fn value and grad are calculated in the same step */
	      if (uf)
		phi_val =
		  splint_comb_dir(&calc_pot, xi, neigh->slot[0], neigh->shift[0], neigh->step[0], &phi_grad);
	      else
		phi_val = splint_dir(&calc_pot, xi, neigh->slot[0], neigh->shift[0], neigh->step[0]);



/* unittest */
/* Replace the spline-based force calculation algorithom in potfit by an analytic 
	Lennard-Jones potential to check that the difference in optimized parameters
	between KIM and potfit is due to the method how the forces are calcualted. In 
	KIM, analytic function is used, but potfit uses splines */
/*
 AnalyticForce(xi_opt[2], xi_opt[3], calc_pot.end[neigh->col[0]],
 								neigh->r, &phi_val, &phi_grad);
*/
/* unittest ends*/



	      /* avoid double counting if atom is interacting with a copy of itself */
	      if (self) {
		phi_val *= 0.5;
		phi_grad *= 0.5;
	      }

	      /* add cohesive energy */
	      forces[energy_p + h] += phi_val;

	      /* calculate forces */
	      if (uf) {
		tmp_force.x = neigh->dist_r.x * phi_grad;
		tmp_force.y = neigh->dist_r.y * phi_grad;
		tmp_force.z = neigh->dist_r.z * phi_grad;
		forces[n_i + 0] += tmp_force.x;
		forces[n_i + 1] += tmp_force.y;
		forces[n_i + 2] += tmp_force.z;
		/* actio = reactio */
		n_j = 3 * neigh->nr;
		forces[n_j + 0] -= tmp_force.x;
		forces[n_j + 1] -= tmp_force.y;
		forces[n_j + 2] -= tmp_force.z;

	

/* unittest (to use this, print the same thing in the model (model driver), 
   then compare. Note that one needs to `make' `make install' and 
	 `make isntall-set-dafault-to-vx' for kim. vx is the main version of kim )
*/
/* print atom number, neighbor number, neigh distance, cutoff
printf("%3d %3d %15.5e %15.5e ", n_i, j, neigh->r, calc_pot.end[neigh->col[0]] );
*/
/* print accumulated forces, note that the forces include reference forces 
printf("%15.5e\n", forces[n_i + 0]);
*/

/* print components used to calculate forces */
/*
printf("potfit forces each step %3d %f %f %f %f\n",j, phi_grad, neigh->dist_r.x, neigh->dist_r.y, neigh->dist_r.z  );
*/
/* unittest ends */



#ifdef STRESS
		/* also calculate pair stresses */
		if (us) {
		  forces[stresses + 0] -= neigh->dist.x * tmp_force.x;
		  forces[stresses + 1] -= neigh->dist.y * tmp_force.y;
		  forces[stresses + 2] -= neigh->dist.z * tmp_force.z;
		  forces[stresses + 3] -= neigh->dist.x * tmp_force.y;
		  forces[stresses + 4] -= neigh->dist.y * tmp_force.z;
		  forces[stresses + 5] -= neigh->dist.z * tmp_force.x;
		}
#endif /* STRESS */
	      }
	    }			/* neighbor in range */
	 }			/* loop over all neighbors */

	  /* then we can calculate contribution of forces right away */
	  if (uf) {
#ifdef FWEIGHT
	    /* Weigh by absolute value of force */
	    forces[n_i + 0] /= FORCE_EPS + atom->absforce;
	    forces[n_i + 1] /= FORCE_EPS + atom->absforce;
	    forces[n_i + 2] /= FORCE_EPS + atom->absforce;
#endif /* FWEIGHT */

	    /* sum up forces */
#ifdef CONTRIB
	    if (atom->contrib)
#endif /* CONTRIB */
	      tmpsum += conf_weight[h] *
		(dsquare(forces[n_i + 0]) + dsquare(forces[n_i + 1]) + dsquare(forces[n_i + 2]));
	  }
	}			/* second loop over atoms */



/* unittest */
/* used together with `printf("flag kimtmp_sum of forces from kim%f\n", kim_tmpsum);' 
 (located a few lines below, to see that the forces calculated via KIM are the same
 as the potfit built-in ones.  */
/*
printf("flag tmp_sum of forces from potfit %f\n",tmpsum);
*/
/* uinttest ends */


	/* energy contributions */
	forces[energy_p + h] /= (double)inconf[h];
	forces[energy_p + h] -= force_0[energy_p + h];
	tmpsum += conf_weight[h] * eweight * dsquare(forces[energy_p + h]);


/* unittest */
/* used together with `printf("flag kimtmp_sum of energy from kim%f\n", kim_tmpsum);' 
 (located a few lines below, to see that the forces calculated via KIM are the same
 as the potfit built-in ones.  */
/*
printf("flag tmp_sum of energy from potfit %f\n",tmpsum);
*/
/* uinttest ends */


#ifdef STRESS
	/* stress contributions */
	if (uf && us) {
	  for (i = 0; i < 6; i++) {
	    forces[stresses + i] /= conf_vol[h - firstconf];
	    forces[stresses + i] -= force_0[stresses + i];
	    tmpsum += conf_weight[h] * sweight * dsquare(forces[stresses + i]);
	  }
	}
#endif /* STRESS */

/* added */
/***************************************************************************
* 
* Calculate forces from KIM (general forces, including forces, virial and energy)
*
***************************************************************************/
/* the following 20 lines is the same as CalcForce, either one could be used*/
	/*  unpack data from KIM */
/*	KIM_API_getm_data(pkimObj[h], &status, 2*3,
                    "energy",              &kimenergy,              1,
                    "forces",              &kimforce,               1   ,
										"virial",           	 &virial,       			 1 );
  if (KIM_STATUS_OK > status)
  {
		KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);
    return status;
  }	
*/

  /* Call model compute */
/*  status = KIM_API_model_compute(pkimObj[h]);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_compute", status);
    return status;
  }
*/


	status = CalcForce( pkimObj[h], &kimenergy,  &kimforce,  &kimvirial);
	if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__, "KIM: compute forces failed", status);
   	exit(1);
 	}

	/* forces contributation */
	for (i = 0; i < inconf[h]; i++) {
	  n_i = DIM * (cnfstart[h] + i);	
	  if (uf) {
	    forces[n_i + 0] =  kimforce[DIM*i + 0] - force_0[n_i + 0];
	    forces[n_i + 1] =  kimforce[DIM*i + 1] - force_0[n_i + 1];
	    forces[n_i + 2] =  kimforce[DIM*i + 2] - force_0[n_i + 2];
	  } else {
	    forces[n_i + 0] = 0.0;
	    forces[n_i + 1] = 0.0;
	    forces[n_i + 2] = 0.0;
	  }
		kim_tmpsum += conf_weight[h] * (  dsquare(forces[n_i + 0]) 
																		+ dsquare(forces[n_i + 1])
																		+ dsquare(forces[n_i + 2]) );
	}


/* unittest */
/* used together with `printf("flag tmp_sum of forces from potfit %f\n",tmpsum);' 
 (located a few lines below, to see that the forces calculated via KIM are the same
 as the potfit built-in ones.  */
/*
printf("flag kimtmp_sum of forces from kim%f\n", kim_tmpsum);
*/
/* uinttest ends */


	/* energy contributation */
#ifdef APOT
	if (enable_cp)
	  forces[energy_p + h] = chemical_potential(ntypes, na_type[h], xi_opt + cp_start);
#endif /* APOT */

	forces[energy_p + h] = *kimenergy / (double)inconf[h];
	forces[energy_p + h] -= force_0[energy_p + h];
	kim_tmpsum += conf_weight[h] * eweight * dsquare(forces[energy_p + h]);


/* unittest */
/* used together with `printf("flag tmp_sum of energy from potfit %f\n",tmpsum);' 
 (located a few lines below, to see that the forces calculated via KIM are the same
 as the potfit built-in ones.  */
/*
printf("flag kimtmp_sum of energy from kim%f\n", kim_tmpsum);
*/
/* uinttest ends */


#ifdef STRESS
	/* stress contributions */
	if (uf && us) {
  for (i = 0; i < 6; i++) {
		forces[stresses + i] 	= kimvirial[i]; 
		forces[stresses + i] /= conf_vol[h - firstconf];
	  forces[stresses + i] -= force_0[stresses + i];
    kim_tmpsum += conf_weight[h] * sweight * dsquare(forces[stresses + i]);
	  }
	}	
#endif /* STRESS */
/* added ends */

      }				/* loop over configurations */
    }				/* parallel region */



/* added */
/***************************************************************************
*
* Replace sum of errors calculated from potfit by that from KIM 
*
***************************************************************************/

	tmpsum = kim_tmpsum;

/* added ends */


    /* dummy constraints (global) */
#ifdef APOT
    /* add punishment for out of bounds (mostly for powell_lsq) */
    if (0 == myid) {
      tmpsum += apot_punish(xi_opt, forces);
    }
#endif /* APOT */

#ifdef MPI
    /* reduce global sum */
    sum = 0.0;
    MPI_Reduce(&tmpsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    /* gather forces, energies, stresses */
    if (0 == myid) {		/* root node already has data in place */
      /* forces */
      MPI_Gatherv(MPI_IN_PLACE, myatoms, MPI_VECTOR, forces,
	atom_len, atom_dist, MPI_VECTOR, 0, MPI_COMM_WORLD);
      /* energies */
      MPI_Gatherv(MPI_IN_PLACE, myconf, MPI_DOUBLE, forces + energy_p,
	conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#ifdef STRESS
      /* stresses */
      MPI_Gatherv(MPI_IN_PLACE, myconf, MPI_STENS, forces + stress_p,
	conf_len, conf_dist, MPI_STENS, 0, MPI_COMM_WORLD);
#endif /* STRESS */
    } else {
      /* forces */
      MPI_Gatherv(forces + firstatom * 3, myatoms, MPI_VECTOR,
	forces, atom_len, atom_dist, MPI_VECTOR, 0, MPI_COMM_WORLD);
      /* energies */
      MPI_Gatherv(forces + energy_p + firstconf, myconf, MPI_DOUBLE,
	forces + energy_p, conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#ifdef STRESS
      /* stresses */
      MPI_Gatherv(forces + stress_p + 6 * firstconf, myconf, MPI_STENS,
	forces + stress_p, conf_len, conf_dist, MPI_STENS, 0, MPI_COMM_WORLD);
#endif /* STRESS */
    }
#else
    sum = tmpsum;		/* global sum = local sum  */
#endif /* MPI */

    /* root process exits this function now */
    if (0 == myid) {
      fcalls++;			/* Increase function call counter */
      if (isnan(sum)) {
#ifdef DEBUG
	printf("\n--> Force is nan! <--\n\n");
#endif /* DEBUG */
	return 10e10;
      } else
	return sum;
    }
  }				/* end of infinite loop */

  /* once a non-root process arrives here, all is done. */
  return -1.0;
}

#endif /* PAIR */



/* added */	
/***************************************************************************
* 
* Create KIM object	 
*
* transferred varialbes:
* pkim: KIM object pointer
* Natoms: number of atoms in this configuration  
* Nspecies: number of atom species in this configuraiton 
*	start: index of the first atom of this configuration in potfit atom array
*
* golbal potfit varialbes, used but not transferred: 
* atoms 
* elements
* box_side_len  (this is defined by us, not potfit) 
***************************************************************************/

int CreateKIMObj(void* pkim, int Natoms, int Nspecies, int start)
{
  /* model inputs */
  int* numberOfParticles;
  int* numberOfSpecies;
  int* particleSpecies;
  double* coords;	
	int* numberContrib;
	NeighObjectType* NeighObject;
	const char* NBCstr;
	int NBC;
	/* local variables */
  int status;	
	int species_code;	
	int halfflag;
	int	Ncontrib;	/* number of contriburting particles */	
	int i, j, k;
	/* We have to allocate additional memory for neighbor list, since potfit 
		assigns the neighbors of each particle to atoms[i].neigh[k].nr (nr is
		the index neighbor k of atom i.). The memory is not continuous.
		The following variable (neighListLength), together with `BeginIdx' in
		the NeighObject, are used to gather the neighbor info to continuous memory. 
	*/	
	int neighListLength; /* total length of neighList */


	/* set value */
	Ncontrib = Natoms;

	/* allocate memory for NeighObject ( simulator should be in charge of
	allocating memory for neighbor objects, not kim model ) */
	NeighObject = (NeighObjectType*) malloc(sizeof(NeighObjectType));
  if (NULL == NeighObject) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}	

	/* Allocate memory via the KIM system */
  KIM_API_allocate(pkim, Natoms, Nspecies, &status);
  if (KIM_STATUS_OK > status)
  {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_allocate", status);
    return(status);
  }

	/* register for neighObject */
	KIM_API_setm_data(pkim, &status, 1*4,
				   					 "neighObject",     1,   NeighObject,   1);
  if (KIM_STATUS_OK > status) {
   KIM_API_report_error(__LINE__, __FILE__,"KIM_API_setm_data",status);
    return(status);		
	 }

	/* register for get_neigh */
  status = KIM_API_set_method(pkim, "get_neigh", 1, (func_ptr) &get_neigh);
  if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__,"KIM_API_set_method",status);
    return(status);				
	}	
	/* call Model's init routine */
  status = KIM_API_model_init(pkim);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_init", status);
    return(status);		
  }

/* Determine which neighbor list type to use */
  halfflag = KIM_API_is_half_neighbors(pkim, &status);
  if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__,"is_half_neighbors", status);
    return(status);		
	}

  /* Unpack data from KIM object */
  KIM_API_getm_data(pkim, &status, 5*3,
                    "numberOfParticles",   &numberOfParticles,   1,
                    "numberOfSpecies",     &numberOfSpecies,     1,
                    "particleSpecies",     &particleSpecies,     1,
                    "coordinates",         &coords,              1,
					 "numberContributingParticles",  &numberContrib, (1==halfflag) );
  if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);
  	return status;
	}

  /* Set values */
  *numberOfParticles = Natoms;
  *numberOfSpecies 	 = Nspecies;	
	*numberContrib		 = Ncontrib;

	/* set boxSideLengths if MI_OPBC is used */
	/* determine which NBC is used */
  status = KIM_API_get_NBC_method(pkim, &NBCstr);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method", status);
    return status; }
  if ((!strcmp("NEIGH_RVEC_H",NBCstr)) || (!strcmp("NEIGH_RVEC_F",NBCstr))) {
    NBC = 0; }
  else if ((!strcmp("NEIGH_PURE_H",NBCstr)) || (!strcmp("NEIGH_PURE_F",NBCstr))) {
    NBC = 1; }
  else if ((!strcmp("MI_OPBC_H",NBCstr)) || (!strcmp("MI_OPBC_F",NBCstr))) {
    NBC = 2; }
  else if (!strcmp("CLUSTER",NBCstr)) {
    NBC = 3; }
  else {
    status = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "Unknown NBC method", status);
    return status; }

	if (NBC == 2) {
		/* define local varialbe */
		double* boxSideLen;
		int which_conf;		/* which config we are in? */
	
		which_conf = atoms[start].conf; 

  	/* Unpack data from KIM object */
  	KIM_API_getm_data(pkim, &status, 1*3,
					 					"boxSideLengths",      &boxSideLen,			1 );
  	if (KIM_STATUS_OK > status)
 	  {
			KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);
 		 	return status;
		}
 
	  /* Set values */
		boxSideLen[0]	= box_side_len[DIM*which_conf + 0];
		boxSideLen[1]	= box_side_len[DIM*which_conf + 1];
		boxSideLen[2] = box_side_len[DIM*which_conf + 2];
	}		


  /* set the species types */
	/* to use this, the #C Header line in configuration file has to be included.
		Also the order of elements name after #C should be in accordane	with the
		species code in the first column of configuration data, ranging from 0 to
		(ntypes-1).  
		Here, assume potfit ensures this.
		e.g. if we have `#C Ar Ne', then the code for Ar and Ne should be 0 and 1,
		respectively. */
		
	/* QUESTION: in KIM_API_get_species_code, doest the second argument has to
		 to be totally the same as that kim descriptor file? e.g. Upper case or lower
		 case matters or not? */

	for (i = 0; i < *numberOfParticles; i++) { 
		/* in potfit, atom types range from 0 to (ntype-1) */
		if (atoms[start+i].type < 0 || atoms[start+i].type >= *numberOfSpecies) {
 			status = KIM_STATUS_FAIL;			
			KIM_API_report_error(__LINE__, __FILE__,
                           "Unexpected species detected", status);
  		return status;
		}
		else {
			j = atoms[start+i].type;
			species_code =  KIM_API_get_species_code(pkim, elements[j], &status);	
		  if (KIM_STATUS_OK > status)
  		{
   			KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_species_code," 
														"the species names need to be exactly the same"
														"as that in KIM standard file", status);
 			 	return status;
 			}
 	    particleSpecies[i] = species_code;
		}
	}

	/* set coords values */
	for (i = 0; i < *numberOfParticles; i++) {
		coords[DIM*i]   = atoms[start+i].pos.x;
		coords[DIM*i+1] = atoms[start+i].pos.y;
		coords[DIM*i+2] = atoms[start+i].pos.z;
	}

	/* calcualte the length of neighborList */
	neighListLength = 0;
	for (i = 0; i < *numberOfParticles; i++) {
		neighListLength += atoms[start+i].num_neigh;	
	}

  /* allocate memory for NeighObject members */
  NeighObject->NNeighbors = (int*) malloc((*numberOfParticles) * sizeof(int));
  if (NULL==NeighObject->NNeighbors) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}	
  NeighObject->neighborList = (int*) malloc(neighListLength * sizeof(int));
  if (NULL==NeighObject->neighborList) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);		
	}	
 NeighObject->RijList = (double*) malloc((DIM*neighListLength) * sizeof(double));
  if (NULL==NeighObject->RijList) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);		
	}	
	NeighObject->BeginIdx	= (int*) malloc((*numberOfParticles) * sizeof(int));
	if (NULL==NeighObject->BeginIdx) {
		KIM_API_report_error(__LINE__,__FILE__,"malloc unsuccessful", -1);
		exit(1);		
	}

	/* copy the number of neighbors to NeighObject.NNeighbors */
	status = KIM_STATUS_FAIL; /* assume staus fails */
	for (i = 0; i < *numberOfParticles; i++) {
		NeighObject->NNeighbors[i] = atoms[start+i].num_neigh;	
	}
	if (i == *numberOfParticles){
		status = KIM_STATUS_OK;
	}
  if (KIM_STATUS_OK > status)	{
   	KIM_API_report_error(__LINE__, __FILE__,
												"copy number of neighbors to NNeighbors failed", status);
  return status;
 	}

	/* copy neighborlist from distributed memory locations to 
		continuous ones.
		copy atoms[i].neigh[j].nr to NeighObject.neighborList,
		copy atoms[i].neigh[k].dist.x ... to 	NeighObject.Rij[DIM*k] ... */	
	status = KIM_STATUS_FAIL; /* assume staus fails */	
	k = 0;
	for (i = 0; i < *numberOfParticles; i++) {
		NeighObject->BeginIdx[i] = k;
		for (j = 0; j < NeighObject->NNeighbors[i]; j++) {
			NeighObject->neighborList[k]  =	atoms[start+i].neigh[j].nr - start;
			NeighObject->RijList[DIM*k]   =	atoms[start+i].neigh[j].dist.x;
			NeighObject->RijList[DIM*k+1] =	atoms[start+i].neigh[j].dist.y;
			NeighObject->RijList[DIM*k+2] =	atoms[start+i].neigh[j].dist.z;
			k++;

/* unittest */
/* used together with to see that the stuff in the neighbor lista are correct
	set  start == 0 , will check for the first config. We set j==0, and j==last
	atom in the neighbor list of an atom, since we don't want that verbose info.
*/
/*
if (start != 0 && (j== 0 || j == atoms[start+i].num_neigh-1 )) {
printf("last neighbor: %d\n", atoms[start+ *numberOfParticles-1].num_neigh);
printf("which atom: %d\n",i);
printf("%d %f %f %f\n", atoms[start+i].neigh[j].nr,
												atoms[start+i].neigh[j].dist.x,
												atoms[start+i].neigh[j].dist.y,
												atoms[start+i].neigh[j].dist.z );
}
*/
/* unittest ends*/
		}
	}
	
	if (i == *numberOfParticles && k == neighListLength){
		status = KIM_STATUS_OK;
	}
  if (KIM_STATUS_OK > status)	{
   	KIM_API_report_error(__LINE__, __FILE__,
												"copy neighbor list failed", status);
  	return status; } 

	/* If the number of neighbors of an atom is zero, set the BeginIdx to the 
		 last position in neghbor list. Actually, the main purpose is to ensure 
		 that the BeginIdx of the last atom will not go beyond limit of 
		 neighborlist length.
		 e.g. there are 128 atoms in the config, and we use half neighbor list,
		 then the 128th atom will have no neighbors. Then the the begin index
		 for the last atom, BeginIdx[127] will go beyond the limit of Neighbor
		 list, which may result in segfault. So we need to do something to avoid
		 this.
		 I'm sure, there are better ways to do this.
	*/	
	for (i = 0; i < *numberOfParticles; i++) {
		if( NeighObject->NNeighbors[i] == 0) {
			NeighObject->BeginIdx[i] = k-1;
		}
	}

 	return KIM_STATUS_OK;
}
/* added ends */


/* added */
/***************************************************************************
* 
* get_neigh 
*
***************************************************************************/

int get_neigh(void* kimmdl, int *mode, int *request, int* part,
                       int* numnei, int** nei1part, double** Rij)
{
   /* local variables */
  intptr_t* pkim = *((intptr_t**) kimmdl);
  int partToReturn;
  int status;
  int* numberOfParticles;
	int idx; /* index of first neighbor of each particle*/
  NeighObjectType* nl;

  /* initialize numnei */
  *numnei = 0;

	/* unpack neighbor list object */
  numberOfParticles = (int*) KIM_API_get_data(pkim, "numberOfParticles", &status);
  if (KIM_STATUS_OK > status) {
  KIM_API_report_error(__LINE__, __FILE__,"get_data", status);
	}

  nl = (NeighObjectType*) KIM_API_get_data(pkim, "neighObject", &status);
  if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__,"get_data", status);
	}

  /* check mode and request */
  if (0 == *mode) /* iterator mode */
  {
		if (0 == *request) /* reset iterator */
    {
	    (*nl).iteratorId = -1;
      return KIM_STATUS_NEIGH_ITER_INIT_OK;
		}
    else if (1 == *request) /* increment iterator */
    {
  	  (*nl).iteratorId++;
      if ((*nl).iteratorId >= *numberOfParticles)
      {
				return KIM_STATUS_NEIGH_ITER_PAST_END;
      }
      else
      {
				partToReturn = (*nl).iteratorId;
      }
    }
    else /* invalid request value */
    {
 	    KIM_API_report_error(__LINE__, __FILE__,"Invalid request in get_periodic_neigh",
													 KIM_STATUS_NEIGH_INVALID_REQUEST);
      return KIM_STATUS_NEIGH_INVALID_REQUEST;
		}
 	}
  else if (1 == *mode) /* locator mode */
  {
    if ((*request >= *numberOfParticles) || (*request < 0)) /* invalid id */
    {
			KIM_API_report_error(__LINE__, __FILE__,"Invalid part ID in get_periodic_neigh",
	 												KIM_STATUS_PARTICLE_INVALID_ID);
			return KIM_STATUS_PARTICLE_INVALID_ID;
    }
    else
    {
      partToReturn = *request;
    }
  }
  else /* invalid mode */
  {
		KIM_API_report_error(__LINE__, __FILE__,"Invalid mode in get_periodic_neigh",
		 											KIM_STATUS_NEIGH_INVALID_MODE);
    return KIM_STATUS_NEIGH_INVALID_MODE;
  }

	/* index of the first neigh of each particle */
	idx = (*nl).BeginIdx[partToReturn];

  /* set the returned part */
  *part = partToReturn;

  /* set the returned number of neighbors for the returned part */
  *numnei = (*nl).NNeighbors[partToReturn];

  /* set the location for the returned neighbor list */
  *nei1part = &((*nl).neighborList[idx]);

  /* set the pointer to Rij to appropriate value */
  *Rij = &((*nl).RijList[DIM*idx]);

	return KIM_STATUS_OK;
}

/* added ends */



/* added */
/***************************************************************************
*
* Publish KIM parameter ((update epsilon and sigma of LJ))
*
* transferred varialbes:
*	pkkm: KIM ojbect
*	PotTable: potential table where the potential paramenters are stored
*	CutOff: cutoff radius   
*
***************************************************************************/
/*Don't forget to publish params for all KIM objest */

int PublishParam(void* pkim, double* PotTable, double CutOff)
{	
  int status;
	/* published param */
	double* param_cutoff;
	double* param_epsilon;
	double* param_sigma;
	
	/* Access to free parameters, also unpack other variables as needed */
	KIM_API_getm_data(pkim, &status, 3*3,
										"PARAM_FREE_cutoff",	 &param_cutoff,				 1,
										"PARAM_FREE_epsilon",	 &param_epsilon,			 1,
										"PARAM_FREE_sigma",		 &param_sigma,				 1 );
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);
    return status;
  }

	/* Update parameters */
	*param_cutoff  = CutOff; 
	/*first two entry of parameter table stores gradients, so we need to 
	read parameters starting from index 2 */
	*param_epsilon = PotTable[2]; 
	*param_sigma	 = PotTable[3]; 
 	status = KIM_API_model_reinit(pkim);
  if (KIM_STATUS_OK > status)	
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_reinit", status);
    return status;
  }
	return KIM_STATUS_OK;
}
/* added ends */



/* added */
/***************************************************************************
* 
* Calculate force from KIM (general force, including force, virial, energy)
*
* transferred variaves 
*	Pkim: KIM object
* the following three are also model output 
*	energy;
* force;
*	virial;
*
***************************************************************************/

int CalcForce(void* pkim, double** energy, double** force, double** virial)
{	
	/* local variables */
  int status;

	/* Access to free parameters, also unpack other variables as needed */
	KIM_API_getm_data(pkim, &status, 2*3,
                    "energy",              energy,              1,
                    "forces",              force,               1 /*  ,
										"virial",           	 virial,       			 1 */);
  if (KIM_STATUS_OK > status)
  {
		KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);
    return status;
  }

  /* Call model compute */
  status = KIM_API_model_compute(pkim);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_compute", status);
    return status;
  }

	return KIM_STATUS_OK;
}
/* added ends */



/* added */
/***************************************************************************
* 
* Calculate analytic Lennard-Jones potential value and gradient
*
* This is just a test, we don't need it once it is checked
*
*	transferred variable:
* phi_val: potential value
* phi_grad: gradient of potential 
***************************************************************************/

int AnalyticForce(double epsilon, double sigma, double cutoff,
									double r, double* phi_val, double* phi_grad)
{	
	/* local variables */
	double sor;
	double sor6;
	double sor12;

	sor   = sigma/(double)r;
	sor6  = pow(sor,6);
	sor12 = pow(sor6,2);

	if( r > cutoff) {
		*phi_val  = 0.0;
		*phi_grad = 0.0;
	} else {
		*phi_val  = 4.0*epsilon*(sor12 - sor6);
		*phi_grad = 24.0*epsilon*(-2.0*sor12 + sor6)/(double)r;
	}
	return 1;
}

