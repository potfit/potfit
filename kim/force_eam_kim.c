/****************************************************************
 *
 * force_kim_eam.c: Routine used for calculating eam forces/energies
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
#include "kim.h"
#if defined EAM && !defined COULOMB

#include "../functions.h"
#include "../potential.h"
#include "../splines.h"
#include "../utils.h"


/****************************************************************
 *
 *  compute forces using eam potentials with spline interpolation
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
 * flag is an integer controlling the behaviour of calc_forces_eam.
 *    flag == 1 will cause all processes to exit calc_forces_eam after
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
  int   first, col, i = flag;
  double tmpsum = 0.0, sum = 0.0;
  double *xi = NULL;

  double rho_sum_loc = 0.0, rho_sum = 0.0;

  /* TBEAM: additional s-band contribution */
#ifdef TBEAM
  double rho_s_sum_loc = 0.0, rho_s_sum = 0.0;
#endif /* TBEAM */

  atom_t *atom;
  int   h, j;
  int   n_i, n_j;
  int   self;
  int   uf;
#ifdef APOT
  double temp_eng;
#endif /* APOT */
#ifdef STRESS
  int   us, stresses;
#endif /* STRESS */

  /* pointer for neighbor table */
  neigh_t *neigh;

  /* pair variables */
  double phi_val, phi_grad, r;
  vector tmp_force;

  /* EAM variables */
  int   col_F;
  double eam_force;
  double rho_val, rho_grad, rho_grad_j;
#ifdef TBEAM
  int   col_F_s;
  double rho_s_val, rho_s_grad, rho_s_grad_j;
#endif /* TBEAM */


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
* Declear(define) KIM local variables  
*
***************************************************************************/
	int status;	
	double kim_tmpsum = 0.0;   /* similar to tmp_sum, used by kim */
	double* kimenergy;
	double* kimforce; 			
	double* kimvirial;				 
	int kim_uf = 0;							/* do we need to calculate force? */
	int kim_us = 0;							/* do we need to calculate stress? */
	int numOfconf = nconf;		/* number of configurations in reference data */



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
		status = publish_param(pkimObj[i], &FreeParamAllConfig[i], xi);	
		if (KIM_STATUS_OK > status) {
			KIM_API_report_error(__LINE__, __FILE__, 
														"KIM: publish parameters failed", status);
    	exit(1);
 	 	}
	}

/*
KIM_API_print(pkimObj[0], &status);
*/



/* added ends */


  /* This is the start of an infinite loop */
  while (1) {
    tmpsum = 0.0;		/* sum of squares of local process */
    rho_sum_loc = 0.0;


/* added */
		kim_tmpsum = 0.0;
/* added ends */

   
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


/* added */
/***************************************************************************
* 
* Calculate forces from KIM (general forces, including forces, virial and energy)
*
***************************************************************************/

	kim_uf = uf;
#ifdef STREEE
	kim_us = us;
#endif

	status = calc_force_KIM( pkimObj[h], &kimenergy,  &kimforce,  &kimvirial, kim_uf, kim_us);
	if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__, "KIM: compute forces failed", status);
   	exit(1);
 	}

	/* forces contributation */
  for (i = 0; i < inconf[h]; i++) {
    atom = conf_atoms + i + cnfstart[h] - firstatom;
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
#ifdef FWEIGHT
    /* Weigh by absolute value of force */
    forces[n_i + 0] /= FORCE_EPS + atom->absforce;
    forces[n_i + 1] /= FORCE_EPS + atom->absforce;
    forces[n_i + 2] /= FORCE_EPS + atom->absforce;
#endif /* FWEIGHT */

#ifdef CONTRIB
    if (atom->contrib) 
#endif /* CONTRIB */
		  kim_tmpsum += conf_weight[h] * (dsquare(forces[n_i + 0]) 
			  															+ dsquare(forces[n_i + 1])
				  														+ dsquare(forces[n_i + 2]) );
	}
        
	/* energy contributation */
	forces[energy_p + h] = *kimenergy / (double)inconf[h];
	forces[energy_p + h] -= force_0[energy_p + h];
	kim_tmpsum += conf_weight[h] * eweight * dsquare(forces[energy_p + h]);

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

      }				/* loop over configurations */
    }				/* parallel region */



	tmpsum = kim_tmpsum;

#ifndef KIM

/* added ends */













#ifdef MPI
    /* Reduce rho_sum */
    rho_sum = 0.0;
    MPI_Reduce(&rho_sum_loc, &rho_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef TBEAM
    rho_s_sum = 0.0;
    MPI_Reduce(&rho_s_sum_loc, &rho_s_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif /* TBEAM */
#else /* MPI */
    rho_sum = rho_sum_loc;
#ifdef TBEAM
    rho_s_sum = rho_s_sum_loc;
#endif /* TBEAM */
#endif /* MPI */

    /* dummy constraints (global) */
#ifdef APOT
    /* add punishment for out of bounds (mostly for powell_lsq) */
    if (0 == myid) {
      tmpsum += apot_punish(xi_opt, forces);
    }
#endif /* APOT */











#ifndef NOPUNISH
    if (myid == 0) {
      int   g;
      for (g = 0; g < ntypes; g++) {
#ifndef RESCALE
	/* clear field */
	forces[dummy_p + ntypes + g] = 0.0;	/* Free end... */
	/* NEW: Constraint on U': U'(1.0)=0.0; */
	forces[dummy_p + g] = DUMMY_WEIGHT * splint_grad(&calc_pot, xi, paircol + ntypes + g, 1.0);
#ifdef TBEAM
	/* clear field */
	forces[dummy_p + 3 * ntypes + g] = 0.0;	/* Free end... */
	/* NEW: Constraint on U': U'(1.0)=0.0; */
	forces[dummy_p + 2 * ntypes + g] =
	  DUMMY_WEIGHT * splint_grad(&calc_pot, xi, paircol + 3 * ntypes + g, 1.0);
#endif /* TBEAM */
#else /* !RESCALE */
	forces[dummy_p + ntypes + g] = 0.0;	/* Free end... */
	/* constraints on U`(n) */
	forces[dummy_p + g] =
	  DUMMY_WEIGHT * splint_grad(&calc_pot, xi, paircol + ntypes + g,
	  0.5 * (calc_pot.begin[paircol + ntypes + g] + calc_pot.end[paircol + ntypes + g]))
	  - force_0[dummy_p + g];
#ifdef TBEAM
	forces[dummy_p + 3 * ntypes + g] = 0.0;	/* Free end... */
	/* constraints on U`(n) */
	forces[dummy_p + 2 * ntypes + g] =
	  DUMMY_WEIGHT * splint_grad(&calc_pot, xi, paircol + 3 * ntypes + g,
	  0.5 * (calc_pot.begin[paircol + 3 * ntypes + g] + calc_pot.end[paircol + 3 * ntypes + g]))
	  - force_0[dummy_p + 2 * ntypes + g];
#endif /* TBEAM */
#endif /* !RESCALE */

	/* add punishments to total error sum */
	tmpsum += dsquare(forces[dummy_p + g]);
	tmpsum += dsquare(forces[dummy_p + ntypes + g]);
#ifdef TBEAM
	tmpsum += dsquare(forces[dummy_p + 2 * ntypes + g]);
	tmpsum += dsquare(forces[dummy_p + 3 * ntypes + g]);
#endif /* TBEAM */
      }				/* loop over types */

#ifndef RESCALE
      /* NEW: Constraint on n: <n>=1.0 ONE CONSTRAINT ONLY */
      if (rho_sum > 0.0) {
	/* Calculate averages */
	rho_sum /= (double)natoms;
	/* ATTN: if there are invariant potentials, things might be problematic */
	forces[dummy_p + ntypes] = DUMMY_WEIGHT * (rho_sum - 1.0);
	tmpsum += dsquare(forces[dummy_p + ntypes]);
      }
#ifdef TBEAM
      if (rho_s_sum > 0.0) {
	/* Calculate averages */
	rho_s_sum /= (double)natoms;
	/* ATTN: if there are invariant potentials, things might be problematic */
	forces[dummy_p + 3 * ntypes] = DUMMY_WEIGHT * (rho_s_sum - 1.0);
	tmpsum += dsquare(forces[dummy_p + 3 * ntypes]);
      }
#endif /* TBEAM */
#endif /* !RESCALE */
    }				/* only root process */
#endif /* !NOPUNISH */

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
#ifdef RESCALE
      /* punishment constraints */
      MPI_Gatherv(MPI_IN_PLACE, myconf, MPI_DOUBLE, forces + limit_p,
	conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif /* RESCALE */
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
#ifdef RESCALE
      /* punishment constraints */
      MPI_Gatherv(forces + limit_p + firstconf, myconf, MPI_DOUBLE,
	forces + limit_p, conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif /* RESCALE */
    }
    /* no need to pick up dummy constraints - they are already @ root */
#else
    sum = tmpsum;		/* global sum = local sum  */
#endif /* MPI */


#endif /* !KIM */

/* added */
    sum = tmpsum;
/* ends */



    /* root process exits this function now */
    if (0 == myid) {
      fcalls++;			/* increase function call counter */
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

#endif /* EAM && !COULOMB */


