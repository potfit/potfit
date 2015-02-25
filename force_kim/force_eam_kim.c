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
#include "force_kim.h"
#include "KIM_API_C.h"
#if defined EAM && !defined COULOMB

#include "../functions.h"
#include "../potential.h"
#include "../splines.h"
#include "../utils.h"


/* added */
#define DIM 3
/* added ends */


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
		modelname = "EAM_Dynamo_test_model_for_potfit";
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
		status = PublishParam(pkimObj[i], xi);	
		if (KIM_STATUS_OK > status) {
			KIM_API_report_error(__LINE__, __FILE__, 
														"KIM: publish parameters failed", status);
    	exit(1);
 	 	}
	}


/* make sure that the derivatives does not change */
xi[0] = 0.0;
calc_pot.table[0]=0.0;

/* added ends */


  /* This is the start of an infinite loop */
  while (1) {
    tmpsum = 0.0;		/* sum of squares of local process */
    rho_sum_loc = 0.0;


/* added */
		kim_tmpsum = 0.0;
/* added ends */


#ifdef TBEAM
    rho_s_sum_loc = 0.0;
#endif /* TBEAM */

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

    /* [0, ...,  paircol - 1] = pair potentials */
    /* [paircol, ..., paircol + ntypes - 1] = transfer function */
    for (col = 0; col < paircol + ntypes; col++) {
      first = calc_pot.first[col];
      if (0 == format || 3 == format)
	spline_ed(calc_pot.step[col], xi + first,
	  calc_pot.last[col] - first + 1, *(xi + first - 2), 0.0, calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
	  calc_pot.last[col] - first + 1, *(xi + first - 2), 0.0, calc_pot.d2tab + first);
    }
    /* [paircol + ntypes, ..., paircol + 2 * ntypes - 1] = embedding function */
    for (col = paircol + ntypes; col < paircol + 2 * ntypes; col++) {
      first = calc_pot.first[col];
      /* gradient at left boundary matched to square root function,
         when 0 not in domain(F), else natural spline */
      if (0 == format || 3 == format)
	spline_ed(calc_pot.step[col], xi + first, calc_pot.last[col] - first + 1,
	  *(xi + first - 2), *(xi + first - 1), calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first, calc_pot.last[col] - first + 1,
	  *(xi + first - 2), *(xi + first - 1), calc_pot.d2tab + first);
    }

#ifdef TBEAM
    /* [paircol + 2 * ntypes, ..., paircol + 3 * ntypes - 1] = s-band transfer function */
    for (col = paircol + 2 * ntypes; col < paircol + 3 * ntypes; col++) {
      first = calc_pot.first[col];
      if (0 == format || 3 == format)
	spline_ed(calc_pot.step[col], xi + first,
	  calc_pot.last[col] - first + 1, *(xi + first - 2), 0.0, calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
	  calc_pot.last[col] - first + 1, *(xi + first - 2), 0.0, calc_pot.d2tab + first);
    }

    /* [paircol + 3 * ntypes, ..., paircol + 4 * ntypes - 1] = s-band embedding function */
    for (col = paircol + 3 * ntypes; col < paircol + 4 * ntypes; col++) {
      first = calc_pot.first[col];
      /* gradient at left boundary matched to square root function,
         when 0 not in domain(F), else natural spline */
      if (0 == format || 3 == format)
	spline_ed(calc_pot.step[col], xi + first, calc_pot.last[col] - first + 1,
	  *(xi + first - 2), *(xi + first - 1), calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first, calc_pot.last[col] - first + 1,
	  *(xi + first - 2), *(xi + first - 1), calc_pot.d2tab + first);
    }
#endif /* TBEAM */

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

#ifdef RESCALE
	/* set limiting constraints */
	forces[limit_p + h] = -force_0[limit_p + h];
#endif /* RESCALE */

	/* first loop over atoms: reset forces, densities */
	for (i = 0; i < inconf[h]; i++) {
	  n_i = 3 * (cnfstart[h] + i);
	  if (uf) {
	    forces[n_i + 0] = -force_0[n_i + 0];
	    forces[n_i + 1] = -force_0[n_i + 1];
	    forces[n_i + 2] = -force_0[n_i + 2];
	  } else {
	    forces[n_i + 0] = 0.0;
	    forces[n_i + 1] = 0.0;
	    forces[n_i + 2] = 0.0;
	  }
	  /* reset atomic density */
	  conf_atoms[cnfstart[h] - firstatom + i].rho = 0.0;
#ifdef TBEAM
	  conf_atoms[cnfstart[h] - firstatom + i].rho_s = 0.0;
#endif /* TBEAM */
	}
	/* end of first loop */

	/* 2nd loop: calculate pair forces and energies, atomic densities. */
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
	      }			/* uf */
	    }

	    /* neighbor in range */
	    /* calculate atomic densities */
	    if (atom->type == neigh->type) {
	      /* then transfer(a->b)==transfer(b->a) */
	      if (neigh->r < calc_pot.end[neigh->col[1]]) {
		rho_val = splint_dir(&calc_pot, xi, neigh->slot[1], neigh->shift[1], neigh->step[1]);
		atom->rho += rho_val;

		/* avoid double counting if atom is interacting with a copy of itself */
		if (!self) {
		  conf_atoms[neigh->nr - firstatom].rho += rho_val;
		}
	      }
#ifdef TBEAM
	      if (neigh->r < calc_pot.end[neigh->col[2]]) {
		rho_s_val = splint_dir(&calc_pot, xi, neigh->slot[2], neigh->shift[2], neigh->step[2]);
		atom->rho_s += rho_s_val;
		/* avoid double counting if atom is interacting with a copy of itself */
		if (!self) {
		  conf_atoms[neigh->nr - firstatom].rho_s += rho_s_val;
		}
	      }
#endif /* TBEAM */
	    } else {
	      /* transfer(a->b)!=transfer(b->a) */
	      if (neigh->r < calc_pot.end[neigh->col[1]]) {
		atom->rho += splint_dir(&calc_pot, xi, neigh->slot[1], neigh->shift[1], neigh->step[1]);
	      }
	      /* cannot use slot/shift to access splines */
	      if (neigh->r < calc_pot.end[paircol + atom->type]) {
		conf_atoms[neigh->nr - firstatom].rho +=
		  splint(&calc_pot, xi, paircol + atom->type, neigh->r);
	      }
#ifdef TBEAM
	      if (neigh->r < calc_pot.end[neigh->col[2]]) {
		atom->rho_s += splint_dir(&calc_pot, xi, neigh->slot[2], neigh->shift[2], neigh->step[2]);
	      }
	      /* cannot use slot/shift to access splines */
	      if (neigh->r < calc_pot.end[paircol + 2 * ntypes + atom->type]) {
		conf_atoms[neigh->nr - firstatom].rho_s +=
		  splint(&calc_pot, xi, paircol + 2 * ntypes + atom->type, neigh->r);
	      }
#endif /* TBEAM */
	    }
	  }			/* loop over all neighbors */

	  /* column of F */
	  col_F = paircol + ntypes + atom->type;
#ifdef TBEAM
	  /* column of F of the s-band */
	  col_F_s = col_F + 2 * ntypes;
#endif /* TBEAM */

#ifdef RESCALE
	  /* we punish the potential for bad behavior:
	   * if the density of one atom is smaller or greater than we have the
	   * embedding function tabulated a punishment is added */

	  if (atom->rho > calc_pot.end[col_F]) {
	    /* then punish target function -> bad potential */
	    forces[limit_p + h] += DUMMY_WEIGHT * 10.0 * dsquare(atom->rho - calc_pot.end[col_F]);
	    atom->rho = calc_pot.end[col_F];
	  }

	  if (atom->rho < calc_pot.begin[col_F]) {
	    /* then punish target function -> bad potential */
	    forces[limit_p + h] += DUMMY_WEIGHT * 10.0 * dsquare(calc_pot.begin[col_F] - atom->rho);
	    atom->rho = calc_pot.begin[col_F];
	  }
#ifdef TBEAM
	  if (atom->rho_s > calc_pot.end[col_F_s]) {
	    /* then punish target function -> bad potential */
	    forces[limit_p + h] += DUMMY_WEIGHT * 10.0 * dsquare(atom->rho_s - calc_pot.end[col_F_s]);
	    atom->rho_s = calc_pot.end[col_F_s];
	  }

	  if (atom->rho_s < calc_pot.begin[col_F_s]) {
	    /* then punish target function -> bad potential */
	    forces[limit_p + h] += DUMMY_WEIGHT * 10.0 * dsquare(calc_pot.begin[col_F_s] - atom->rho_s);
	    atom->rho_s = calc_pot.begin[col_F_s];
	  }
#endif /* TBEAM */
#endif /* RESCALE */

	  /* embedding energy, embedding gradient */
	  /* contribution to cohesive energy is F(n) */

#ifndef RESCALE
	  if (atom->rho < calc_pot.begin[col_F]) {
#ifdef APOT
	    /* calculate analytic value explicitly */
	    apot_table.fvalue[col_F] (atom->rho, xi_opt + opt_pot.first[col_F], &temp_eng);
	    atom->gradF = apot_grad(atom->rho, xi_opt + opt_pot.first[col_F], apot_table.fvalue[col_F]);
	    forces[energy_p + h] += temp_eng;
#else
	    /* linear extrapolation left */
	    rho_val = splint_comb(&calc_pot, xi, col_F, calc_pot.begin[col_F], &atom->gradF);
	    forces[energy_p + h] += rho_val + (atom->rho - calc_pot.begin[col_F]) * atom->gradF;
#endif /* APOT */
	  } else if (atom->rho > calc_pot.end[col_F]) {
#ifdef APOT
	    /* calculate analytic value explicitly */
	    apot_table.fvalue[col_F] (atom->rho, xi_opt + opt_pot.first[col_F], &temp_eng);
	    atom->gradF = apot_grad(atom->rho, xi_opt + opt_pot.first[col_F], apot_table.fvalue[col_F]);
	    forces[energy_p + h] += temp_eng;
#else
	    /* and right */
	    rho_val =
	      splint_comb(&calc_pot, xi, col_F, calc_pot.end[col_F] - 0.5 * calc_pot.step[col_F],
	      &atom->gradF);
	    forces[energy_p + h] += rho_val + (atom->rho - calc_pot.end[col_F]) * atom->gradF;
#endif /* APOT */
	  } else {		/* and in-between */
#ifdef APOT
	    /* calculate small values directly */
	    if (atom->rho < 0.1) {
	      apot_table.fvalue[col_F] (atom->rho, xi_opt + opt_pot.first[col_F], &temp_eng);
	      atom->gradF = apot_grad(atom->rho, xi_opt + opt_pot.first[col_F], apot_table.fvalue[col_F]);
	      forces[energy_p + h] += temp_eng;
	    } else
#endif
	      forces[energy_p + h] += splint_comb(&calc_pot, xi, col_F, atom->rho, &atom->gradF);
	  }
#else
	  forces[energy_p + h] += splint_comb(&calc_pot, xi, col_F, atom->rho, &atom->gradF);
#endif /* !RESCALE */

	  /* sum up rho */
	  rho_sum_loc += atom->rho;

#ifdef TBEAM
#ifndef RESCALE
	  if (atom->rho_s < calc_pot.begin[col_F_s]) {
#ifdef APOT
	    /* calculate analytic value explicitly */
	    apot_table.fvalue[col_F_s] (atom->rho_s, xi_opt + opt_pot.first[col_F_s], &temp_eng);
	    atom->gradF_s =
	      apot_grad(atom->rho_s, xi_opt + opt_pot.first[col_F_s], apot_table.fvalue[col_F_s]);
	    forces[energy_p + h] += temp_eng;
#else
	    /* linear extrapolation left */
	    rho_s_val = splint_comb(&calc_pot, xi, col_F_s, calc_pot.begin[col_F_s], &atom->gradF_s);
	    forces[energy_p + h] += rho_s_val + (atom->rho_s - calc_pot.begin[col_F_s]) * atom->gradF_s;
#endif /* APOT */
	  } else if (atom->rho_s > calc_pot.end[col_F_s]) {
#ifdef APOT
	    /* calculate analytic value explicitly */
	    apot_table.fvalue[col_F_s] (atom->rho_s, xi_opt + opt_pot.first[col_F_s], &temp_eng);
	    atom->gradF_s =
	      apot_grad(atom->rho_s, xi_opt + opt_pot.first[col_F_s], apot_table.fvalue[col_F_s]);
	    forces[energy_p + h] += temp_eng;
#else
	    /* and right */
	    rho_s_val =
	      splint_comb(&calc_pot, xi, col_F_s, calc_pot.end[col_F_s] - 0.5 * calc_pot.step[col_F_s],
	      &atom->gradF_s);
	    forces[energy_p + h] += rho_s_val + (atom->rho_s - calc_pot.end[col_F_s]) * atom->gradF_s;
#endif /* APOT */
	  }
	  /* and in-between */
	  else {
#ifdef APOT
	    /* calculate small values directly */
	    if (atom->rho_s < 0.1) {
	      apot_table.fvalue[col_F_s] (atom->rho_s, xi_opt + opt_pot.first[col_F_s], &temp_eng);
	      atom->gradF_s =
		apot_grad(atom->rho_s, xi_opt + opt_pot.first[col_F_s], apot_table.fvalue[col_F_s]);
	      forces[energy_p + h] += temp_eng;
	    } else
#endif
	      forces[energy_p + h] += splint_comb(&calc_pot, xi, col_F_s, atom->rho_s, &atom->gradF_s);
	  }
#else
	  forces[energy_p + h] += splint_comb(&calc_pot, xi, col_F_s, atom->rho_s, &atom->gradF_s);
#endif /* !RESCALE */

	  /* sum up rho_s */
	  rho_s_sum_loc += atom->rho_s;
#endif /* TBEAM */
	}			/* second loop over atoms */

	/* 3rd loop over atom: EAM force */
	if (uf) {		/* only required if we calc forces */
	  for (i = 0; i < inconf[h]; i++) {
	    atom = conf_atoms + i + cnfstart[h] - firstatom;
	    n_i = 3 * (cnfstart[h] + i);
	    for (j = 0; j < atom->num_neigh; j++) {
	      /* loop over neighbors */
	      neigh = atom->neigh + j;
	      /* In small cells, an atom might interact with itself */
	      self = (neigh->nr == i + cnfstart[h]) ? 1 : 0;
	      col_F = paircol + ntypes + atom->type;	/* column of F */
#ifdef TBEAM
	      col_F_s = col_F + 2 * ntypes;
#endif /* TBEAM */
	      r = neigh->r;
	      /* are we within reach? */
	      if ((r < calc_pot.end[neigh->col[1]]) || (r < calc_pot.end[col_F - ntypes])) {
		rho_grad =
		  (r < calc_pot.end[neigh->col[1]]) ? splint_grad_dir(&calc_pot, xi, neigh->slot[1],
		  neigh->shift[1], neigh->step[1]) : 0.0;
		if (atom->type == neigh->type)	/* use actio = reactio */
		  rho_grad_j = rho_grad;
		else
		  rho_grad_j =
		    (r < calc_pot.end[col_F - ntypes]) ? splint_grad(&calc_pot, xi, col_F - ntypes, r) : 0.0;
		/* now we know everything - calculate forces */
		eam_force = (rho_grad * atom->gradF + rho_grad_j * conf_atoms[(neigh->nr) - firstatom].gradF);

#ifdef TBEAM			/* s-band contribution to force for TBEAM */
		if ((r < calc_pot.end[neigh->col[2]]) || (r < calc_pot.end[col_F_s - ntypes])) {
		  rho_s_grad =
		    (r < calc_pot.end[neigh->col[2]]) ? splint_grad_dir(&calc_pot, xi, neigh->slot[2],
		    neigh->shift[2], neigh->step[2]) : 0.0;
		  if (atom->type == neigh->type) {	/* use actio = reactio */
		    rho_s_grad_j = rho_s_grad;
		  } else {
		    rho_s_grad_j = (r < calc_pot.end[col_F_s - ntypes]) ?
		      splint_grad(&calc_pot, xi, col_F_s - ntypes, r) : 0.0;
		  }
		  /* now we know everything - calculate forces */
		  eam_force +=
		    (rho_s_grad * atom->gradF_s + rho_s_grad_j * conf_atoms[(neigh->nr) - firstatom].gradF_s);
		}
#endif /* TBEAM */

		/* avoid double counting if atom is interacting with a copy of itself */
		if (self)
		  eam_force *= 0.5;
		tmp_force.x = neigh->dist_r.x * eam_force;
		tmp_force.y = neigh->dist_r.y * eam_force;
		tmp_force.z = neigh->dist_r.z * eam_force;
		forces[n_i + 0] += tmp_force.x;
		forces[n_i + 1] += tmp_force.y;
		forces[n_i + 2] += tmp_force.z;
		/* actio = reactio */
		n_j = 3 * neigh->nr;
		forces[n_j + 0] -= tmp_force.x;
		forces[n_j + 1] -= tmp_force.y;
		forces[n_j + 2] -= tmp_force.z;

#ifdef STRESS
		/* and stresses */
		if (us) {
		  forces[stresses + 0] -= neigh->dist.x * tmp_force.x;
		  forces[stresses + 1] -= neigh->dist.y * tmp_force.y;
		  forces[stresses + 2] -= neigh->dist.z * tmp_force.z;
		  forces[stresses + 3] -= neigh->dist.x * tmp_force.y;
		  forces[stresses + 4] -= neigh->dist.y * tmp_force.z;
		  forces[stresses + 5] -= neigh->dist.z * tmp_force.x;
		}
#endif /* STRESS */
	      }			/* within reach */
	    }			/* loop over neighbours */

#ifdef FWEIGHT
	    /* Weigh by absolute value of force */
	    forces[n_i + 0] /= FORCE_EPS + atom->absforce;
	    forces[n_i + 1] /= FORCE_EPS + atom->absforce;
	    forces[n_i + 2] /= FORCE_EPS + atom->absforce;
#endif /* FWEIGHT */

	    /* sum up forces  */
#ifdef CONTRIB
	    if (atom->contrib)
#endif /* CONTRIB */
	      tmpsum += conf_weight[h] *
		(dsquare(forces[n_i + 0]) + dsquare(forces[n_i + 1]) + dsquare(forces[n_i + 2]));
	  }			/* third loop over atoms */
	}

	/* use forces */
	/* energy contributions */
	forces[energy_p + h] /= (double)inconf[h];
	forces[energy_p + h] -= force_0[energy_p + h];
	tmpsum += conf_weight[h] * eweight * dsquare(forces[energy_p + h]);

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
#ifdef RESCALE
	/* limiting constraints per configuration */
	tmpsum += conf_weight[h] * dsquare(forces[limit_p + h]);
#endif /* RESCALE */



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
* Publish parameters 
*
* KIM global variabls:
* embeddingData
* densityData
* rPhiData
*
* potfit global variables: 
* ntypes
* opt_pot.first
* opt_pot.last
* opt_pot.ncols   # number of total potentials 
* rcutmax
* opt_pot.step		# step size of potentials 
* potfit global variables: 
* ntypes
* paircol         # number of pair potentials(\Phi) in EAM
* opt_pot.table   # potential table
* opt_pot.first
* opt_pot.last
* opt_pot.ncols   # number of total potentials 

***************************************************************************/

int PublishParam(void* pkim, double* PotTable)
{	
	int status;
	/* published param */
	double* param_cutoff;
	double* param_deltaRho;
	double* param_deltaR;
  double* param_embedding;
	double* param_density;
	double* param_rPhi;
	/* local variables */
	int potIdx;				/* which potfit potential we are visiting? ranging from 0 to ncols.
										 	 The data are arranged as \Phi(r),..,\rho(r),...,U(\rho) in 
											 opt_pot.table. For example, if there are 2 species i j, the data 
											 will be arranged as \Phi(ii), \Phi(ij), \Phi(jj), \rho(i), \rho(j)
											 U(i),U(j)*/
	int tableIdx;			/* where we are in the potfit potential table? */			
	int i,j,k;
	int shape[3];



/* Access to free parameters */
	KIM_API_getm_data(pkim, &status, 3*3,
										"PARAM_FREE_cutoff",	   		 &param_cutoff,	 		 1,
										"PARAM_FREE_deltaRho",	 		 &param_deltaRho,		 1,
										"PARAM_FREE_deltaR",	 			 &param_deltaR,			 1);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);
    return status;
  }

	/* Update parameters */
	/* The data are arranged as \Phi(r),..,\rho(r),...,U(\rho) in opt_pot.table
	so the first ``0'' is realatd to r data, and the last ``ncols - 1'' is 
	related to rho data */
	*param_cutoff    = rcutmax;
	*param_deltaRho  = opt_pot.step[opt_pot.ncols - 1];
	*param_deltaR    = opt_pot.step[0];


	/* Updata potential data */	
	/* shape[0] : number of species
		 shape[1] : number of U(\rho) point for embedding data
		 shape[1] : number of species for density and rPhi data
		 shape[2] : number of potential data points for density and rPhi data
	*/
  /* copy embedding data */	
	param_embedding = (double*) GetShapeAndPointer(pkim,"PARAM_FREE_embeddingData", 
																								 shape, &status);

	for (i = 0; i<shape[0]; i++) {	
		potIdx = shape[0]*(shape[0] + 1)/2 + shape[0] + i;		
		for (k = 0; k < shape[1]; k++) {
			tableIdx = opt_pot.first[potIdx] + k;
			param_embedding[i*shape[1] + k] = PotTable[tableIdx];
		}
		
	}

	/* copy density data */
	param_density   = (double*) GetShapeAndPointer(pkim,"PARAM_FREE_densityData", 
																								 shape, &status);
	for (i = 0; i<shape[0]; i++) {		
		potIdx = shape[0]*(shape[0] + 1)/2 + i;
		for (k = 0; k < shape[2]; k++) {
			tableIdx = opt_pot.first[potIdx] + k;
			param_density[i*shape[1] + 0*shape[2] + k] = PotTable[tableIdx];
		}	
		/* fill in remaining columns of density data */
		for (j = 1; j < shape[1]; j++) {
			for (k = 0; k < shape[2]; k++) {
				param_density[i*shape[1] + j*shape[2] + k] = param_density[i*shape[1] + 0*shape[2] + k];
		 	}
		}
	}


	/* copy Phi data to upper-triangular part of rPhiData*/
	param_rPhi = (double*) GetShapeAndPointer(pkim,"PARAM_FREE_rPhiData", 
																								shape, &status);

	for (i = 0; i<shape[0]; i++) {	
		for ( j = i; j < shape[1]; j++) {
			potIdx = (i*shape[0] - i*(i - 1)/2) + (j-i);				
			for (k = 0; k < shape[2]; k++) {
				tableIdx = opt_pot.first[potIdx] + k;
				param_rPhi[i*shape[1] + j*shape[2] + k] = PotTable[tableIdx];	
				/* this following is used in replace of the above line to publish 
				r*Phi data, not just Phi data. /*
				/*param_rPhi[i*shape[1] + j*shape[2] + k] = PotTable[tableIdx]*k*(*param_deltaR);				
				*/			
			}
		}
		/* filling in lower-triangular part of rPhiData */
    for (j = 0; j < i; j++)
    {
      for (k = 0; k < shape[2]; k++)
      {
        param_rPhi[i*shape[1] + j*shape[2] + k] =	param_rPhi[j*shape[1] + i*shape[2] + k] ;
      }
    }
	}


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
* Get Shape and Pointer
*
* Shape stores the number of species and number of potentials data poitnter
* Pointer points to the data
*
***************************************************************************/

void* GetShapeAndPointer(void* const pkim, char const* const argName, 
												 int* const shape, int* const status)
{
   void* retValue;
   KIM_API_get_shape((void*) pkim, (char*) argName, (int*) shape, status);
   if (KIM_STATUS_OK > *status)
   {
      KIM_API_report_error(__LINE__, __FILE__,"KIM_API_get_shape", *status);
      KIM_API_model_destroy(pkim);
      KIM_API_free((void**)&pkim, status);
      return 0;
   }
   retValue = KIM_API_get_data((void*) pkim, (char*) argName, status);
   if (KIM_STATUS_OK > *status)
   {
      KIM_API_report_error(__LINE__, __FILE__,"KIM_API_get_data", *status);
      KIM_API_model_destroy(pkim);
      KIM_API_free((void**)&pkim, status);
      return 0;
   }

   return retValue;
}



