/****************************************************************
*
* force.c: Routines used for calculating pair/monopole/dipole
*     forces/energies in various interpolation schemes.
*
*****************************************************************/
/*
*   Copyright 2002-2010 Peter Brommer, Franz G"ahler, Daniel Schopf, 
                        Philipp Beck
*             Institute for Theoretical and Applied Physics
*             University of Stuttgart, D-70550 Stuttgart, Germany
*             http://www.itap.physik.uni-stuttgart.de/
*
*****************************************************************/
/*
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
*   along with potfit; if not, write to the Free Software
*   Foundation, Inc., 51 Franklin St, Fifth Floor,
*   Boston, MA  02110-1301  USA
*/
/****************************************************************
* $Revision: 1.2 $
* $Date: 2010/04/22 13:29:40 $
*****************************************************************/

#ifdef DIPOLE
#include "potfit.h"

/*****************************************************************************
*
*  compute forces using pair potentials with spline interpolation
*
*  returns sum of squares of differences between calculated and reference
*     values
*
*  arguments: *xi and *xi_d - pointer to potentials
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
******************************************************************************/

real calc_forces_dipole(real *xi_opt, real *forces, int flag)
{
  real  tmpsum, sum = 0.;
  int   first, col, ne, size, i;
  real *xi = NULL;
  real *xi_d = NULL;
  apot_table_t *apt = &apot_table;

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

  xi_d = calc_pot.table_dipole;
  ne =  apot_table.total_ne_par;
  size = apt->number;

  /* This is the start of an infinite loop */
  while (1) {
    tmpsum = 0.;		/* sum of squares of local process */
#ifndef APOT
    if (format > 4 && myid == 0)
      update_calc_table(xi_opt, xi, 0);
#endif /* APOT */

#if defined APOT && !defined MPI
    if (format == 0) {
      apot_check_params(xi_opt);
      update_calc_table(xi_opt, xi, 0);
    }
#endif

#ifdef MPI
    /* exchange potential and flag value */
#ifndef APOT
    MPI_Bcast(xi, calc_pot.len, REAL, 0, MPI_COMM_WORLD);
#endif /* APOT */
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef APOT
    if (myid == 0)
      apot_check_params(xi_opt);
    MPI_Bcast(xi_opt, ndimtot, REAL, 0, MPI_COMM_WORLD);
    if (format == 0)
      update_calc_table(xi_opt, xi, 0);
    if (flag == 1)
      break;
#else /* APOT */

    /* if flag==2 then the potential parameters have changed -> sync */
    if (flag == 2)
      potsync();
    /* non root processes hang on, unless...  */
    if (flag == 1)
      break;			/* Exception: flag 1 means clean up */
#endif /* APOT */
#endif /* MPI */

    /* init second derivatives for splines */
    for (col = 0; col < paircol; col++) {
      first = calc_pot.first[col];
      if (format == 3 || format == 0)
	spline_ed(calc_pot.step[col], xi + first,
		  calc_pot.last[col] - first + 1,
		  *(xi + first - 2), 0.0, calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
		  calc_pot.last[col] - first + 1,
		  *(xi + first - 2), 0.0, calc_pot.d2tab + first);
    }
    spline_ed(calc_pot.step[0], xi_d + calc_pot.first[0],
	      calc_pot.last[0] - calc_pot.first[0] + 1,
	      *(xi_d + calc_pot.first[0] - 2), 0.0, calc_pot.d2tab_dipole + calc_pot.first[0]);
    

#ifndef MPI
    myconf = nconf;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif

    /* region containing loop over configurations,
       also OMP-parallelized region */
    {
      int   self;
      vector tmp_force;
      int   h, j, k, l, typ1, typ2, uf, us, stresses;	// config
      real  fnval, grad, fnval_tail, grad_tail, eval_i, eval_j, p_sr_tail;
      atom_t *atom;
      neigh_t *neigh;

#ifdef _OPENMP
#pragma omp for reduction(+:tmpsum,rho_sum_loc)
#endif
      
      /* reset dipoles and fields: LOOP Z E R O */
      for (i = 0; i < natoms; i++) {
	atoms[i].E_stat.x = 0;
	atoms[i].E_stat.y = 0;
	atoms[i].E_stat.z = 0;
	atoms[i].p_sr.x = 0;
	atoms[i].p_sr.y = 0;
	atoms[i].p_sr.z = 0;
	atoms[i].E_ind.x = 0;
	atoms[i].E_ind.y = 0;
	atoms[i].E_ind.z = 0;
	atoms[i].p_ind.x = 0;
	atoms[i].p_ind.y = 0;
	atoms[i].p_ind.z = 0;
	atoms[i].E_old.x = 0;
	atoms[i].E_old.y = 0;
	atoms[i].E_old.z = 0;
	atoms[i].E_temp.x = 0;
	atoms[i].E_temp.y = 0;
	atoms[i].E_temp.z = 0;
      }

      /* loop over configurations: M A I N LOOP CONTAINING ALL ATOM-LOOPS */
      for (h = firstconf; h < firstconf + myconf; h++) {
	uf = conf_uf[h - firstconf];
	us = conf_us[h - firstconf];
	/* reset energies and stresses */
	forces[energy_p + h] = 0.;
	for (i = 0; i < 6; i++)
	  forces[stress_p + 6 * h + i] = 0.;

	/* F I R S T LOOP OVER ATOMS: reset forces, dipoles */
	for (i = 0; i < inconf[h]; i++) {    //atoms
	  if (uf) {
	    k = 3 * (cnfstart[h] + i);
	    forces[k] = -force_0[k];
	    forces[k + 1] = -force_0[k + 1];
	    forces[k + 2] = -force_0[k + 2];
	  } else {
	    k = 3 * (cnfstart[h] + i);
	    forces[k] = 0.;
	    forces[k + 1] = 0.;
	    forces[k + 2] = 0.;
	  }
	}
	/* end F I R S T LOOP */

	/* S E C O N D loop: calculate short-range and monopole forces,
	   calculate static field- and dipole-contributions */
	for (i = 0; i < inconf[h]; i++) {       //atoms
	  atom = conf_atoms + i + cnfstart[h] - firstatom;
	  typ1 = atom->typ;
	  k = 3 * (cnfstart[h] + i);
	  for (j = 0; j < atom->n_neigh; j++) {    //neighbors
	    neigh = atom->neigh + j;
	    typ2 = neigh->typ;
	    /* In small cells, an atom might interact with itself */
	    self = (neigh->nr == i + cnfstart[h]) ? 1 : 0;

	      /* calculate short-range forces */
	      if (neigh->r < calc_pot.end[neigh->col[0]]) {
		/* fn value and grad are calculated in the same step */
		if (uf) {
		  fnval =
		    splint_comb_dir(&calc_pot, xi, neigh->slot[0],
				    neigh->shift[0], neigh->step[0], &grad);
		} else {
		  fnval =
		    splint_dir(&calc_pot, xi, neigh->slot[0], neigh->shift[0],
			       neigh->step[0]);
		}

		/* avoid double counting if atom is interacting with a
		   copy of itself */
		if (self) {
		  fnval *= 0.5;
		  grad *= 0.5;
		}
		forces[energy_p + h] += fnval;

		if (uf) {
		  tmp_force.x = neigh->dist.x * grad;
		  tmp_force.y = neigh->dist.y * grad;
		  tmp_force.z = neigh->dist.z * grad;
		  forces[k] += tmp_force.x;
		  forces[k + 1] += tmp_force.y;
		  forces[k + 2] += tmp_force.z;
		  /* actio = reactio */
		  l = 3 * neigh->nr;
		  forces[l] -= tmp_force.x;
		  forces[l + 1] -= tmp_force.y;
		  forces[l + 2] -= tmp_force.z;

#ifdef STRESS
		  /* calculate pair stresses */
		  if (us) {
		    tmp_force.x *= neigh->r;
		    tmp_force.y *= neigh->r;
		    tmp_force.z *= neigh->r;
		    stresses = stress_p + 6 * h;
		    forces[stresses] -= neigh->dist.x * tmp_force.x;
		    forces[stresses + 1] -= neigh->dist.y * tmp_force.y;
		    forces[stresses + 2] -= neigh->dist.z * tmp_force.z;
		    forces[stresses + 3] -= neigh->dist.x * tmp_force.y;
		    forces[stresses + 4] -= neigh->dist.y * tmp_force.z;
		    forces[stresses + 5] -= neigh->dist.z * tmp_force.x;
		  }
#endif /* STRESS */
		}
	      }

	      /* calculate monopole forces */
	      if (neigh->r < dp_cut) { 

		if (uf) {
		  fnval_tail = splint_comb_dir_dipole(&calc_pot, xi_d, 
						      neigh->slot[0] - APOT_STEPS * neigh->col[0],
						      neigh->shift[0],
						      neigh->step[0], &grad_tail);
		} else {
		  fnval_tail = splint_dir(&calc_pot, xi_d, 
					  neigh->slot[0] - APOT_STEPS * neigh->col[0],
					  neigh->shift[0],
					  neigh->step[0]);
		}	       
		
		eval_i = xi_opt[idx[ne] + typ2] * fnval_tail;
		if(typ1 == typ2) {
		  eval_j = eval_i;
		} else {
		  eval_j = xi_opt[idx[ne] + typ1] * fnval_tail;
		}
		fnval =  xi_opt[idx[ne] + typ1] * eval_i;
		grad =  xi_opt[idx[ne] + typ1] * xi_opt[idx[ne] + typ2] * grad_tail;
		
		if (self) {
		  eval_i *= 0.5;
		  eval_j *= 0.5;  
		  fnval *= 0.5;
		  grad *= 0.5;
		}
		
		forces[energy_p + h] += fnval;

		if (uf) {
		  tmp_force.x = neigh->dist.x * grad;
		  tmp_force.y = neigh->dist.y * grad;
		  tmp_force.z = neigh->dist.z * grad;
		  forces[k] += tmp_force.x;
		  forces[k + 1] += tmp_force.y;
		  forces[k + 2] += tmp_force.z;
		  /* actio = reactio */
		  l = 3 * neigh->nr;
		  forces[l] -= tmp_force.x;
		  forces[l + 1] -= tmp_force.y;
		  forces[l + 2] -= tmp_force.z;

#ifdef STRESS
		  /* calculate pair stresses */
		  if (us) {
		    tmp_force.x *= neigh->r;
		    tmp_force.y *= neigh->r;
		    tmp_force.z *= neigh->r;
		    stresses = stress_p + 6 * h;
		    forces[stresses] -= neigh->dist.x * tmp_force.x;
		    forces[stresses + 1] -= neigh->dist.y * tmp_force.y;
		    forces[stresses + 2] -= neigh->dist.z * tmp_force.z;
		    forces[stresses + 3] -= neigh->dist.x * tmp_force.y;
		    forces[stresses + 4] -= neigh->dist.y * tmp_force.z;
		    forces[stresses + 5] -= neigh->dist.z * tmp_force.x;
		  }
#endif /* STRESS */
		}

		/* calculate static field-contributions */
		atom->E_stat.x += neigh->dist.x * eval_i;
		atom->E_stat.y += neigh->dist.y * eval_i;
		atom->E_stat.z += neigh->dist.z * eval_i; 
		atoms[neigh->nr].E_stat.x += neigh->dist.x * eval_j;  
		atoms[neigh->nr].E_stat.y += neigh->dist.y * eval_j;  
		atoms[neigh->nr].E_stat.z += neigh->dist.z * eval_j;  

		/* calculate short-range dipoles */
		if( (xi_opt[idx[ne + ntypes] + typ1]) && (xi_opt[idx[ne + 2*ntypes] + neigh->col[0]]) && (xi_opt[idx[ne + 2*ntypes + size] + neigh->col[0]]) ) {	
		p_sr_tail = shortrange_value(neigh->r, &xi_opt[idx[ne + ntypes] + typ1], 
					       &xi_opt[idx[ne + 2*ntypes] + neigh->col[0]], &xi_opt[idx[ne + 2*ntypes + size] + neigh->col[0]]);
		atom->p_sr.x += xi_opt[idx[ne] + typ2] * neigh->dist.x * p_sr_tail;
		atom->p_sr.y += xi_opt[idx[ne] + typ2] * neigh->dist.y * p_sr_tail;
		atom->p_sr.z += xi_opt[idx[ne] + typ2] * neigh->dist.z * p_sr_tail;
		}
		if( (xi_opt[idx[ne + ntypes] + typ2]) && (xi_opt[idx[ne + 2*ntypes] + neigh->col[0]]) && (xi_opt[idx[ne + 2*ntypes + size] + neigh->col[0]]) && !self) {
		  p_sr_tail = shortrange_value(neigh->r, &xi_opt[idx[ne + ntypes] + typ2], 
						 &xi_opt[idx[ne + 2*ntypes] + neigh->col[0]], &xi_opt[idx[ne + 2*ntypes + size] + neigh->col[0]]);
		  atoms[neigh->nr].p_sr.x += xi_opt[idx[ne] + typ1] * neigh->dist.x * p_sr_tail;
		  atoms[neigh->nr].p_sr.y += xi_opt[idx[ne] + typ1] * neigh->dist.y * p_sr_tail;
		  atoms[neigh->nr].p_sr.z += xi_opt[idx[ne] + typ1] * neigh->dist.z * p_sr_tail;
		}
		
		}
	      
	  }			/* loop over neighbours */

	  /*then we can calculate contribution of forces right away */
	  if (uf) {
#ifdef FWEIGHT
	    /* Weigh by absolute value of force */
	    forces[k] /= FORCE_EPS + atom->absforce;
	    forces[k + 1] /= FORCE_EPS + atom->absforce;
	    forces[k + 2] /= FORCE_EPS + atom->absforce;
#endif /* FWEIGHT */
	    /* Returned force is difference between */
	    /* calculated and input force */
	    tmpsum +=
	      conf_weight[h] * (SQR(forces[k]) + SQR(forces[k + 1]) +
				SQR(forces[k + 2]));
	  }
	}         /* end S E C O N D loop over atoms */


	/* energy contributions from short-range and monopole interactions */
	forces[energy_p + h] *= eweight / (real)inconf[h];
	forces[energy_p + h] -= force_0[energy_p + h];
	tmpsum += conf_weight[h] * SQR(forces[energy_p + h]);
#ifdef STRESS
	/* stress contributions from short-range and monopole interactions */
	if (uf && us) {
	  for (i = 0; i < 6; i++) {
	    forces[stress_p + 6 * h + i] *= sweight / conf_vol[h - firstconf];
	    forces[stress_p + 6 * h + i] -= force_0[stress_p + 6 * h + i];
	    tmpsum += conf_weight[h] * SQR(forces[stress_p + 6 * h + i]);
	  }
	}
#endif /* STRESS */

	/* T H I R D loop: calculate whole dipole moment for every atom */
	real rp, dp_sum;
	int dp_converged = 0, dp_it = 0;
	real max_diff = 100000.;

	while(dp_converged == 0) {
	  dp_sum = 0;
	  for (i = 0; i < inconf[h]; i++) {    //atoms	  
	    atom = conf_atoms + i + cnfstart[h] - firstatom;
	    typ1 = atom->typ;
	
	      if(xi_opt[idx[ne + ntypes] + typ1]) {

		atom->p_ind.x = xi_opt[idx[ne + ntypes] + typ1] * (atom->E_stat.x + atom->E_ind.x) + atom->p_sr.x;
		atom->p_ind.y = xi_opt[idx[ne + ntypes] + typ1] * (atom->E_stat.y + atom->E_ind.y) + atom->p_sr.x;
		atom->p_ind.z = xi_opt[idx[ne + ntypes] + typ1] * (atom->E_stat.z + atom->E_ind.z) + atom->p_sr.x;
	      
		atom->E_temp.x = 0;
		atom->E_temp.y = 0;
		atom->E_temp.z = 0;

		atom->E_old.x = atom->E_ind.x;
		atom->E_old.y = atom->E_ind.y;
		atom->E_old.z = atom->E_ind.z;
	      } 
	  }
    
	  
	  for (i = 0; i < inconf[h]; i++) {     //atoms	  
	    atom = conf_atoms + i + cnfstart[h] - firstatom;
	    typ1 = atom->typ;
	    for (j = 0; j < atom->n_neigh; j++) {  //neighbors
	      neigh = atom->neigh + j;
	      typ2 = neigh->typ;
	      /* In small cells, an atom might interact with itself */
	      self = (neigh->nr == i + cnfstart[h]) ? 1 : 0;
	      
	      if (neigh->r < dp_cut) { 
		if(xi_opt[idx[ne + ntypes] + typ1]) {
		  rp = SPROD(atoms[neigh->nr].p_ind, neigh->dist);	      
		  atom->E_temp.x += dp_eps * ( 3 * rp * neigh->dist.x / neigh->r5 - atoms[neigh->nr].p_ind.x / neigh->r3);
		  atom->E_temp.y += dp_eps * ( 3 * rp * neigh->dist.y / neigh->r5 - atoms[neigh->nr].p_ind.y / neigh->r3);
		  atom->E_temp.z += dp_eps * ( 3 * rp * neigh->dist.z / neigh->r5 - atoms[neigh->nr].p_ind.z / neigh->r3);
		}
		if(xi_opt[idx[ne + ntypes] + typ2] && !self) {
		  rp = SPROD(atom->p_ind, neigh->dist);
		  atoms[neigh->nr].E_temp.x += dp_eps * ( 3 * rp * neigh->dist.x / neigh->r5 - atom->p_ind.x / neigh->r3);
		  atoms[neigh->nr].E_temp.y += dp_eps * ( 3 * rp * neigh->dist.y / neigh->r5 - atom->p_ind.y / neigh->r3);
		  atoms[neigh->nr].E_temp.z += dp_eps * ( 3 * rp * neigh->dist.z / neigh->r5 - atom->p_ind.z / neigh->r3);
		}
	      }
	    }
	  }

	  for (i = 0; i < inconf[h]; i++) {    //atoms	  
	    atom = conf_atoms + i + cnfstart[h] - firstatom;
	    typ1 = atom->typ;

	    if(xi_opt[idx[ne + ntypes] + typ1]) {

	      atom->E_ind.x = (1 - dp_mix) * atom->E_temp.x + dp_mix * atom->E_old.x;
	      atom->E_ind.y = (1 - dp_mix) * atom->E_temp.y + dp_mix * atom->E_old.y;
	      atom->E_ind.z = (1 - dp_mix) * atom->E_temp.z + dp_mix * atom->E_old.z;
	     
	      dp_sum += SQR(atom->E_old.x - atom->E_ind.x);
	      dp_sum += SQR(atom->E_old.y - atom->E_ind.y);
	      dp_sum += SQR(atom->E_old.z - atom->E_ind.z);
	    }
	  }
	  
	  dp_sum = sqrt(dp_sum);
	  dp_sum /= inconf[h];
	  if (dp_it) {
	    if( (dp_sum > max_diff) || (dp_it > 50) ) {
	      //  printf("Convergence error in dipole iteration loop %d\n", dp_it);
	      dp_converged = 1;
	      for (i = 0; i < inconf[h]; i++) {    //atoms	  
		atom = conf_atoms + i + cnfstart[h] - firstatom;
		typ1 = atom->typ;
		if(xi_opt[idx[ne + ntypes] + typ1]) {
		  atom->p_ind.x = xi_opt[idx[ne + ntypes] + typ1] * atom->E_stat.x + atom->p_sr.x;
		  atom->p_ind.y = xi_opt[idx[ne + ntypes] + typ1] * atom->E_stat.y + atom->p_sr.x;
		  atom->p_ind.z = xi_opt[idx[ne + ntypes] + typ1] * atom->E_stat.z + atom->p_sr.x;
		  atom->E_ind.x = atom->E_stat.x;
		  atom->E_ind.y = atom->E_stat.y;
		  atom->E_ind.z = atom->E_stat.z;
		}
	      }
	    }
	  } 

	  if (dp_sum < dp_tol)
	    dp_converged = 1;
	  
	  dp_it++;
	}        	/* end T H I R D loop over atoms */  
	
	

	/* F O U R T H  loop: calculate monopole-dipole and dipole-dipole forces */     
	for (i = 0; i < inconf[h]; i++) {    //atoms	  
	  atom = conf_atoms + i + cnfstart[h] - firstatom;
	  typ1 = atom->typ;

	  if(xi_opt[idx[ne + ntypes] + typ1]) {
	   
	    fnval =  - SPROD(atom->p_ind, atom->E_ind);
	    forces[energy_p + h] += fnval;

		if (uf) {
		  tmp_force.x = - atom->E_ind.x;
		  tmp_force.y = - atom->E_ind.y;
		  tmp_force.z = - atom->E_ind.z;
		  forces[k] += tmp_force.x;
		  forces[k + 1] += tmp_force.y;
		  forces[k + 2] += tmp_force.z;
		  /* actio = reactio */
		  l = 3 * neigh->nr;
		  forces[l] -= tmp_force.x;
		  forces[l + 1] -= tmp_force.y;
		  forces[l + 2] -= tmp_force.z;

#ifdef STRESS
		  /* calculate pair stresses */
		  if (us) {
		    tmp_force.x *= neigh->r;
		    tmp_force.y *= neigh->r;
		    tmp_force.z *= neigh->r;
		    stresses = stress_p + 6 * h;
		    forces[stresses] -= neigh->dist.x * tmp_force.x;
		    forces[stresses + 1] -= neigh->dist.y * tmp_force.y;
		    forces[stresses + 2] -= neigh->dist.z * tmp_force.z;
		    forces[stresses + 3] -= neigh->dist.x * tmp_force.y;
		    forces[stresses + 4] -= neigh->dist.y * tmp_force.z;
		    forces[stresses + 5] -= neigh->dist.z * tmp_force.x;
		  }
#endif /* STRESS */
		}
 
	  }
	}  /* end F O U R T H loop over atoms */ 


      }				/* end M A I N loop over configurations */
    }				/* parallel region */
 
    /* dummy constraints (global) */
#ifdef APOT
    /* add punishment for out of bounds (mostly for powell_lsq) */
    if (myid == 0) {
      tmpsum += apot_punish(xi_opt, forces);
    }
#endif

    sum = tmpsum;		/* global sum = local sum  */

#ifdef MPI
    /* reduce global sum */
    sum = 0.;
    MPI_Reduce(&tmpsum, &sum, 1, REAL, MPI_SUM, 0, MPI_COMM_WORLD);
    /* gather forces, energies, stresses */
    MPI_Gatherv(forces + firstatom * 3, myatoms, MPI_VEKTOR,	/* forces */
		forces, atom_len, atom_dist, MPI_VEKTOR, 0, MPI_COMM_WORLD);
    MPI_Gatherv(forces + natoms * 3 + firstconf, myconf, REAL,	/* energies */
		forces
		+ natoms * 3, conf_len, conf_dist, REAL, 0, MPI_COMM_WORLD);
    /* stresses */
    MPI_Gatherv(forces + natoms * 3 + nconf +
		6 * firstconf, myconf, MPI_STENS,
		forces + natoms * 3 + nconf,
		conf_len, conf_dist, MPI_STENS, 0, MPI_COMM_WORLD);
#endif /* MPI */

    /* root process exits this function now */
    if (myid == 0) {
      fcalls++;			/* Increase function call counter */
      if (isnan(sum)) {
#ifdef DEBUG
	printf("\n--> Force is nan! <--\n\n");
#endif
	return 10e10;
      } else
	return sum;
    }

  }

  /* once a non-root process arrives here, all is done. */
  return -1.;
}

#endif /* DIPOLE */
