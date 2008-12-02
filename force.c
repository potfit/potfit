/****************************************************************
* 
* force.c: Routines used for calculating forces/energies in various 
*     interpolation schemes. 
*
*****************************************************************/
/*
*   Copyright 2002-2008 Peter Brommer, Franz G"ahler, Daniel Schopf
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
* $Revision: 1.55 $
* $Date: 2008/12/02 15:19:53 $
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

real calc_forces_pair(real *xi_opt, real *forces, int flag)
{
  real  tmpsum, sum = 0.;
  int   first, col1;
  real  grad0, y0, y1, x0, x1;
  real *xi;
#ifdef EAM
  static real rho_sum_loc, rho_sum;
  rho_sum_loc = rho_sum = 0.;
#endif /* EAM */
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
  /* This is the start of an infinite loop */
  while (1) {
    tmpsum = 0.;		/* sum of squares of local process */
    if ((format == 0 || format > 4) && myid == 0)
#ifdef APOT
      update_calc_table(xi_opt, xi, 0);
#else
      update_calc_table(xi_opt, xi);
#endif

#ifdef MPI
    /* exchange potential and flag value */
    MPI_Bcast(xi, calc_pot.len, REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef APOT
    if (myid == 0 && do_smooth)
      potsync_apot();
    if (flag == 1)
      break;
#endif

#ifdef EAM
    /* if flag==2 then the potential parameters have changed -> sync */
    if (flag == 2)
      potsync();
    /* non root processes hang on, unless...  */
    if (flag == 1)
      break;			/* Exception: flag 1 means clean up */
#endif /* EAM || APOT */

#endif /* MPI */

    /* init second derivatives for splines */
    for (col1 = 0; col1 < paircol; col1++) {	/* just pair potentials */
      first = calc_pot.first[col1];
      if (format == 3 || format == 0)
	spline_ed(calc_pot.step[col1], xi + first,
		  calc_pot.last[col1] - first + 1,
		  *(xi + first - 2), 0.0, calc_pot.d2tab + first);
      else			/* format == 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
		  calc_pot.last[col1] - first + 1,
		  *(xi + first - 2), 0.0, calc_pot.d2tab + first);
    }

#ifdef EAM
    for (col1 = paircol; col1 < paircol + ntypes; col1++) {	/* rho */
      first = calc_pot.first[col1];
      if (format == 3)
	spline_ed(calc_pot.step[col1], xi + first,
		  calc_pot.last[col1] - first + 1,
		  *(xi + first - 2), 0.0, calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
		  calc_pot.last[col1] - first + 1,
		  *(xi + first - 2), 0.0, calc_pot.d2tab + first);
    }
#ifndef PARABEL
/* if we have parabolic interpolation, we don't need that */
    for (col1 = paircol + ntypes; col1 < paircol + 2 * ntypes; col1++) {
      /* F */
      first = calc_pot.first[col1];
      /* gradient at left boundary matched to square root function, 
         when 0 not in domain(F), else natural spline */
      if (format == 3 || format == 0)
	spline_ed(calc_pot.step[col1], xi + first,
		  calc_pot.last[col1] - first + 1,
#ifdef WZERO
		  ((calc_pot.begin[col1] <= 0.) ? *(xi + first - 2)
		   : .5 / xi[first]),
		  ((calc_pot.end[col1] >= 0.) ? *(xi + first - 1)
		   : -.5 / xi[calc_pot.last[col1]]),
#else /* WZERO: F is natural spline in any case */
		  *(xi + first - 2), *(xi + first - 1),
#endif /* WZERO */
		  calc_pot.d2tab + first);	/* XXX */
      else			/* format == 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
		  calc_pot.last[col1] - first + 1,
#ifdef WZERO
		  (calc_pot.begin[col1] <= 0. ? *(xi + first - 2)
		   : .5 / xi[first]),
		  (calc_pot.end[col1] >= 0. ? *(xi + first - 1)
		   : -.5 / xi[calc_pot.last[col1]]),
#else /* WZERO */
		  *(xi + first - 2), *(xi + first - 1),
#endif /* WZERO */
		  calc_pot.d2tab + first);
    }
#endif /* PARABEL */
#endif /* EAM */

#ifndef MPI
    myconf = nconf;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif

    /* region containing loop over configurations, 
       also OMP-parallelized region */
    {

#ifdef EAM
      vektor d;
      int   col2;
      real  gradF, gradF2, grad2, r, eamforce;
      atom_t *atom2;
#endif
      vektor tmp_force;
      int   h, k, i, l, j, typ1, typ2, col, config, stresses, uf, us;
      real  fnval, grad;
      atom_t *atom;

      neigh_t *neigh;


#ifdef _OPENMP
#pragma omp for reduction(+:tmpsum,rho_sum_loc)
#endif
      /* loop over configurations */
      for (h = firstconf; h < firstconf + myconf; h++) {
	config = 3 * natoms + h;	/* slot for energies */
	stresses = 3 * natoms + nconf + 6 * h;	/* slot for stresses */
	uf = conf_uf[h - firstconf];
	us = conf_us[h - firstconf];
	/* reset forces */
	forces[config] = 0.;
	for (i = stresses; i < stresses + 6; i++)
	  forces[i] = 0.;

#ifdef EAM
	/* set dummy constraints */
	forces[config + 7 * nconf] = -force_0[config + 7 * nconf];
#endif
	/* first loop over atoms: reset forces, densities */
	for (i = 0; i < inconf[h]; i++) {
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

#ifdef EAM
	  conf_atoms[cnfstart[h] - firstatom + i].rho = 0.0;
#endif
	}			/* end first loop */

	/* 2nd loop: calculate pair forces and energies, atomic densities. */
	for (i = 0; i < inconf[h]; i++) {
	  atom = conf_atoms + i + cnfstart[h] - firstatom;
	  typ1 = atom->typ;
	  k = 3 * (cnfstart[h] + i);

	  /* loop over neighbours */
	  for (j = 0; j < atom->n_neigh; j++) {
	    neigh = atom->neigh + j;
	    /* only use neigbours with higher numbers, 
	       others are calculated by actio=reactio */
	    if (neigh->nr > i + cnfstart[h]) {
	      typ2 = neigh->typ;
	      /* find correct column */
	      col = (typ1 <= typ2) ?
		typ1 * ntypes + typ2 - ((typ1 * (typ1 + 1)) / 2)
		: typ2 * ntypes + typ1 - ((typ2 * (typ2 + 1)) / 2);
	      if (neigh->r < calc_pot.end[col]) {
		/* fn value and grad are calculated in the same step */

		if (uf)
		  fnval = splint_comb_dir(&calc_pot, xi, col,
					  neigh->slot[0], neigh->shift[0],
					  neigh->step[0], &grad);
		else
		  fnval = splint_dir(&calc_pot, xi, col,
				     neigh->slot[0], neigh->shift[0],
				     neigh->step[0]);
		forces[config] += fnval;
/* not real force: cohesive energy */
		if (uf) {
		  tmp_force.x = neigh->dist.x * grad;
		  tmp_force.y = neigh->dist.y * grad;
		  tmp_force.z = neigh->dist.z * grad;
		  forces[k] += tmp_force.x;
		  forces[k + 1] += tmp_force.y;
		  forces[k + 2] += tmp_force.z;
		  l = 3 * neigh->nr;	/* actio = reactio */
		  forces[l] -= tmp_force.x;
		  forces[l + 1] -= tmp_force.y;
		  forces[l + 2] -= tmp_force.z;
#ifdef STRESS
		  /* also calculate pair stresses */
		  if (us) {
		    tmp_force.x *= neigh->r;
		    tmp_force.y *= neigh->r;
		    tmp_force.z *= neigh->r;
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
#ifdef EAM
	      /* calculate atomic densities */
	      col2 = paircol + typ2;
	      if (typ2 == typ1) {
/* then transfer(a->b)==transfer(b->a) */
		if (neigh->r < calc_pot.end[col2]) {
		  fnval = splint_dir(&calc_pot, xi, col2,
				     neigh->slot[1], neigh->shift[1],
				     neigh->step[1]);
		  atom->rho += fnval;
		  conf_atoms[neigh->nr - firstatom].rho += fnval;
		}
	      } else {		/* transfer(a->b)!=transfer(b->a) */
		col = paircol + typ1;
		if (neigh->r < calc_pot.end[col2]) {
		  atom->rho += splint_dir(&calc_pot, xi, col2,
					  neigh->slot[1], neigh->shift[1],
					  neigh->step[1]);
		}
		/* cannot use slot/shift to access splines */
		if (neigh->r < calc_pot.end[col])
		  conf_atoms[neigh->nr - firstatom].rho +=
		    splint(&calc_pot, xi, col, neigh->r);
	      }

#endif /* EAM */
	    }			/*  neighbours with bigger atom nr */
	  }			/* loop over neighbours */
#ifndef EAM
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
	      SQR(forces[k]) + SQR(forces[k + 1]) + SQR(forces[k + 2]);
	  }
#else /* EAM */

	  col2 = paircol + ntypes + typ1;	/* column of F */
#ifndef NORESCALE
	  if (atom->rho > calc_pot.end[col2]) {
	    /* then punish target function -> bad potential */
	    forces[config + 7 * nconf] +=
	      1000 * SQR(atom->rho - calc_pot.end[col2]);
#ifndef PARABEL
/* then we use the final value, with PARABEL: extrapolate */
	    atom->rho = calc_pot.end[col2];
#endif /* PARABEL */
	  }

	  if (atom->rho < calc_pot.begin[col2]) {
	    /* then punish target function -> bad potential */
	    forces[config + 7 * nconf] +=
	      1000 * SQR(calc_pot.begin[col2] - atom->rho);
#ifndef PARABEL
/* then we use the final value, with PARABEL: extrapolate */
	    atom->rho = calc_pot.begin[col2];
#endif /* PARABEL */
	  }
#endif /* NOT NORESCALE */
	  /* embedding energy, embedding gradient */
	  /* contribution to cohesive energy is F(n) */
#ifdef PARABEL
	  forces[config] +=
	    parab_comb(&calc_pot, xi, col2, atom->rho, &atom->gradF);
#elif defined(NORESCALE)
	  if (atom->rho < calc_pot.begin[col2]) {
	    /* linear extrapolation left */
	    fnval =
	      splint_comb(&calc_pot, xi, col2, calc_pot.begin[col2],
			  &atom->gradF);
	    forces[config] +=
	      fnval + (atom->rho - calc_pot.begin[col2]) * atom->gradF;
	  } else if (atom->rho > calc_pot.end[col2]) {
	    /* and right */
	    fnval =
	      splint_comb(&calc_pot, xi, col2, calc_pot.end[col2],
			  &atom->gradF);
	    forces[config] +=
	      fnval + (atom->rho - calc_pot.end[col2]) * atom->gradF;
	  }
	  /* and in-between */
	  else {
	    forces[config] +=
	      splint_comb(&calc_pot, xi, col2, atom->rho, &atom->gradF);
	  }
#else
	  forces[config] +=
	    splint_comb(&calc_pot, xi, col2, atom->rho, &atom->gradF);
#endif
	  /* sum up rho */
	  rho_sum_loc += atom->rho;
#endif /* EAM */
	}			/* second loop over atoms */

#ifdef EAM			/* if we don't have EAM, we're done */
	/* 3rd loop over atom: EAM force */
	if (uf) {		/* only required if we calc forces */
	  for (i = 0; i < inconf[h]; i++) {
	    atom = conf_atoms + i + cnfstart[h] - firstatom;
	    typ1 = atom->typ;
	    k = 3 * (cnfstart[h] + i);
	    col = paircol + ntypes + typ1;	/* column of F */


	    for (j = 0; j < atom->n_neigh; j++) {
	      /* loop over neighbours */
	      neigh = atom->neigh + j;
	      /* only neigbours higher than current atom are of interest */
	      if (neigh->nr > i + cnfstart[h]) {
		typ2 = neigh->typ;
		col2 = paircol + typ2;
		r = neigh->r;

		/* are we within reach? */
		if ((r < calc_pot.end[col2])
		    || (r < calc_pot.end[col - ntypes])) {
		  grad =
		    (r < calc_pot.end[col2]) ?
		    splint_grad_dir(&calc_pot, xi, col2, neigh->slot[1],
				    neigh->shift[1], neigh->step[1])
		    : 0.;

		  if (typ2 == typ1)	/* use actio = reactio */
		    grad2 = grad;
		  else
		    grad2 = (r < calc_pot.end[col - ntypes]) ?
		      splint_grad(&calc_pot, xi, col - ntypes, r) : 0.;

		  /* now we know everything - calculate forces */
		  eamforce = (grad * atom->gradF +
			      grad2 *
			      conf_atoms[(neigh->nr) - firstatom].gradF);

		  tmp_force.x = neigh->dist.x * eamforce;
		  tmp_force.y = neigh->dist.y * eamforce;
		  tmp_force.z = neigh->dist.z * eamforce;
		  forces[k] += tmp_force.x;
		  forces[k + 1] += tmp_force.y;
		  forces[k + 2] += tmp_force.z;
		  l = 3 * neigh->nr;
		  forces[l] -= tmp_force.x;
		  forces[l + 1] -= tmp_force.y;
		  forces[l + 2] -= tmp_force.z;
#ifdef STRESS
		  /* and stresses */
		  if (us) {
		    tmp_force.x *= neigh->r;
		    tmp_force.y *= neigh->r;
		    tmp_force.z *= neigh->r;
		    forces[stresses] -= neigh->dist.x * tmp_force.x;
		    forces[stresses + 1] -= neigh->dist.y * tmp_force.y;
		    forces[stresses + 2] -= neigh->dist.z * tmp_force.z;
		    forces[stresses + 3] -= neigh->dist.x * tmp_force.y;
		    forces[stresses + 4] -= neigh->dist.y * tmp_force.z;
		    forces[stresses + 5] -= neigh->dist.z * tmp_force.x;
		  }
#endif /* STRESS */
		}		/* within reach */
	      }			/* higher neigbours */
	    }			/* loop over neighbours */
#ifdef FWEIGHT
	    /* Weigh by absolute value of force */
	    forces[k] /= FORCE_EPS + atom->absforce;
	    forces[k + 1] /= FORCE_EPS + atom->absforce;
	    forces[k + 2] /= FORCE_EPS + atom->absforce;
#endif /* FWEIGHT */

	    /* sum up forces  */
	    tmpsum +=
	      SQR(forces[k]) + SQR(forces[k + 1]) + SQR(forces[k + 2]);
	  }			/* third loop over atoms */
	}			/* use forces */
#endif /* EAM */

	/* energy contributions */
	forces[config] *= eweight / (real)inconf[h];
	forces[config] -= force_0[config];
	tmpsum += SQR(forces[config]);

#ifdef STRESS
	/* stress contributions */
	if (uf && us) {
	  for (i = stresses; i < stresses + 6; i++) {
	    forces[i] *= sweight / conf_vol[h - firstconf];
	    forces[i] -= force_0[i];
	    tmpsum += SQR(forces[i]);

	  }
	}
#endif /* STRESS */
	/* dummy constraints per configuration */
#ifdef EAM
	tmpsum += SQR(forces[config + 7 * nconf]);
#endif
      }				/* loop over configurations */
    }				/* parallel region */
#ifdef MPI
    /* Reduce rho_sum */
    MPI_Reduce(&rho_sum_loc, &rho_sum, 1, REAL, MPI_SUM, 0, MPI_COMM_WORLD);
#else /* MPI */
#ifdef EAM
    rho_sum = rho_sum_loc;
#endif /* EAM */
#endif /* MPI */
    /* dummy constraints (global) */
#ifdef APOT
    int   i, j;
    real  x;
    /* *INDENT-OFF* */ 
    /* turn off indent because of a bug in GNU indent 2.2.10 */
    if (myid == 0) {
      for (i = 0; i < ndim; i++) {
	forces[mdim - apot_table.total_par + i] = 0;
	if (x = xi_opt[idx[i]] -
	    apot_table.pmin[apot_table.idxpot[i]][apot_table.idxparam[i]],
	    x <= 0) {
	  forces[mdim - apot_table.total_par + i] = APOT_PUNISH * x * x;
	} else if (x = xi_opt[idx[i]] -
		   apot_table.pmax[apot_table.
				   idxpot[i]][apot_table.idxparam[i]],
		   x >= 0) {
	  forces[mdim - apot_table.total_par + i] = APOT_PUNISH * x * x;
	}
	tmpsum += forces[mdim - apot_table.total_par + i];
      }
    }
    /* *INDENT-ON* */
#endif
#ifdef EAM
    if (myid == 0) {
      int   g, i;
      for (g = 0; g < ntypes; g++) {
	/* PARABEL, WZERO, NORESC - different behaviour */
#ifdef PARABEL
/* constraints on U(n) */
	forces[mdim - ntypes + g] = DUMMY_WEIGHT *
	  parab(&calc_pot, xi, paircol + ntypes + g, 0.)
	  - force_0[mdim - ntypes + g];
/* constraints on U`(n) */
	forces[mdim - 2 * ntypes + g] = DUMMY_WEIGHT *
	  parab_grad(&calc_pot, xi, paircol + ntypes + g,
		     .5 * (calc_pot.begin[paircol + ntypes + g] +
			   calc_pot.end[paircol + ntypes + g]))
	  - force_0[mdim - 2 * ntypes + g];
#elif defined(WZERO)
	if (calc_pot.begin[paircol + ntypes + g] <= 0.)
	  /* 0 in domain of U(n) */
/* constraints on U(n) */
	  forces[mdim - ntypes + g] = DUMMY_WEIGHT *
	    splint(&calc_pot, xi, paircol + ntypes + g, 0.)
	    - force_0[mdim - ntypes + g];
	else
	  /* 0 not in domain of U(n) */
	  forces[mdim - ntypes + g] = 0.;	/* Free end... */
/* constraints on U`(n) */
	forces[mdim - 2 * ntypes + g] = DUMMY_WEIGHT *
	  splint_grad(&calc_pot, xi, paircol + ntypes + g,
		      .5 * (calc_pot.begin[paircol + ntypes + g] +
			    calc_pot.end[paircol + ntypes + g]))
	  - force_0[mdim - 2 * ntypes + g];
#elif defined(NORESCALE)
	/* clear field */
	forces[mdim - ntypes + g] = 0.;	/* Free end... */
	/* NEW: Constraint on U': U'(1.)=0; */
	forces[mdim - 2 * ntypes + g] = DUMMY_WEIGHT *
	  splint_grad(&calc_pot, xi, paircol + ntypes + g, 1.);
#else /* NOTHING */
	forces[mdim - ntypes + g] = 0.;	/* Free end... */
/* constraints on U`(n) */
	forces[mdim - 2 * ntypes + g] = DUMMY_WEIGHT *
	  splint_grad(&calc_pot, xi, paircol + ntypes + g,
		      .5 * (calc_pot.begin[paircol + ntypes + g] +
			    calc_pot.end[paircol + ntypes + g]))
	  - force_0[mdim - 2 * ntypes + g];

#endif /* Dummy constraints */
	tmpsum += SQR(forces[mdim - ntypes + g]);
	tmpsum += SQR(forces[mdim - 2 * ntypes + g]);
      }				/* loop over types */
#ifdef NORESCALE
      /* NEW: Constraint on n: <n>=1. ONE CONSTRAINT ONLY */
      /* Calculate averages */
      rho_sum /= (real)natoms;
      /* ATTN: if there are invariant potentials, things might be problematic */
      forces[mdim - ntypes] = DUMMY_WEIGHT * (rho_sum - 1.);
      tmpsum += SQR(forces[mdim - ntypes]);
#endif
    }				/* only root process */
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
		forces + natoms * 3, conf_len, conf_dist, REAL, 0,
		MPI_COMM_WORLD);
    /* stresses */
    MPI_Gatherv(forces + natoms * 3 + nconf + 6 * firstconf, myconf,
		MPI_STENS, forces + natoms * 3 + nconf, conf_len, conf_dist,
		MPI_STENS, 0, MPI_COMM_WORLD);
#ifdef EAM
    /* punishment constraints */
    MPI_Gatherv(forces + natoms * 3 + 7 * nconf + firstconf, myconf, REAL,
		forces + natoms * 3 + 7 * nconf, conf_len, conf_dist, REAL,
		0, MPI_COMM_WORLD);
    /* no need to pick up dummy constraints - are already @ root */
#endif /* EAM */
#endif /* MPI */

    /* root process exits this function now */
    if (myid == 0) {
      fcalls++;			/* Increase function call counter */
      return sum;
    }

  }
  /* once a non-root process arrives here, all is done. */
  return -1.;
}
