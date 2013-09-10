/****************************************************************
 *
 * force_eam.c: Routine used for calculating eam forces/energies
 *
 ****************************************************************
 *
 * Copyright 2002-2013
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

#include "potfit.h"

#if defined EAM && !defined COULOMB

#include "functions.h"
#include "potential.h"
#include "splines.h"
#include "utils.h"

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

double calc_forces_eam(double *xi_opt, double *forces, int flag)
{
  int   first, col, i;
  double tmpsum = 0.0, sum = 0.0;
  double *xi = NULL;

  static double rho_sum_loc, rho_sum;
  rho_sum_loc = rho_sum = 0.0;

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
    tmpsum = 0.0;		/* sum of squares of local process */
    rho_sum_loc = 0.0;

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
#ifndef PARABOLA
    /* if we have parabolic interpolation, we don't need that */
    for (col = paircol + ntypes; col < paircol + 2 * ntypes; col++) {
      first = calc_pot.first[col];
      /* gradient at left boundary matched to square root function,
         when 0 not in domain(F), else natural spline */
      if (0 == format || 3 == format)
	spline_ed(calc_pot.step[col], xi + first, calc_pot.last[col] - first + 1,
#ifdef WZERO
	  ((calc_pot.begin[col] <= 0.0) ? *(xi + first - 2)
	    : 0.5 / xi[first]), ((calc_pot.end[col] >= 0.0) ? *(xi + first - 1)
	    : -0.5 / xi[calc_pot.last[col]]),
#else /* WZERO: F is natural spline in any case */
	  *(xi + first - 2), *(xi + first - 1),
#endif /* WZERO */
	  calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first, calc_pot.last[col] - first + 1,
#ifdef WZERO
	  (calc_pot.begin[col] <= 0.0 ? *(xi + first - 2)
	    : 0.5 / xi[first]), (calc_pot.end[col] >= 0.0 ? *(xi + first - 1)
	    : -0.5 / xi[calc_pot.last[col]]),
#else /* WZERO */
	  *(xi + first - 2), *(xi + first - 1),
#endif /* WZERO */
	  calc_pot.d2tab + first);
    }
#endif /* PARABOLA */

#ifndef MPI
    myconf = nconf;
#endif /* MPI */

    /* region containing loop over configurations */
    {
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
      double phi_val, phi_grad;
      double r;
      vector tmp_force;

      /* eam variables */
      int   col_F;
      double eam_force;
      double rho_val, rho_grad, rho_grad_j;

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

	/* set limiting constraints */
	forces[limit_p + h] = -force_0[limit_p + h];

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
	      }
	    }

	    /* neighbor in range */
	    /* calculate atomic densities */
	    if (atom->type == neigh->type) {
	      /* then transfer(a->b)==transfer(b->a) */
	      if (neigh->r < calc_pot.end[neigh->col[1]]) {
		rho_val = splint_dir(&calc_pot, xi, neigh->slot[1], neigh->shift[1], neigh->step[1]);
		atom->rho += rho_val;
		/* avoid double counting if atom is interacting with a
		   copy of itself */
		if (!self) {
		  conf_atoms[neigh->nr - firstatom].rho += rho_val;
		}
	      }
	    } else {
	      /* transfer(a->b)!=transfer(b->a) */
	      if (neigh->r < calc_pot.end[neigh->col[1]]) {
		atom->rho += splint_dir(&calc_pot, xi, neigh->slot[1], neigh->shift[1], neigh->step[1]);
	      }
	      /* cannot use slot/shift to access splines */
	      if (neigh->r < calc_pot.end[paircol + atom->type])
		conf_atoms[neigh->nr - firstatom].rho +=
		  splint(&calc_pot, xi, paircol + atom->type, neigh->r);
	    }
	  }			/* loop over all neighbors */

	  col_F = paircol + ntypes + atom->type;	/* column of F */
#ifndef NORESCALE
	  if (atom->rho > calc_pot.end[col_F]) {
	    /* then punish target function -> bad potential */
	    forces[limit_p + h] += DUMMY_WEIGHT * 10.0 * dsquare(atom->rho - calc_pot.end[col_F]);
#ifndef PARABOLA
	    /* then we use the final value, with PARABOLA: extrapolate */
	    atom->rho = calc_pot.end[col_F];
#endif /* PARABOLA */
	  }

	  if (atom->rho < calc_pot.begin[col_F]) {
	    /* then punish target function -> bad potential */
	    forces[limit_p + h] += DUMMY_WEIGHT * 10.0 * dsquare(calc_pot.begin[col_F] - atom->rho);
#ifndef PARABOLA
	    /* then we use the final value, with PARABOLA: extrapolate */
	    atom->rho = calc_pot.begin[col_F];
#endif /* PARABOLA */
	  }
#endif /* !NORESCALE */

	  /* embedding energy, embedding gradient */
	  /* contribution to cohesive energy is F(n) */

#ifdef PARABOLA
	  forces[energy_p + h] += parab_comb(&calc_pot, xi, col_F, atom->rho, &atom->gradF);
#elif defined(NORESCALE)
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
	  }
	  /* and in-between */
	  else {
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
#endif /* NORESCALE */
	  /* sum up rho */
	  rho_sum_loc += atom->rho;
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
		    (r < calc_pot.end[col_F - ntypes]) ? splint_grad(&calc_pot, xi, col_F - ntypes, r) : 0.;
		/* now we know everything - calculate forces */
		eam_force = (rho_grad * atom->gradF + rho_grad_j * conf_atoms[(neigh->nr) - firstatom].gradF);
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
	/* limiting constraints per configuration */
	tmpsum += conf_weight[h] * dsquare(forces[limit_p + h]);
      }				/* loop over configurations */
    }				/* parallel region */
#ifdef MPI
    /* Reduce rho_sum */
    rho_sum = 0.0;
    MPI_Reduce(&rho_sum_loc, &rho_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else /* MPI */
    rho_sum = rho_sum_loc;
#endif /* MPI */

    /* dummy constraints (global) */
#ifdef APOT
    /* add punishment for out of bounds (mostly for powell_lsq) */
    if (0 == myid) {
      tmpsum += apot_punish(xi_opt, forces);
    }
#endif /* APOT */

#ifndef NOPUNISH
    if (0 == myid) {
      int   g;
      for (g = 0; g < ntypes; g++) {
	/* PARABOLA, WZERO, NORESC - different behaviour */
#ifdef PARABOLA
/* constraints on U(n) */
	forces[dummy_p + ntypes + g] = DUMMY_WEIGHT * parab(&calc_pot, xi, paircol + ntypes + g, 0.0)
	  - force_0[dummy_p + ntypes + g];
/* constraints on U`(n) */
	forces[dummy_p + g] =
	  DUMMY_WEIGHT * parab_grad(&calc_pot, xi, paircol + ntypes + g,
	  .5 * (calc_pot.begin[paircol + ntypes + g] + calc_pot.end[paircol + ntypes + g])) -
	  force_0[dummy_p + g];
#elif defined(WZERO)
	if (calc_pot.begin[paircol + ntypes + g] <= 0.0)
	  /* 0 in domain of U(n) */
/* constraints on U(n) */
	  forces[dummy_p + ntypes + g] = DUMMY_WEIGHT * splint(&calc_pot, xi, paircol + ntypes + g, 0.0)
	    - force_0[dummy_p + ntypes + g];
	else
	  /* 0 not in domain of U(n) */
	  forces[dummy_p + ntypes + g] = 0.0;	/* Free end... */
/* constraints on U`(n) */
	forces[dummy_p + g] =
	  DUMMY_WEIGHT * splint_grad(&calc_pot, xi, paircol + ntypes + g,
	  0.5 * (calc_pot.begin[paircol + ntypes + g] + calc_pot.end[paircol + ntypes + g]))
	  - force_0[dummy_p + g];
#elif defined(NORESCALE)
	/* clear field */
	forces[dummy_p + ntypes + g] = 0.0;	/* Free end... */
	/* NEW: Constraint on U': U'(1.)=0; */
	forces[dummy_p + g] = DUMMY_WEIGHT * splint_grad(&calc_pot, xi, paircol + ntypes + g, 1.0);
#else /* NOTHING */
	forces[dummy_p + ntypes + g] = 0.0;	/* Free end... */
/* constraints on U`(n) */
	forces[dummy_p + g] =
	  DUMMY_WEIGHT * splint_grad(&calc_pot, xi, paircol + ntypes + g,
	  0.5 * (calc_pot.begin[paircol + ntypes + g] + calc_pot.end[paircol + ntypes + g]))
	  - force_0[dummy_p + g];
#endif /* Dummy constraints */
	tmpsum += dsquare(forces[dummy_p + ntypes + g]);
	tmpsum += dsquare(forces[dummy_p + g]);
      }				/* loop over types */
#ifdef NORESCALE
      /* NEW: Constraint on n: <n>=1. ONE CONSTRAINT ONLY */
      /* Calculate averages */
      rho_sum /= (double)natoms;
      /* ATTN: if there are invariant potentials, things might be problematic */
      forces[dummy_p + ntypes] = DUMMY_WEIGHT * (rho_sum - 1.0);
      tmpsum += dsquare(forces[dummy_p + ntypes]);
#endif /* NORESCALE */
    }				/* only root process */
#endif /* !NOPUNISH */

#ifdef MPI
    /* reduce global sum */
    sum = 0.0;
    MPI_Reduce(&tmpsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    /* gather forces, energies, stresses */
    if (0 == myid) {		/* root node already has data in place */
      /* forces */
      MPI_Gatherv(MPI_IN_PLACE, myatoms, MPI_VECTOR, forces, atom_len,
	atom_dist, MPI_VECTOR, 0, MPI_COMM_WORLD);
      /* energies */
      MPI_Gatherv(MPI_IN_PLACE, myconf, MPI_DOUBLE, forces + natoms * 3,
	conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      /* stresses */
      MPI_Gatherv(MPI_IN_PLACE, myconf, MPI_STENS, forces + natoms * 3 + nconf,
	conf_len, conf_dist, MPI_STENS, 0, MPI_COMM_WORLD);
      /* punishment constraints */
      MPI_Gatherv(MPI_IN_PLACE, myconf, MPI_DOUBLE, forces + natoms * 3 + 7 * nconf,
	conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
      /* forces */
      MPI_Gatherv(forces + firstatom * 3, myatoms, MPI_VECTOR, forces, atom_len,
	atom_dist, MPI_VECTOR, 0, MPI_COMM_WORLD);
      /* energies */
      MPI_Gatherv(forces + natoms * 3 + firstconf, myconf, MPI_DOUBLE,
	forces + natoms * 3, conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      /* stresses */
      MPI_Gatherv(forces + natoms * 3 + nconf + 6 * firstconf, myconf, MPI_STENS,
	forces + natoms * 3 + nconf, conf_len, conf_dist, MPI_STENS, 0, MPI_COMM_WORLD);
      /* punishment constraints */
      MPI_Gatherv(forces + natoms * 3 + 7 * nconf + firstconf, myconf, MPI_DOUBLE,
	forces + natoms * 3 + 7 * nconf, conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    /* no need to pick up dummy constraints - are already @ root */
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

#endif /* EAM */
