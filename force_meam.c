/****************************************************************
 *
 * force_meam.c: Routine used for calculating meam forces/energies
 *
 ****************************************************************
 *
 * Copyright 2002-2013
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.sourceforge.net/
 *
 *	MEAM potential: Jeremy Nicklas
 *		Ohio State University, Physics Department
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

#ifdef MEAM

#include "potfit.h"

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

/****************************************************************
 *
 * DOES NOT WORK WITH THE FOLLOWING FLAGS:
 * PARABOLA, WZERO
 *
 ****************************************************************/

double calc_forces_meam(double *xi_opt, double *forces, int flag)
{
  int   first, col, i;
  double *xi = NULL;

  /* Some useful temp variables */
  static double tmpsum = 0.0, sum = 0.0;
  static double rho_sum = 0.0, rho_sum_loc = 0.0;

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

    /* Reset tmpsum and rho_sum_loc
       tmpsum = Sum of all the forces, energies and constraints
       rho_sum_loc = Sum of density, rho, for all atoms */
    tmpsum = 0.;
    rho_sum_loc = 0.;

#if defined APOT && !defined MPI
    if (0 == format) {
      apot_check_params(xi_opt);
      update_calc_table(xi_opt, xi, 0);
    }
#endif /* APOT && !MPI */

#ifdef MPI
    /* exchange potential and flag value */
#ifndef APOT
    MPI_Bcast(xi, calc_pot.len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif /* APOT */
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (flag == 1)
      break;			/* Exception: flag 1 means clean up */

#ifdef APOT
    if (myid == 0)
      apot_check_params(xi_opt);
    MPI_Bcast(xi_opt, ndimtot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    update_calc_table(xi_opt, xi, 0);
#else
    /* if flag==2 then the potential parameters have changed -> sync */
    if (flag == 2)
      potsync();
#endif /* APOT */
#endif /* MPI */

    /* First step is to initialize 2nd derivatives for splines */

    /* Pair potential (phi), density (rho), embedding funtion (F)
       where paircol is number of pair potential columns
       and ntypes is number of rho columns
       and ntypes is number of F columns */
    for (col = 0; col < 2 * paircol + 3 * ntypes; col++) {
      /* Pointer to first entry */
      first = calc_pot.first[col];

      /* Initialize 2nd derivatives
         step = width of spline knots (known as h)
         xi+first = array with spline values
         calc_pot.last[col1] - first + 1 = num of spline pts
         *(xi + first - 2) = value of endpoint gradient (default: 1e30)
         *(xi + first - 1) = value of other endpoint gradient
         (default: phi=0.0, rho=0.0, F=1e30)
         calc_pot.d2tab + first = array to hold 2nd deriv */
      spline_ed(calc_pot.step[col], xi + first, calc_pot.last[col] - first + 1,
	*(xi + first - 2), *(xi + first - 1), calc_pot.d2tab + first);
    }

#ifndef MPI
    myconf = nconf;
#endif /* MPI */

    /* region containing loop over configurations */
    {
      /* Temp variables */
      atom_t *atom;		// atom type
      int   h, j, k;
      int   n_i, n_j, n_k;
      int   uf;
#ifdef APOT
      double temp_eng;
#endif /* APOT */
#ifdef STRESS
      int   us, stresses;
#endif /* STRESS */

      /* Some useful temp struct variable types */
      /* neighbor pointers */
      neigh_t *neigh_j, *neigh_k;

      /* Pair variables */
      double phi_val, phi_grad;
      vector tmp_force;

      /* EAM variables */
      int   col_F;
      double eam_force;
#if defined NORESCALE && !defined APOT
      double rho_val;
#endif /* NORESCALE && !APOT */

      /* MEAM variables */
      double dV3j, dV3k, V3, vlj, vlk, vv3j, vv3k;
      vector dfj, dfk;
      angl *n_angl;

      /* Loop over configurations */
      for (h = firstconf; h < firstconf + myconf; h++) {
	uf = conf_uf[h - firstconf];
#ifdef STRESS
	us = conf_us[h - firstconf];
#endif /* STRESS */
	/* Reset energies */
	forces[energy_p + h] = 0.0;
#ifdef STRESS
	/* Reset stresses */
	stresses = stress_p + 6 * h;
	for (i = 0; i < 6; ++i)
	  forces[stresses + i] = 0.0;
#endif /* STRESS */

	/* Set limiting constraints */
	forces[limit_p + h] = -force_0[limit_p + h];

	/* FIRST LOOP: Reset forces and densities for each atom */
	for (i = 0; i < inconf[h]; i++) {
	  /* Skip every 3 spots in force array starting from position of first atom */
	  n_i = 3 * (cnfstart[h] + i);
	  if (uf) {
	    /* Set initial forces to negative of user given forces so we can take difference */
	    forces[n_i + 0] = -force_0[n_i + 0];
	    forces[n_i + 1] = -force_0[n_i + 1];
	    forces[n_i + 2] = -force_0[n_i + 2];
	  } else {
	    /* Set initial forces to zero if not using forces */
	    forces[n_i + 0] = 0.0;
	    forces[n_i + 1] = 0.0;
	    forces[n_i + 2] = 0.0;
	  }			/* uf */
	  /* Reset the density for each atom */
	  conf_atoms[cnfstart[h] - firstatom + i].rho = 0.0;
	}			/* i */
	/* END OF FIRST LOOP */

	/* SECOND LOOP: Calculate pair forces and energies, atomic densities */
	for (i = 0; i < inconf[h]; i++) {
	  /* Set pointer to temp atom pointer */
	  atom = conf_atoms + (cnfstart[h] - firstatom + i);
	  /* Skip every 3 spots for force array */
	  n_i = 3 * (cnfstart[h] + i);
	  /* Loop over neighbors */
	  for (j = 0; j < atom->num_neigh; j++) {
	    /* Set pointer to temp neighbor pointer */
	    neigh_j = atom->neigh + j;
	    /* Find the correct column in the potential table for pair potential: phi_ij
	       For Binary Alloy: 0 = phi_AA, 1 = (phi_AB or phi_BA), 2 = phi_BB
	       where typ = A = 0 and typ = B = 1 */
	    /* We need to check that neighbor atom exists inside pair potential's radius */
	    if (neigh_j->r < calc_pot.end[neigh_j->col[0]]) {
	      /* Compute phi and phi' value given radial distance
	         NOTE: slot = spline point index right below radial distance
	         shift = % distance from 'slot' spline pt
	         step = width of spline points (given as 'h' in books)
	         0 means the pair potential columns */
	      /* fn value and grad are calculated in the same step */
	      if (uf)
		phi_val =
		  splint_comb_dir(&calc_pot, xi, neigh_j->slot[0], neigh_j->shift[0], neigh_j->step[0],
		  &phi_grad);
	      else
		phi_val = splint_dir(&calc_pot, xi, neigh_j->slot[0], neigh_j->shift[0], neigh_j->step[0]);

	      /* Add in piece contributed by neighbor to energy */
	      forces[energy_p + h] += 0.5 * phi_val;

	      if (uf) {
		/* Compute tmp force values */
		tmp_force.x = neigh_j->dist_r.x * phi_grad;
		tmp_force.y = neigh_j->dist_r.y * phi_grad;
		tmp_force.z = neigh_j->dist_r.z * phi_grad;
		/* Add in force on atom i from atom j */
		forces[n_i + 0] += tmp_force.x;
		forces[n_i + 1] += tmp_force.y;
		forces[n_i + 2] += tmp_force.z;
#ifdef STRESS
		if (us) {
		  /* also calculate pair stresses */
		  forces[stresses + 0] -= 0.5 * neigh_j->dist.x * tmp_force.x;
		  forces[stresses + 1] -= 0.5 * neigh_j->dist.y * tmp_force.y;
		  forces[stresses + 2] -= 0.5 * neigh_j->dist.z * tmp_force.z;
		  forces[stresses + 3] -= 0.5 * neigh_j->dist.x * tmp_force.y;
		  forces[stresses + 4] -= 0.5 * neigh_j->dist.y * tmp_force.z;
		  forces[stresses + 5] -= 0.5 * neigh_j->dist.z * tmp_force.x;
		}
#endif /* STRESS */
	      }
	    }

	    /* r < cutoff */
	    /* END IF STMNT: NEIGH LIES INSIDE CUTOFF FOR PAIR POTENTIAL */
	    /* Find the correct column in the potential table for atomic density, rho_ij
	       paircol = number of pair potential columns
	       Binary Alloy: paircol = 3 (3 pair potentials with index 0, 1, 2)
	       index of densitiy functions: 3 = rho_A, 4 = rho_B
	       where A, B are atom type for the neighbor */
	    /* Compute rho rho value and sum them up
	       Need to play tricks so that rho values are put in the correct
	       columns if alloy. If atom j is A or B, fn value needs to be
	       in correct rho_A or rho_B respectively, it doesn't depend on atom i. */
	    /* Check that atom j lies inside rho_typ2 */
	    if (neigh_j->r < calc_pot.end[neigh_j->col[1]]) {
	      /* Store gradient in the neighbor for the pair r_ij
	         to be used in the future when computing forces
	         and sum up rho for atom i */
	      atom->rho +=
		splint_comb_dir(&calc_pot, xi, neigh_j->slot[1], neigh_j->shift[1], neigh_j->step[1],
		&neigh_j->drho);
	    } else {
	      /* If the pair distance does not lie inside rho_typ2
	         We set the grad to 0 so it doesn't sum into the net force */
	      neigh_j->drho = 0.0;
	    }			/* r < cutoff */

	    /* Compute the f_ij values and store the fn and grad in each neighbor struct for easy access later */

	    /* Find the correct column in the potential table for "f": f_ij
	       For Binary Alloy: 0 = f_AA, 1 = f_AB, f_BA, 2 = f_BB
	       where typ = A = 0 and typ = B = 1
	       Note: it is "paircol+2*ntypes" spots away in the array */

	    /* Check that atom j lies inside f_col2 */
	    if (neigh_j->r < calc_pot.end[neigh_j->col[2]]) {
	      /* Store the f(r_ij) value and the gradient for future use */
	      neigh_j->f =
		splint_comb_dir(&calc_pot, xi, neigh_j->slot[2], neigh_j->shift[2], neigh_j->step[2],
		&neigh_j->df);
	    } else {
	      /* Store f and f' = 0 if doesn't lie in boundary to be used later when calculating forces */
	      neigh_j->f = 0.0;
	      neigh_j->df = 0.0;
	    }

	    /* END LOOP OVER NEIGHBORS */
	  }

	  /* Find the correct column in the potential table for angle part: g_ijk
	     Binary Alloy: 0 = g_A, 1 = g_B
	     where A, B are atom type for the main atom i
	     Note: it is now "2*paircol+2*ntypes" from beginning column
	     to account for phi(paircol)+rho(nytpes)+F(ntypes)+f(paircol)
	     col2 = 2 * paircol + 2 * ntypes + typ1; */

	  /* Loop over every angle formed by neighbors
	     N(N-1)/2 possible combinations
	     Used in computing angular part g_ijk */

	  /* set n_angl pointer to angl_part of current atom */
	  n_angl = atom->angl_part;

	  for (j = 0; j < atom->num_neigh - 1; j++) {

	    /* Get pointer to neighbor jj */
	    neigh_j = atom->neigh + j;

	    for (k = j + 1; k < atom->num_neigh; k++) {

	      /* Get pointer to neighbor kk */
	      neigh_k = atom->neigh + k;

	      /* The cos(theta) should always lie inside -1 ... 1
	         So store the g and g' without checking bounds */
	      n_angl->g =
		splint_comb_dir(&calc_pot, xi, n_angl->slot, n_angl->shift, n_angl->step, &n_angl->dg);

	      /* Sum up rho piece for atom i caused by j and k
	         f_ij * f_ik * m_ijk */
	      atom->rho += neigh_j->f * neigh_k->f * n_angl->g;

	      /* Increase n_angl pointer */
	      n_angl++;
	    }
	  }

	  /* Column for embedding function, F */
	  col_F = paircol + ntypes + atom->type;

#ifndef NORESCALE
	  /* Compute energy, gradient for embedding function F
	     Check if rho lies short of inner cutoff of F(rho) */
	  if (atom->rho < calc_pot.begin[col_F]) {

	    /* Punish this potential for having rho lie outside of F */
	    forces[limit_p + h] += DUMMY_WEIGHT * 10. * dsquare(calc_pot.begin[col_F] - atom->rho);

	    /* Set the atomic density to the first rho in the spline F */
	    atom->rho = calc_pot.begin[col_F];

	  } else if (atom->rho > calc_pot.end[col_F]) {	/* rho is to the right of the spline */

	    /* Punish this potential for having rho lie outside of F */
	    forces[limit_p + h] += DUMMY_WEIGHT * 10. * dsquare(atom->rho - calc_pot.end[col_F]);

	    /* Set the atomic density to the last rho in the spline F */
	    atom->rho = calc_pot.end[col_F];
	  }
	  /* Compute energy piece from F, and store the gradient for later use */
	  forces[energy_p + h] += splint_comb(&calc_pot, xi, col_F, atom->rho, &atom->gradF);

#else
	  /* Compute energy, gradient for embedding function F
	     Check if rho lies short of inner cutoff of F(rho) */
	  if (atom->rho < calc_pot.begin[col_F]) {
#ifdef APOT
	    /* calculate analytic value explicitly */
	    apot_table.fvalue[col_F] (atom->rho, xi_opt + opt_pot.first[col_F], &temp_eng);
	    atom->gradF = apot_grad(atom->rho, xi_opt + opt_pot.first[col_F], apot_table.fvalue[col_F]);
	    forces[energy_p + h] += temp_eng;
#else
	    /* Linear extrapolate values to left to get F_i(rho)
	       This gets value and grad of initial spline point */
	    rho_val = splint_comb(&calc_pot, xi, col_F, calc_pot.begin[col_F], &atom->gradF);

	    /* Sum this to the total energy for this configuration
	       Linear extrapolate this energy */
	    forces[energy_p + h] += rho_val + (atom->rho - calc_pot.begin[col_F]) * atom->gradF;
#endif /* APOT */
	    /* rho is to the right of the spline */
	  } else if (atom->rho > calc_pot.end[col_F]) {
#ifdef APOT
	    /* calculate analytic value explicitly */
	    apot_table.fvalue[col_F] (atom->rho, xi_opt + opt_pot.first[col_F], &temp_eng);
	    atom->gradF = apot_grad(atom->rho, xi_opt + opt_pot.first[col_F], apot_table.fvalue[col_F]);
	    forces[energy_p + h] += temp_eng;
#else
	    /* Get value and grad at 1/2 the width from the final spline point */
	    rho_val =
	      splint_comb(&calc_pot, xi, col_F,
	      calc_pot.end[col_F] - .5 * calc_pot.step[col_F], &atom->gradF);
	    /* Linear extrapolate to the right to get energy */
	    forces[energy_p + h] += rho_val + (atom->rho - calc_pot.end[col_F]) * atom->gradF;
#endif /* APOT */
	    /* and in-between */
	  } else {
#ifdef APOT
	    /* calculate small values directly */
	    if (atom->rho < 0.1) {
	      apot_table.fvalue[col_F] (atom->rho, xi_opt + opt_pot.first[col_F], &temp_eng);
	      atom->gradF = apot_grad(atom->rho, xi_opt + opt_pot.first[col_F], apot_table.fvalue[col_F]);
	      forces[energy_p + h] += temp_eng;
	    } else
#endif
	      /* Get energy value from within spline and store the grad */
	      forces[energy_p + h] += splint_comb(&calc_pot, xi, col_F, atom->rho, &atom->gradF);
	  }
#endif /* !NORESCALE */

	  /* Sum up rho for future MPI use */
	  rho_sum_loc += atom->rho;

	  /* Calculate remaining forces from embedding function */

	  if (uf) {
	    /* Loop over neighbors */
	    for (j = 0; j < atom->num_neigh; ++j) {

	      /* Set pointer to temp neighbor pointer and record type */
	      neigh_j = atom->neigh + j;

	      /* Check that radial distance between pair is within
	         cutoff distance of either possible rho_A or rho_B
	         for alloys, where A or B stands for atom i
	         WARNING: Double check this!!! May not need this
	         since drho will be 0 otherwise */
	      if (neigh_j->r < calc_pot.end[neigh_j->col[1]]) {

		/* Calculate eam force */
		eam_force = neigh_j->drho * atom->gradF;

		/* Multiply the eamforce with x/r to get real force */
		tmp_force.x = neigh_j->dist_r.x * eam_force;
		tmp_force.y = neigh_j->dist_r.y * eam_force;
		tmp_force.z = neigh_j->dist_r.z * eam_force;

		/* Sum up forces acting on atom i from atom j */
		forces[n_i + 0] += tmp_force.x;
		forces[n_i + 1] += tmp_force.y;
		forces[n_i + 2] += tmp_force.z;

		/* Subtract off forces acting on atom j from atom i */
		n_j = 3 * neigh_j->nr;
		forces[n_j + 0] -= tmp_force.x;
		forces[n_j + 1] -= tmp_force.y;
		forces[n_j + 2] -= tmp_force.z;

#ifdef STRESS
		if (us) {
		  forces[stresses + 0] -= neigh_j->dist.x * tmp_force.x;
		  forces[stresses + 1] -= neigh_j->dist.y * tmp_force.y;
		  forces[stresses + 2] -= neigh_j->dist.z * tmp_force.z;
		  forces[stresses + 3] -= neigh_j->dist.x * tmp_force.y;
		  forces[stresses + 4] -= neigh_j->dist.y * tmp_force.z;
		  forces[stresses + 5] -= neigh_j->dist.z * tmp_force.x;
		}
#endif /* STRESS */
	      }			/* END IF STMT: Inside reach of rho cutoff */
	    }			/* END LOOP OVER NEIGHBORS */

	    /* Compute MEAM Forces */
	    /********************************/

	    /* Loop over every angle formed by neighbors
	       N(N-1)/2 possible combinations
	       Used in computing angular part g_ijk */

	    /* set n_angl pointer to angl_part of current atom */
	    n_angl = atom->angl_part;

	    for (j = 0; j < atom->num_neigh - 1; j++) {

	      /* Get pointer to neighbor j */
	      neigh_j = atom->neigh + j;
	      /* Force location for atom j */
	      n_j = 3 * neigh_j->nr;

	      for (k = j + 1; k < atom->num_neigh; k++) {

		/* Get pointer to neighbor k */
		neigh_k = atom->neigh + k;

		/* Force location for atom k */
		n_k = 3 * neigh_k->nr;

		/* Some tmp variables to clean up force fn below */
		dV3j = n_angl->g * neigh_j->df * neigh_k->f;
		dV3k = n_angl->g * neigh_j->f * neigh_k->df;
		V3 = neigh_j->f * neigh_k->f * n_angl->dg;

		vlj = V3 * neigh_j->inv_r;
		vlk = V3 * neigh_k->inv_r;
		vv3j = dV3j - vlj * n_angl->cos;
		vv3k = dV3k - vlk * n_angl->cos;

		dfj.x = vv3j * neigh_j->dist_r.x + vlj * neigh_k->dist_r.x;
		dfj.y = vv3j * neigh_j->dist_r.y + vlj * neigh_k->dist_r.y;
		dfj.z = vv3j * neigh_j->dist_r.z + vlj * neigh_k->dist_r.z;

		dfk.x = vv3k * neigh_k->dist_r.x + vlk * neigh_j->dist_r.x;
		dfk.y = vv3k * neigh_k->dist_r.y + vlk * neigh_j->dist_r.y;
		dfk.z = vv3k * neigh_k->dist_r.z + vlk * neigh_j->dist_r.z;

		/* Force on atom i from j and k */
		forces[n_i + 0] += atom->gradF * (dfj.x + dfk.x);
		forces[n_i + 1] += atom->gradF * (dfj.y + dfk.y);
		forces[n_i + 2] += atom->gradF * (dfj.z + dfk.z);

		/* Reaction force on atom j from i and k */
		forces[n_j + 0] -= atom->gradF * dfj.x;
		forces[n_j + 1] -= atom->gradF * dfj.y;
		forces[n_j + 2] -= atom->gradF * dfj.z;

		/* Reaction force on atom k from i and j */
		forces[n_k + 0] -= atom->gradF * dfk.x;
		forces[n_k + 1] -= atom->gradF * dfk.y;
		forces[n_k + 2] -= atom->gradF * dfk.z;

#ifdef STRESS
		if (us) {
		  /* Force from j on atom i */
		  tmp_force.x = atom->gradF * dfj.x;
		  tmp_force.y = atom->gradF * dfj.y;
		  tmp_force.z = atom->gradF * dfj.z;
		  forces[stresses + 0] -= neigh_j->dist.x * tmp_force.x;
		  forces[stresses + 1] -= neigh_j->dist.y * tmp_force.y;
		  forces[stresses + 2] -= neigh_j->dist.z * tmp_force.z;
		  forces[stresses + 3] -= neigh_j->dist.x * tmp_force.y;
		  forces[stresses + 4] -= neigh_j->dist.y * tmp_force.z;
		  forces[stresses + 5] -= neigh_j->dist.z * tmp_force.x;

		  /* Force from k on atom i */
		  tmp_force.x = atom->gradF * dfk.x;
		  tmp_force.y = atom->gradF * dfk.y;
		  tmp_force.z = atom->gradF * dfk.z;
		  forces[stresses + 0] -= neigh_k->dist.x * tmp_force.x;
		  forces[stresses + 1] -= neigh_k->dist.y * tmp_force.y;
		  forces[stresses + 2] -= neigh_k->dist.z * tmp_force.z;
		  forces[stresses + 3] -= neigh_k->dist.x * tmp_force.y;
		  forces[stresses + 4] -= neigh_k->dist.y * tmp_force.z;
		  forces[stresses + 5] -= neigh_k->dist.z * tmp_force.x;
		}
#endif // STRESS
		/* Increase n_angl pointer */
		n_angl++;
	      }			/* End inner loop over angles (neighbor atom k) */
	    }			/* End outer loop over angles (neighbor atom j) */
	  }			/* uf */
	}			/* END OF SECOND LOOP OVER ATOM i */

	/* 3RD LOOP OVER ATOM i */
	/* Sum up the square of the forces for each atom
	   then multiply it by the weight for this config */
	for (i = 0; i < inconf[h]; i++) {
	  atom = conf_atoms + i + cnfstart[h] - firstatom;
	  n_i = 3 * (cnfstart[h] + i);
#ifdef FWEIGHT
	  /* Weigh by absolute value of force */
	  forces[n_i + 0] /= FORCE_EPS + atom->absforce;
	  forces[n_i + 1] /= FORCE_EPS + atom->absforce;
	  forces[n_i + 2] /= FORCE_EPS + atom->absforce;
#endif /* FWEIGHT */

#ifdef CONTRIB
	  if (atom->contrib)
#endif /* CONTRIB */
	    tmpsum += conf_weight[h] *
	      (dsquare(forces[n_i + 0]) + dsquare(forces[n_i + 1]) + dsquare(forces[n_i + 2]));
	}			/* END OF THIRD LOOP OVER ATOM i */

	/* Add in the energy per atom and its weight to the sum */
	/* First divide by num atoms */
	forces[energy_p + h] /= (double)inconf[h];

	/* Then subtract off the cohesive energy given to use by user */
	forces[energy_p + h] -= force_0[energy_p + h];

	/* Sum up square of this new energy term for each config
	   multiplied by its respective weight */
	tmpsum += conf_weight[h] * eweight * dsquare(forces[energy_p + h]);

#ifdef STRESS
	/* LOOP OVER STRESSES */
	for (i = 0; i < 6; ++i) {
	  /* Multiply weight to stresses and divide by volume */
	  forces[stresses + i] /= conf_vol[h - firstconf];
	  /* Subtract off user supplied stresses */
	  forces[stresses + i] -= force_0[stresses + i];
	  /* Sum in the square of each stress component with config weight */
	  tmpsum += conf_weight[h] * sweight * dsquare(forces[stresses + i]);
	}
#endif /* STRESS */

#ifndef NORESCALE
	/* Add in the square of the limiting constraints for each config */
	/* This is punishment from going out of bounds for F(rho)
	   if NORESCALE is not defined */
	forces[limit_p + h] *= conf_weight[h];
	tmpsum += dsquare(forces[limit_p + h]);
#endif /* !NORESCALE */
      }				/* END MAIN LOOP OVER CONFIGURATIONS */
    }

#ifdef MPI
    /* Reduce the rho_sum into root node */
    MPI_Reduce(&rho_sum_loc, &rho_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    rho_sum = rho_sum_loc;
#endif // MPI

#ifdef NORESCALE
    if (myid == 0) {
      /* Calculate the average rho_sum per atom
         NOTE: This gauge constraint exists for both EAM and MEAM */
      rho_sum /= (double)natoms;

      /* Another constraint for the gauge conditions
         this sets the avg rho per atom to 1
         Please read the other constraint on gauge conditions
         above. */
      forces[dummy_p + ntypes] = DUMMY_WEIGHT * (rho_sum - 1.);
      tmpsum += dsquare(forces[dummy_p + ntypes]);
    }
#endif /* NORESCALE */

#ifdef MPI
    /* Reduce the global sum from all the tmpsum's */
    sum = 0.0;
    MPI_Reduce(&tmpsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    /* gather forces, energies, stresses */
    if (myid == 0) {		/* root node already has data in place */
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
    /* Set tmpsum to sum - only matters when not running MPI */
    sum = tmpsum;
#endif /* MPI */

    /* Root process only */
    if (myid == 0) {
      /* Increment function calls */
      fcalls++;
      /* If total sum is NAN return large number instead */
      if (isnan(sum)) {
#ifdef DEBUG
	printf("\n--> Force is nan! <--\n\n");
#endif /* DEBUG */
	return 10e10;
      } else
	return sum;
    }
  }				/* END OF INFINITE LOOP */

  /* Kill off other procs */
  return -1.0;
}

#endif /* MEAM */
