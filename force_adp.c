/****************************************************************
*
* force_adp.c: Routine used for calculating adp forces/energies
*
*****************************************************************/
/*
*   Copyright 2010 Daniel Schopf
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
* $Date: 2010/04/20 12:31:21 $
*****************************************************************/

#ifdef ADP
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

real calc_forces_adp(real *xi_opt, real *forces, int flag)
{
  real  tmpsum = 0, sum = 0;
  int   first, col, i;
  real *xi = NULL;

  static real rho_sum_loc, rho_sum;
  rho_sum_loc = rho_sum = 0.;

#if defined DEBUG && defined FORCES
  real  store_punish;
#endif

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
    rho_sum_loc = 0.;

#if !defined APOT
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
#else
    /* if flag==2 then the potential parameters have changed -> sync */
    if (flag == 2)
      potsync();
    /* non root processes hang on, unless...  */
    if (flag == 1)
      break;			/* Exception: flag 1 means clean up */
#endif /* APOT */

#endif /* MPI */

    /* init second derivatives for splines */

    /* pair potentials */
    for (col = 0; col < paircol; col++) {
      first = calc_pot.first[col];
      if (format == 3 || format == 0)
	spline_ed(calc_pot.step[col], xi + first,
		  calc_pot.last[col] - first + 1, *(xi + first - 2), 0.0,
		  calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
		  calc_pot.last[col] - first + 1, *(xi + first - 2), 0.0,
		  calc_pot.d2tab + first);
    }

    /* rho */
    for (col = paircol; col < paircol + ntypes; col++) {
      first = calc_pot.first[col];
      if (format == 3 || format == 0)
	spline_ed(calc_pot.step[col], xi + first,
		  calc_pot.last[col] - first + 1, *(xi + first - 2), 0.0,
		  calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
		  calc_pot.last[col] - first + 1, *(xi + first - 2), 0.0,
		  calc_pot.d2tab + first);
    }

    /* F */
    for (col = paircol + ntypes; col < paircol + 2 * ntypes; col++) {
      first = calc_pot.first[col];
      if (format == 3 || format == 0)
	spline_ed(calc_pot.step[col], xi + first,
		  calc_pot.last[col] - first + 1, *(xi + first - 2),
		  *(xi + first - 1), calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
		  calc_pot.last[col] - first + 1, *(xi + first - 2),
		  *(xi + first - 1), calc_pot.d2tab + first);
    }

    /* u */
    for (col = paircol + 2 * ntypes; col < 2 * paircol + 2 * ntypes; col++) {
      first = calc_pot.first[col];
      if (format == 3 || format == 0)
	spline_ed(calc_pot.step[col], xi + first,
		  calc_pot.last[col] - first + 1, *(xi + first - 2),
		  *(xi + first - 1), calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
		  calc_pot.last[col] - first + 1, *(xi + first - 2),
		  *(xi + first - 1), calc_pot.d2tab + first);
    }

    /* w */
    for (col = 2 * paircol + 2 * ntypes; col < 3 * paircol + 2 * ntypes;
	 col++) {
      first = calc_pot.first[col];
      if (format == 3 || format == 0)
	spline_ed(calc_pot.step[col], xi + first,
		  calc_pot.last[col] - first + 1, *(xi + first - 2),
		  *(xi + first - 1), calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
		  calc_pot.last[col] - first + 1, *(xi + first - 2),
		  *(xi + first - 1), calc_pot.d2tab + first);
    }

#ifndef MPI
    myconf = nconf;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif

    /* region containing loop over configurations,
       also OMP-parallelized region */
    {

      int   j, k, l, h;
      int   uf, us;
      int   typ1, typ2;
      /* columns for potentials 0=phi, 1=rho, 2=F, 3=u, 4=w */
      int   col0, col1, col2, col3, col4;
      int   self, stresses;
      real  r;

      real  phi_val;
      real  phi_grad, rho_grad;
      real  eam_force;
      vector u_force;
      sym_tens w_force;
      real  nu;
      real  rho_grad_j;

      real  tmp, trace;
      real  eng_store;

      vector tmp_force;
      atom_t *atom;
      neigh_t *neigh;

#ifdef _OPENMP
#pragma omp for reduction(+:tmpsum,rho_sum_loc)
#endif
      /* loop over configurations */
      for (h = firstconf; h < firstconf + myconf; h++) {
	uf = conf_uf[h - firstconf];
	us = conf_us[h - firstconf];
	/* reset energies and stresses */
	forces[energy_p + h] = 0.;
	for (i = 0; i < 6; i++)
	  forces[stress_p + 6 * h + i] = 0.;

	/* set limiting constraints */
	forces[limit_p + h] = -force_0[limit_p + h];

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
	  /* reset atomic density and dipole and quadrupol distortions */
	  conf_atoms[cnfstart[h] - firstatom + i].rho = 0.0;
	  conf_atoms[cnfstart[h] - firstatom + i].mu.x = 0.0;
	  conf_atoms[cnfstart[h] - firstatom + i].mu.y = 0.0;
	  conf_atoms[cnfstart[h] - firstatom + i].mu.z = 0.0;
	  conf_atoms[cnfstart[h] - firstatom + i].lambda.xx = 0.0;
	  conf_atoms[cnfstart[h] - firstatom + i].lambda.yy = 0.0;
	  conf_atoms[cnfstart[h] - firstatom + i].lambda.zz = 0.0;
	  conf_atoms[cnfstart[h] - firstatom + i].lambda.xy = 0.0;
	  conf_atoms[cnfstart[h] - firstatom + i].lambda.yz = 0.0;
	  conf_atoms[cnfstart[h] - firstatom + i].lambda.zx = 0.0;
	}
	/* end first loop */

	/* 2nd loop: calculate pair forces and energies, atomic densities. */
	for (i = 0; i < inconf[h]; i++) {
#if defined DEBUG && defined FORCES
	  fprintf(stderr, "\nWorking on atom %d\n", i);
#endif
	  atom = conf_atoms + i + cnfstart[h] - firstatom;
	  typ1 = atom->typ;
	  k = 3 * (cnfstart[h] + i);
	  /* loop over neighbours */
	  for (j = 0; j < atom->n_neigh; j++) {
#if defined DEBUG && defined FORCES
	    fprintf(stderr, "Working on atom %d neighbour %d\n", i, j);
#endif
	    neigh = atom->neigh + j;
	    /* only use neigbours with higher numbers,
	       others are calculated by actio=reactio */
	    if (neigh->nr >= i + cnfstart[h]) {
	      /* In small cells, an atom might interact with itself */
	      self = (neigh->nr == i + cnfstart[h]) ? 1 : 0;
	      typ2 = neigh->typ;

	      /* pair potential part */
	      col0 = (typ1 <= typ2) ?
		typ1 * ntypes + typ2 - ((typ1 * (typ1 + 1)) / 2)
		: typ2 * ntypes + typ1 - ((typ2 * (typ2 + 1)) / 2);
	      if (neigh->r < calc_pot.end[col0]) {
		/* fn value and grad are calculated in the same step */
		if (uf)
		  phi_val =
		    splint_comb_dir(&calc_pot, xi, col0, neigh->slot[0],
				    neigh->shift[0], neigh->step[0],
				    &phi_grad);
		else
		  phi_val =
		    splint_dir(&calc_pot, xi, col0, neigh->slot[0],
			       neigh->shift[0], neigh->step[0]);

		/* avoid double counting if atom is interacting with a
		   copy of itself */
		if (self) {
		  phi_val *= 0.5;
		  phi_grad *= 0.5;
		}

		/* not real force: cohesive energy */
		forces[energy_p + h] += phi_val;
#if defined DEBUG && defined FORCES
		fprintf(stderr, "pair-energy=%f\n", phi_val);
#endif
		if (uf) {
		  tmp_force.x = neigh->dist.x * phi_grad;
		  tmp_force.y = neigh->dist.y * phi_grad;
		  tmp_force.z = neigh->dist.z * phi_grad;
		  forces[k] += tmp_force.x;
		  forces[k + 1] += tmp_force.y;
		  forces[k + 2] += tmp_force.z;
#if defined DEBUG && defined FORCES
		  fprintf(stderr, "forces %f %f %f\n", k, forces[k],
			  forces[k + 1], forces[k + 2]);
#endif
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
	      /* dipole distortion part */
	      col3 = col0 + paircol + 2 * ntypes;
	      if (neigh->r < calc_pot.end[col3]) {
		if (uf)
		  neigh->u_val =
		    splint_comb_dir(&calc_pot, xi, col3, neigh->slot[2],
				    neigh->shift[2], neigh->step[2],
				    &neigh->u_grad);
		else
		  neigh->u_val =
		    splint_dir(&calc_pot, xi, col3, neigh->slot[2],
			       neigh->shift[2], neigh->step[2]);

		if (self)
		  neigh->u_val *= 0.5;
		tmp = neigh->u_val * neigh->rdist.x;
		atom->mu.x += tmp;
		conf_atoms[neigh->nr - firstatom].mu.x -= tmp;
		tmp = neigh->u_val * neigh->rdist.y;
		atom->mu.y += tmp;
		conf_atoms[neigh->nr - firstatom].mu.y -= tmp;
		tmp = neigh->u_val * neigh->rdist.z;
		atom->mu.z += tmp;
		conf_atoms[neigh->nr - firstatom].mu.z -= tmp;
	      }
	      /* quadrupole distortion part */
	      col4 = col0 + paircol + 3 * ntypes;
	      if (neigh->r < calc_pot.end[col3]) {
		if (uf)
		  neigh->w_val =
		    splint_comb_dir(&calc_pot, xi, col4, neigh->slot[3],
				    neigh->shift[3], neigh->step[3],
				    &neigh->w_grad);
		else
		  neigh->w_val =
		    splint_dir(&calc_pot, xi, col4, neigh->slot[3],
			       neigh->shift[3], neigh->step[3]);
		if (self)
		  neigh->w_val *= 0.5;
		/* diagonal elements */
		tmp = neigh->w_val * neigh->sqrdist.xx;
		atom->lambda.xx += tmp;
		conf_atoms[neigh->nr - firstatom].lambda.xx += tmp;
		tmp = neigh->w_val * neigh->sqrdist.yy;
		atom->lambda.yy += tmp;
		conf_atoms[neigh->nr - firstatom].lambda.yy += tmp;
		tmp = neigh->w_val * neigh->sqrdist.zz;
		atom->lambda.zz += tmp;
		conf_atoms[neigh->nr - firstatom].lambda.zz += tmp;
		/* offdiagonal elements */
		tmp = neigh->w_val * neigh->sqrdist.yz;
		atom->lambda.yz += tmp;
		conf_atoms[neigh->nr - firstatom].lambda.yz += tmp;
		tmp = neigh->w_val * neigh->sqrdist.zx;
		atom->lambda.zx += tmp;
		conf_atoms[neigh->nr - firstatom].lambda.zx += tmp;
		tmp = neigh->w_val * neigh->sqrdist.xy;
		atom->lambda.xy += tmp;
		conf_atoms[neigh->nr - firstatom].lambda.xy += tmp;
	      }

	      /* calculate atomic densities */
	      col1 = paircol + typ2;
	      if (typ2 == typ1) {
/* then transfer(a->b)==transfer(b->a) */
		if (neigh->r < calc_pot.end[col1]) {
		  phi_val = splint_dir(&calc_pot, xi, col1,
				       neigh->slot[1],
				       neigh->shift[1], neigh->step[1]);
		  atom->rho += phi_val;
#if defined DEBUG && defined FORCES
		  fprintf(stderr, "rho=%f (added %f) dist=%f\n", atom->rho,
			  phi_val, neigh->r);
#endif
		  /* avoid double counting if atom is interacting with a
		     copy of itself */
		  if (!self) {
		    conf_atoms[neigh->nr - firstatom].rho += phi_val;
		  }
		}
	      } else {		/* transfer(a->b)!=transfer(b->a) */
		col0 = paircol + typ1;
		if (neigh->r < calc_pot.end[col1]) {
		  atom->rho += splint_dir(&calc_pot, xi, col1,
					  neigh->slot[1],
					  neigh->shift[1], neigh->step[1]);
#if defined DEBUG && defined FORCES
		  fprintf(stderr, "rho=%f (added %f)\n", atom->rho,
			  splint_dir(&calc_pot, xi, col1, neigh->slot[1],
				     neigh->shift[1], neigh->step[1]));
#endif
		}
		/* cannot use slot/shift to access splines */
		if (neigh->r < calc_pot.end[col0])
		  conf_atoms[neigh->nr - firstatom].rho +=
		    splint(&calc_pot, xi, col0, neigh->r);
	      }
	    }			/*  neighbours with bigger atom nr */
	  }			/* loop over neighbours */

	  col2 = paircol + ntypes + typ1;	/* column of F */
#ifndef NORESCALE
	  if (atom->rho > calc_pot.end[col2]) {
	    /* then punish target function -> bad potential */
	    forces[limit_p + h] += DUMMY_WEIGHT *
	      10. * SQR(atom->rho - calc_pot.end[col2]);
#ifndef PARABEL
/* then we use the final value, with PARABEL: extrapolate */
	    atom->rho = calc_pot.end[col2];
#endif /* PARABEL */
	  }

	  if (atom->rho < calc_pot.begin[col2]) {
	    /* then punish target function -> bad potential */
	    forces[limit_p + h] += DUMMY_WEIGHT *
	      10. * SQR(calc_pot.begin[col2] - atom->rho);
#ifndef PARABEL
/* then we use the final value, with PARABEL: extrapolate */
	    atom->rho = calc_pot.begin[col2];
#endif /* PARABEL */
	  }
#endif /* NOT NORESCALE */
	  /* embedding energy, embedding gradient */
	  /* contribution to cohesive energy is F(n) */
#ifdef PARABEL
	  forces[energy_p + h] +=
	    parab_comb(&calc_pot, xi, col2, atom->rho, &atom->gradF);
#elif defined(NORESCALE)
	  if (atom->rho < calc_pot.begin[col2]) {
	    /* linear extrapolation left */
	    phi_val =
	      splint_comb(&calc_pot, xi, col2, calc_pot.begin[col2],
			  &atom->gradF);
	    forces[energy_p + h] +=
	      phi_val + (atom->rho - calc_pot.begin[col2]) * atom->gradF;
	  } else if (atom->rho > calc_pot.end[col2]) {
	    /* and right */
	    phi_val =
	      splint_comb(&calc_pot, xi, col2,
			  calc_pot.end[col2] - .5 * calc_pot.step[col2],
			  &atom->gradF);
	    forces[energy_p + h] +=
	      phi_val + (atom->rho - calc_pot.end[col2]) * atom->gradF;
	  }
	  /* and in-between */
	  else {
	    forces[energy_p + h] +=
	      splint_comb(&calc_pot, xi, col2, atom->rho, &atom->gradF);
	  }
#else
	  forces[energy_p + h] +=
	    splint_comb(&calc_pot, xi, col2, atom->rho, &atom->gradF);
#endif
#if defined DEBUG && defined FORCES
	  fprintf(stderr, "total eam energy: %f (rho=%f)\n",
		  forces[energy_p + h], atom->rho);
#endif
	  /* sum up rho */
	  rho_sum_loc += atom->rho;

	  eng_store = 0;
	  /* calculate ADP parameters for atom i */
	  eng_store += SQR(atom->mu.x);
	  eng_store += SQR(atom->mu.y);
	  eng_store += SQR(atom->mu.z);
	  atom->nu = atom->lambda.xx + atom->lambda.yy + atom->lambda.zz;
	  trace = atom->nu / 3.0;
	  eng_store += SQR(atom->lambda.xx - trace);
	  eng_store += SQR(atom->lambda.yy - trace);
	  eng_store += SQR(atom->lambda.zz - trace);
	  eng_store += SQR(atom->lambda.yz) * 2.0;
	  eng_store += SQR(atom->lambda.zx) * 2.0;
	  eng_store += SQR(atom->lambda.xy) * 2.0;
	  eng_store *= 0.5;
	  forces[energy_p + h] += eng_store;
	}			/* second loop over atoms */

	/* 3rd loop over atom: ADP force */
	if (uf) {		/* only required if we calc forces */
	  for (i = 0; i < inconf[h]; i++) {
	    atom = conf_atoms + i + cnfstart[h] - firstatom;
	    typ1 = atom->typ;
	    k = 3 * (cnfstart[h] + i);
	    col2 = paircol + ntypes + typ1;	/* column of F */
	    for (j = 0; j < atom->n_neigh; j++) {
	      /* loop over neighbours */
	      neigh = atom->neigh + j;
	      /* only neigbours higher than current atom are of interest */
	      if (neigh->nr >= i + cnfstart[h]) {
		/* In small cells, an atom might interact with itself */
		self = (neigh->nr == i + cnfstart[h]) ? 1 : 0;
		typ2 = neigh->typ;
		col1 = paircol + typ2;
		r = neigh->r;
		/* are we within reach? */
		if ((r < calc_pot.end[col1])
		    || (r < calc_pot.end[col2 - ntypes])) {
		  rho_grad =
		    (r < calc_pot.end[col1]) ?
		    splint_grad_dir(&calc_pot, xi, col1,
				    neigh->slot[1], neigh->shift[1],
				    neigh->step[1]) : 0.;
		  if (typ2 == typ1)	/* use actio = reactio */
		    rho_grad_j = rho_grad;
		  else
		    rho_grad_j = (r < calc_pot.end[col2 - ntypes]) ?
		      splint_grad(&calc_pot, xi, col2 - ntypes, r) : 0.;
		  /* now we know everything - calculate forces */
		  eam_force = (rho_grad * atom->gradF +
			       rho_grad_j *
			       conf_atoms[(neigh->nr) - firstatom].gradF);
#if defined DEBUG && defined FORCES
		  fprintf(stderr,
			  "eam_force %f grad %f gradF %f grad2 %f gradF %f\n",
			  eam_force, phi_grad, atom->gradF, grad2,
			  conf_atoms[(neigh->nr) - firstatom].gradF);
#endif
		  /* avoid double counting if atom is interacting with a
		     copy of itself */
		  if (self)
		    eam_force *= 0.5;
		  tmp_force.x = neigh->dist.x * eam_force;
		  tmp_force.y = neigh->dist.y * eam_force;
		  tmp_force.z = neigh->dist.z * eam_force;
		  forces[k] += tmp_force.x;
		  forces[k + 1] += tmp_force.y;
		  forces[k + 2] += tmp_force.z;
#if defined DEBUG && defined FORCES
		  fprintf(stderr, "EAM k=%d dist %f %f %f eam_force %f\n", k,
			  neigh->dist.x, neigh->dist.y, neigh->dist.z,
			  eam_force);
		  fprintf(stderr, "EAM k=%d forces %f %f %f\n", k, forces[k],
			  forces[k + 1], forces[k + 2]);
#endif
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
		    stresses = stress_p + 6 * h;
		    forces[stresses] -= neigh->dist.x * tmp_force.x;
		    forces[stresses + 1] -= neigh->dist.y * tmp_force.y;
		    forces[stresses + 2] -= neigh->dist.z * tmp_force.z;
		    forces[stresses + 3] -= neigh->dist.x * tmp_force.y;
		    forces[stresses + 4] -= neigh->dist.y * tmp_force.z;
		    forces[stresses + 5] -= neigh->dist.z * tmp_force.x;
		  }
#endif /* STRESS */
		}		/* within reach */
		col3 = col0 + paircol + 2 * ntypes;
		if (neigh->r < calc_pot.end[col3]) {
		  u_force.x =
		    (atom->mu.x - conf_atoms[(neigh->nr) - firstatom].mu.x);
		  u_force.y =
		    (atom->mu.y - conf_atoms[(neigh->nr) - firstatom].mu.y);
		  u_force.z =
		    (atom->mu.z - conf_atoms[(neigh->nr) - firstatom].mu.z);
		  /* avoid double counting if atom is interacting with a
		     copy of itself */
		  if (self) {
		    u_force.x *= 0.5;
		    u_force.y *= 0.5;
		    u_force.z *= 0.5;
		  }
		  tmp = SPROD(u_force, neigh->dist) * neigh->u_grad;
		  tmp_force.x =
		    u_force.x * neigh->u_val + tmp * neigh->rdist.x;
		  tmp_force.y =
		    u_force.y * neigh->u_val + tmp * neigh->rdist.y;
		  tmp_force.z =
		    u_force.z * neigh->u_val + tmp * neigh->rdist.z;
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
		col4 = col0 + paircol + 3 * ntypes;
		if (neigh->r < calc_pot.end[col3]) {
		  w_force.xx =
		    (atom->lambda.xx +
		     conf_atoms[(neigh->nr) - firstatom].lambda.xx);
		  w_force.yy =
		    (atom->lambda.yy +
		     conf_atoms[(neigh->nr) - firstatom].lambda.yy);
		  w_force.zz =
		    (atom->lambda.zz +
		     conf_atoms[(neigh->nr) - firstatom].lambda.zz);
		  w_force.yz =
		    (atom->lambda.yz +
		     conf_atoms[(neigh->nr) - firstatom].lambda.yz);
		  w_force.zx =
		    (atom->lambda.zx +
		     conf_atoms[(neigh->nr) - firstatom].lambda.zx);
		  w_force.xy =
		    (atom->lambda.xy +
		     conf_atoms[(neigh->nr) - firstatom].lambda.xy);
		  nu = (atom->nu + conf_atoms[(neigh->nr)-firstatom].nu) / 3.0;
		  /* avoid double counting if atom is interacting with a
		     copy of itself */
		  if (self) {
		    w_force.xx *= 0.5;
		    w_force.yy *= 0.5;
		    w_force.zz *= 0.5;
		    w_force.yz *= 0.5;
		    w_force.zx *= 0.5;
		    w_force.xy *= 0.5;
		  }
		  tmp = SPROD(u_force, neigh->dist) * neigh->u_grad;
		  tmp_force.x =
		    u_force.x * neigh->u_val + tmp * neigh->rdist.x;
		  tmp_force.y =
		    u_force.y * neigh->u_val + tmp * neigh->rdist.y;
		  tmp_force.z =
		    u_force.z * neigh->u_val + tmp * neigh->rdist.z;
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
	      conf_weight[h] * (SQR(forces[k]) + SQR(forces[k + 1]) +
				SQR(forces[k + 2]));
#if defined DEBUG && defined FORCES
	    fprintf(stderr, "k=%d forces %f %f %f tmpsum=%f\n", k, forces[k],
		    forces[k + 1], forces[k + 2], tmpsum);
#endif
	  }			/* third loop over atoms */
	}

	/* use forces */
	/* energy contributions */
	forces[energy_p + h] *= eweight / (real)inconf[h];
	forces[energy_p + h] -= force_0[energy_p + h];
	tmpsum += conf_weight[h] * SQR(forces[energy_p + h]);
#if defined DEBUG && defined FORCES
	fprintf(stderr, "energy: tmpsum=%f energy=%f\n", tmpsum,
		forces[energy_p + h]);
#endif
#ifdef STRESS
	/* stress contributions */
	if (uf && us) {
	  for (i = 0; i < 6; i++) {
	    forces[stress_p + 6 * h + i] *= sweight / conf_vol[h - firstconf];
	    forces[stress_p + 6 * h + i] -= force_0[stress_p + 6 * h + i];
	    tmpsum += SQR(conf_weight[h] * forces[stress_p + 6 * h + i]);
	  }
	}
#endif /* STRESS */
	/* limiting constraints per configuration */
	tmpsum += conf_weight[h] * SQR(forces[limit_p + h]);

      }				/* loop over configurations */
    }				/* parallel region */
#ifdef MPI
    /* Reduce rho_sum */
    MPI_Reduce(&rho_sum_loc, &rho_sum, 1, REAL, MPI_SUM, 0, MPI_COMM_WORLD);
#else /* MPI */
    rho_sum = rho_sum_loc;
#endif /* MPI */

    /* dummy constraints (global) */
#ifdef APOT
    /* add punishment for out of bounds (mostly for powell_lsq) */
    if (myid == 0) {
      tmpsum += apot_punish(xi_opt, forces);
#if defined DEBUG && defined FORCES
      fprintf(stderr, "\napot punishments: tmpsum=%f punish=%f\n", tmpsum,
	      apot_punish(xi_opt, forces));
      store_punish = tmpsum;
#endif
    }
#endif

    if (myid == 0) {
      int   g;
      for (g = 0; g < ntypes; g++) {
	/* PARABEL, WZERO, NORESC - different behaviour */
#ifdef PARABEL
/* constraints on U(n) */
	forces[dummy_p + ntypes + g] = DUMMY_WEIGHT *
	  parab(&calc_pot, xi, paircol + ntypes + g, 0.)
	  - force_0[dummy_p + ntypes + g];
/* constraints on U`(n) */
	forces[dummy_p + g] = DUMMY_WEIGHT *
	  parab_grad(&calc_pot, xi, paircol + ntypes + g, .5 *
		     (calc_pot.begin[paircol + ntypes + g] +
		      calc_pot.end[paircol + ntypes + g])) -
	  force_0[dummy_p + g];
#elif defined(WZERO)
	if (calc_pot.begin[paircol + ntypes + g] <= 0.)
	  /* 0 in domain of U(n) */
/* constraints on U(n) */
	  forces[dummy_p + ntypes + g] = DUMMY_WEIGHT *
	    splint(&calc_pot, xi, paircol + ntypes + g, 0.)
	    - force_0[dummy_p + ntypes + g];
	else
	  /* 0 not in domain of U(n) */
	  forces[dummy_p + ntypes + g] = 0.;	/* Free end... */
/* constraints on U`(n) */
	forces[dummy_p + g] = DUMMY_WEIGHT *
	  splint_grad(&calc_pot, xi, paircol + ntypes + g, .5 *
		      (calc_pot.begin[paircol + ntypes + g] +
		       calc_pot.end[paircol + ntypes + g]))
	  - force_0[dummy_p + g];
#elif defined(NORESCALE)
	/* clear field */
	forces[dummy_p + ntypes + g] = 0.;	/* Free end... */
	/* NEW: Constraint on U': U'(1.)=0; */
	forces[dummy_p + g] = DUMMY_WEIGHT *
	  splint_grad(&calc_pot, xi, paircol + ntypes + g, 1.);
#else /* NOTHING */
	forces[dummy_p + ntypes + g] = 0.;	/* Free end... */
/* constraints on U`(n) */
	forces[dummy_p + g] = DUMMY_WEIGHT *
	  splint_grad(&calc_pot, xi, paircol + ntypes + g, .5 *
		      (calc_pot.begin[paircol + ntypes + g] +
		       calc_pot.end[paircol + ntypes + g]))
	  - force_0[dummy_p + g];
#endif /* Dummy constraints */
	tmpsum += SQR(forces[dummy_p + ntypes + g]);
	tmpsum += SQR(forces[dummy_p + g]);
#if defined DEBUG && defined FORCES
	fprintf(stderr, "dummy constraints on U: tmpsum=%f punish=%f\n",
		tmpsum, forces[dummy_p + ntypes + g]);
	fprintf(stderr, "dummy constraints on U': tmpsum=%f punish=%f\n",
		tmpsum, forces[dummy_p + g]);
#endif

      }				/* loop over types */
#ifdef NORESCALE
      /* NEW: Constraint on n: <n>=1. ONE CONSTRAINT ONLY */
      /* Calculate averages */
      rho_sum /= (real)natoms;
      /* ATTN: if there are invariant potentials, things might be problematic */
      forces[dummy_p + ntypes] = DUMMY_WEIGHT * (rho_sum - 1.);
      tmpsum += SQR(forces[dummy_p + ntypes]);
#endif
#if defined DEBUG && defined FORCES
/*      fprintf(stderr, "limiting constraints: tmpsum=%f punish=%f\n", tmpsum,*/
/*              forces[limit_p]);*/
      fprintf(stderr, "total EAM punishments: tmpsum=%f punish^2=%f\n\n",
	      tmpsum, tmpsum - store_punish);
#endif
    }				/* only root process */
    sum = tmpsum;		/* global sum = local sum  */
#ifdef MPI
    /* reduce global sum */
    sum = 0.;
    MPI_Reduce(&tmpsum, &sum, 1, REAL, MPI_SUM, 0, MPI_COMM_WORLD);
    /* gather forces, energies, stresses */
    /* forces */
    MPI_Gatherv(forces + firstatom * 3, myatoms, MPI_VEKTOR,
		forces, atom_len, atom_dist, MPI_VEKTOR, 0, MPI_COMM_WORLD);
    /* energies */
    MPI_Gatherv(forces + natoms * 3 + firstconf, myconf, REAL,
		forces
		+ natoms * 3, conf_len, conf_dist, REAL, 0, MPI_COMM_WORLD);
    /* stresses */
    MPI_Gatherv(forces + natoms * 3 + nconf +
		6 * firstconf, myconf, MPI_STENS,
		forces + natoms * 3 + nconf,
		conf_len, conf_dist, MPI_STENS, 0, MPI_COMM_WORLD);
    /* punishment constraints */
    MPI_Gatherv(forces + natoms * 3 + 7 * nconf +
		firstconf, myconf, REAL,
		forces + natoms * 3 + 7 * nconf,
		conf_len, conf_dist, REAL, 0, MPI_COMM_WORLD);
    /* no need to pick up dummy constraints - are already @ root */
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

#endif /* ADP */
