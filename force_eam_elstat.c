/****************************************************************
 *
 * force_eam_elstat.c: Routine used for combination of calculating 
 * pair/monopole/dipole forces/energies in various interpolation schemes.
 * AND eam forces/energies
 ****************************************************************
 *
 * Copyright 2010-2011
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.itap.physik.uni-stuttgart.de/
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

#if defined COULOMB && defined EAM

#include "potfit.h"

#include "functions.h"
#include "potential.h"
#include "splines.h"
#include "utils.h"

/****************************************************************
 *
 *  compute forces using pair potentials with spline interpolation
 *
 *  returns sum of squares of differences between calculated and reference
 *     values
 *
 *  arguments: *xi - pointer to short-range potential
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

double calc_forces_eam_elstat(double *xi_opt, double *forces, int flag)
{
  double  tmpsum, sum = 0.;
  int   first, col, ne, size, i;
  double *xi = NULL;
  apot_table_t *apt = &apot_table;
  double  charge[ntypes];
  double  sum_charges;
  double  dp_kappa;
  int ncols;

  ncols = ntypes * (ntypes + 1) / 2;

#ifdef DIPOLE
  double  dp_alpha[ntypes];
  double  dp_b[ncols];
  double  dp_c[ncols];
#endif /* DIPOLE */

  static double rho_sum_loc, rho_sum;
  rho_sum_loc = rho_sum = 0.;

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

  ne = apot_table.total_ne_par;
  size = apt->number;

  /* This is the start of an infinite loop */
  while (1) {
    tmpsum = 0.;		/* sum of squares of local process */
    rho_sum_loc = 0.;

#if defined APOT && !defined MPI
    if (format == 0) {
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
    if (format == 0)
      update_calc_table(xi_opt, xi, 0);
#else /* APOT */
    /* if flag==2 then the potential parameters have changed -> sync */
    if (flag == 2)
      potsync();
#endif /* APOT */
#endif /* MPI */

    /* local arrays for electrostatic parameters */
    sum_charges = 0;
    for (i = 0; i < ntypes - 1; i++) {
      if (xi_opt[2 * size + ne + i]) {
	charge[i] = xi_opt[2 * size + ne + i];
	sum_charges += apt->ratio[i] * charge[i];
      } else {
	charge[i] = 0.;
      }
    }
    apt->last_charge = -sum_charges / apt->ratio[ntypes - 1];
    charge[ntypes - 1] = apt->last_charge;
    if (xi_opt[2 * size + ne + ntypes - 1]) {
      dp_kappa = xi_opt[2 * size + ne + ntypes - 1];
    } else {
      dp_kappa = 0.;
    }

    /*int k;
    for(k = 0; k < 71; k++) {
      printf("%d\t%f\n", k, xi_opt[k]);
      }
      error(1, "");*/

#ifdef DIPOLE
    for (i = 0; i < ntypes; i++) {
      if (xi_opt[2 * size + ne + ntypes + i]) {
	dp_alpha[i] = xi_opt[2 * size + ne + ntypes + i];
      } else {
	dp_alpha[i] = 0.;
      }
    }
    for (i = 0; i < ncols; i++) {
      if (xi_opt[2 * size + ne + 2 * ntypes + i]) {
	dp_b[i] = xi_opt[2 * size + ne + 2 * ntypes + i];
      } else {
	dp_b[i] = 0.;
      }
      if (xi_opt[2 * size + ne + 2 * ntypes + ncols + i]) {
	dp_c[i] = xi_opt[2 * size + ne + 2 * ntypes + ncols + i];
      } else {
	dp_c[i] = 0.;
      }
    }
#endif /* DIPOLE */

    /* init second derivatives for splines */

    /* pair potentials */
    for (col = 0; col < paircol; col++) {
      first = calc_pot.first[col];
      if (format == 3 || format == 0) {
	spline_ed(calc_pot.step[col], xi + first,
	  calc_pot.last[col] - first + 1, *(xi + first - 2), 0.0,
	  calc_pot.d2tab + first);
      } else {			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
	  calc_pot.last[col] - first + 1, *(xi + first - 2), 0.0,
	  calc_pot.d2tab + first);
      }
    }

    /* rho */
    for (col = paircol; col < paircol + ntypes; col++) {
      first = calc_pot.first[col];
      if (format == 0 || format == 3)
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
      /* gradient at left boundary matched to square root function,
         when 0 not in domain(F), else natural spline */
      if (format == 0 || format == 3)
	spline_ed(calc_pot.step[col], xi + first,
	  calc_pot.last[col] - first + 1,
	  *(xi + first - 2), *(xi + first - 1),
	  calc_pot.d2tab + first);
      else			/* format >= 4 ! */
	spline_ne(calc_pot.xcoord + first, xi + first,
	  calc_pot.last[col] - first + 1,
	  *(xi + first - 2), *(xi + first - 1),
	  calc_pot.d2tab + first);
    }


#ifndef MPI
    myconf = nconf;
#endif /* MPI */

    /* region containing loop over configurations,
       also OMP-parallelized region */
    {
      int   self;
      vector tmp_force;
      int   h, j, k, l, typ1, typ2, uf, us, stresses;	/* config */
      double  fnval, grad, fnval_tail, grad_tail, grad_i, grad_j, p_sr_tail;
      atom_t *atom;
      neigh_t *neigh;
      double  r;
      int   col_F;
      double  eam_force;
      double  rho_val, rho_grad, rho_grad_j;

      /* loop over configurations: M A I N LOOP CONTAINING ALL ATOM-LOOPS */
      for (h = firstconf; h < firstconf + myconf; h++) {
	uf = conf_uf[h - firstconf];
	us = conf_us[h - firstconf];
	/* reset energies and stresses */
	forces[energy_p + h] = 0.;
	for (i = 0; i < 6; i++)
	  forces[stress_p + 6 * h + i] = 0.;

	/* set limiting constraints */
	forces[limit_p + h] = -force_0[limit_p + h];

#ifdef DIPOLE
	/* reset dipoles and fields: LOOP Z E R O */
	for (i = 0; i < inconf[h]; i++) {
	  atom = conf_atoms + i + cnfstart[h] - firstatom;
	  atom->E_stat.x = 0.;
	  atom->E_stat.y = 0.;
	  atom->E_stat.z = 0.;
	  atom->p_sr.x = 0.;
	  atom->p_sr.y = 0.;
	  atom->p_sr.z = 0.;
	}
#endif /* DIPOLE */

	/* F I R S T LOOP OVER ATOMS: reset forces, dipoles */
	for (i = 0; i < inconf[h]; i++) {	/* atoms */
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
	  /* reset atomic density */
	  conf_atoms[cnfstart[h] - firstatom + i].rho = 0.0;
	}			/* end F I R S T LOOP */

	/* S E C O N D loop: calculate short-range and monopole forces,
	   calculate static field- and dipole-contributions,
	   calculate atomic densities */
	for (i = 0; i < inconf[h]; i++) {	/* atoms */
	  atom = conf_atoms + i + cnfstart[h] - firstatom;
	  typ1 = atom->typ;
	  k = 3 * (cnfstart[h] + i);
	  for (j = 0; j < atom->n_neigh; j++) {	/* neighbors */
	    neigh = atom->neigh + j;
	    typ2 = neigh->typ;
	    col = neigh->col[0];

	    /* updating tail-functions - only necessary with variing kappa */
	    if (!apt->sw_kappa) {
	      elstat_shift(neigh->r, dp_kappa, &neigh->fnval_el,
		&neigh->grad_el, &neigh->ggrad_el);
	    }

	    /* In small cells, an atom might interact with itself */
	    self = (neigh->nr == i + cnfstart[h]) ? 1 : 0;

	    /* calculate short-range forces */
	    if (neigh->r < calc_pot.end[col]) {


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
	    if (neigh->r < dp_cut && (charge[typ1]
		|| charge[typ2])) {

	      fnval_tail = neigh->fnval_el;
	      grad_tail = neigh->grad_el;

	      grad_i = charge[typ2] * grad_tail;
	      if (typ1 == typ2) {
		grad_j = grad_i;
	      } else {
		grad_j = charge[typ1] * grad_tail;
	      }
	      fnval = charge[typ1] * charge[typ2] * fnval_tail;
	      grad = charge[typ1] * grad_i;

	      if (self) {
		grad_i *= 0.5;
		grad_j *= 0.5;
		fnval *= 0.5;
		grad *= 0.5;
	      }

	      forces[energy_p + h] += fnval;

	      if (uf) {
		tmp_force.x = neigh->dist.x * grad * neigh->r;
		tmp_force.y = neigh->dist.y * grad * neigh->r;
		tmp_force.z = neigh->dist.z * grad * neigh->r;
		forces[k] += tmp_force.x;
		forces[k + 1] += tmp_force.y;
		forces[k + 2] += tmp_force.z;
		/* actio = reactio */
		l = 3 * neigh->nr;
		forces[l] -= tmp_force.x;
		forces[l + 1] -= tmp_force.y;
		forces[l + 2] -= tmp_force.z;

#ifdef STRESS
		/* calculate coulomb stresses */
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
#ifdef DIPOLE
	      /* calculate static field-contributions */
	      atom->E_stat.x += neigh->dist.x * neigh->r * grad_i;
	      atom->E_stat.y += neigh->dist.y * neigh->r * grad_i;
	      atom->E_stat.z += neigh->dist.z * neigh->r * grad_i;

	      conf_atoms[neigh->nr - firstatom].E_stat.x -=
		neigh->dist.x * neigh->r * grad_j;
	      conf_atoms[neigh->nr - firstatom].E_stat.y -=
		neigh->dist.y * neigh->r * grad_j;
	      conf_atoms[neigh->nr - firstatom].E_stat.z -=
		neigh->dist.z * neigh->r * grad_j;

	      /* calculate short-range dipoles */
	      if (dp_alpha[typ1] && dp_b[col] && dp_c[col]) {
		p_sr_tail =
		  grad_tail * neigh->r * shortrange_value(neigh->r,
		  dp_alpha[typ1], dp_b[col], dp_c[col]);
		atom->p_sr.x += charge[typ2] * neigh->dist.x * p_sr_tail;
		atom->p_sr.y += charge[typ2] * neigh->dist.y * p_sr_tail;
		atom->p_sr.z += charge[typ2] * neigh->dist.z * p_sr_tail;
	      }
	      if (dp_alpha[typ2] && dp_b[col] && dp_c[col] && !self) {
		p_sr_tail =
		  grad_tail * neigh->r * shortrange_value(neigh->r,
		  dp_alpha[typ2], dp_b[col], dp_c[col]);
		conf_atoms[neigh->nr - firstatom].p_sr.x -=
		  charge[typ1] * neigh->dist.x * p_sr_tail;
		conf_atoms[neigh->nr - firstatom].p_sr.y -=
		  charge[typ1] * neigh->dist.y * p_sr_tail;
		conf_atoms[neigh->nr - firstatom].p_sr.z -=
		  charge[typ1] * neigh->dist.z * p_sr_tail;
	      }
#endif /* DIPOLE */

	    }

	    /* calculate atomic densities */
	    if (atom->typ == neigh->typ) {
	      /* then transfer(a->b)==transfer(b->a) */
	      if (neigh->r < calc_pot.end[neigh->col[1]]) {
		rho_val =
		  splint_dir(&calc_pot, xi, neigh->slot[1], neigh->shift[1],
		  neigh->step[1]);
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
		atom->rho +=
		  splint_dir(&calc_pot, xi, neigh->slot[1], neigh->shift[1],
		  neigh->step[1]);
	      }
	      /* cannot use slot/shift to access splines */
	      if (neigh->r < calc_pot.end[paircol + atom->typ])
		conf_atoms[neigh->nr - firstatom].rho +=
		  splint(&calc_pot, xi, paircol + atom->typ, neigh->r);
	    }

	  }			/* loop over neighbours */

	  col_F = paircol + ntypes + atom->typ;	/* column of F */
	  if (atom->rho > calc_pot.end[col_F]) {
	    /* then punish target function -> bad potential */
	    forces[limit_p + h] +=
	      DUMMY_WEIGHT * 10. * dsquare(atom->rho - calc_pot.end[col_F]);
	    atom->rho = calc_pot.end[col_F];
	  }

	  if (atom->rho < calc_pot.begin[col_F]) {
	    /* then punish target function -> bad potential */
	    forces[limit_p + h] +=
	      DUMMY_WEIGHT * 10. * dsquare(calc_pot.begin[col_F] - atom->rho);
	    atom->rho = calc_pot.begin[col_F];
	  }

	  /* embedding energy, embedding gradient */
	  /* contribution to cohesive energy is F(n) */
#ifdef NORESCALE
	  if (atom->rho < calc_pot.begin[col_F]) {
	    /* linear extrapolation left */
	    rho_val =
	      splint_comb(&calc_pot, xi, col_F, calc_pot.begin[col_F],
	      &atom->gradF);
	    forces[energy_p + h] +=
	      rho_val + (atom->rho - calc_pot.begin[col_F]) * atom->gradF;
#ifdef APOT
	    forces[limit_p + h] +=
	      DUMMY_WEIGHT * 10. * dsquare(calc_pot.begin[col_F] - atom->rho);
#endif /* APOT */
	  } else if (atom->rho > calc_pot.end[col_F]) {
	    /* and right */
	    rho_val =
	      splint_comb(&calc_pot, xi, col_F,
	      calc_pot.end[col_F] - .5 * calc_pot.step[col_F], &atom->gradF);
	    forces[energy_p + h] +=
	      rho_val + (atom->rho - calc_pot.end[col_F]) * atom->gradF;
#ifdef APOT
	    forces[limit_p + h] +=
	      DUMMY_WEIGHT * 10. * dsquare(atom->rho - calc_pot.end[col_F]);
#endif /* APOT */
	  }
	  /* and in-between */
	  else {
	    forces[energy_p + h] +=
	      splint_comb(&calc_pot, xi, col_F, atom->rho, &atom->gradF);
	  }
#else
	  forces[energy_p + h] +=
	    splint_comb(&calc_pot, xi, col_F, atom->rho, &atom->gradF);
#endif /* NORESCALE */
	  /* sum up rho */
	  rho_sum_loc += atom->rho;

	}			/* end S E C O N D loop over atoms */

#ifdef DIPOLE
	/* T H I R D loop: calculate whole dipole moment for every atom */
	double  rp, dp_sum;
	int   dp_converged = 0, dp_it = 0;
	double  max_diff = 10;

	while (dp_converged == 0) {
	  dp_sum = 0;
	  for (i = 0; i < inconf[h]; i++) {	/* atoms */
	    atom = conf_atoms + i + cnfstart[h] - firstatom;
	    typ1 = atom->typ;
	    if (dp_alpha[typ1]) {

	      if (dp_it) {
		/* note: mixing parameter is different from that on in IMD */
		atom->E_tot.x =
		  (1 - dp_mix) * atom->E_ind.x + dp_mix * atom->E_old.x +
		  atom->E_stat.x;
		atom->E_tot.y =
		  (1 - dp_mix) * atom->E_ind.y + dp_mix * atom->E_old.y +
		  atom->E_stat.y;
		atom->E_tot.z =
		  (1 - dp_mix) * atom->E_ind.z + dp_mix * atom->E_old.z +
		  atom->E_stat.z;
	      } else {
		atom->E_tot.x = atom->E_ind.x + atom->E_stat.x;
		atom->E_tot.y = atom->E_ind.y + atom->E_stat.y;
		atom->E_tot.z = atom->E_ind.z + atom->E_stat.z;
	      }

	      atom->p_ind.x = dp_alpha[typ1] * atom->E_tot.x + atom->p_sr.x;
	      atom->p_ind.y = dp_alpha[typ1] * atom->E_tot.y + atom->p_sr.y;
	      atom->p_ind.z = dp_alpha[typ1] * atom->E_tot.z + atom->p_sr.z;

	      atom->E_old.x = atom->E_ind.x;
	      atom->E_old.y = atom->E_ind.y;
	      atom->E_old.z = atom->E_ind.z;

	      atom->E_ind.x = 0.;
	      atom->E_ind.y = 0.;
	      atom->E_ind.z = 0.;
	    }
	  }

	  for (i = 0; i < inconf[h]; i++) {	/* atoms */
	    atom = conf_atoms + i + cnfstart[h] - firstatom;
	    typ1 = atom->typ;
	    for (j = 0; j < atom->n_neigh; j++) {	/* neighbors */
	      neigh = atom->neigh + j;
	      typ2 = neigh->typ;
	      col = neigh->col[0];
	      /* In small cells, an atom might interact with itself */
	      self = (neigh->nr == i + cnfstart[h]) ? 1 : 0;

	      if (neigh->r < dp_cut && dp_alpha[typ1] && dp_alpha[typ2]) {

		rp =
		  SPROD(conf_atoms[neigh->nr - firstatom].p_ind, neigh->dist);
		atom->E_ind.x +=
		  neigh->grad_el * (3 * rp * neigh->dist.x -
		  conf_atoms[neigh->nr - firstatom].p_ind.x);
		atom->E_ind.y +=
		  neigh->grad_el * (3 * rp * neigh->dist.y -
		  conf_atoms[neigh->nr - firstatom].p_ind.y);
		atom->E_ind.z +=
		  neigh->grad_el * (3 * rp * neigh->dist.z -
		  conf_atoms[neigh->nr - firstatom].p_ind.z);

		if (!self) {
		  rp = SPROD(atom->p_ind, neigh->dist);
		  conf_atoms[neigh->nr - firstatom].E_ind.x +=
		    neigh->grad_el * (3 * rp * neigh->dist.x - atom->p_ind.x);
		  conf_atoms[neigh->nr - firstatom].E_ind.y +=
		    neigh->grad_el * (3 * rp * neigh->dist.y - atom->p_ind.y);
		  conf_atoms[neigh->nr - firstatom].E_ind.z +=
		    neigh->grad_el * (3 * rp * neigh->dist.z - atom->p_ind.z);
		}
	      }
	    }
	  }

	  for (i = 0; i < inconf[h]; i++) {	/* atoms */
	    atom = conf_atoms + i + cnfstart[h] - firstatom;
	    typ1 = atom->typ;
	    if (dp_alpha[typ1]) {
	      dp_sum +=
		dsquare(dp_alpha[typ1] * (atom->E_old.x - atom->E_ind.x));
	      dp_sum +=
		dsquare(dp_alpha[typ1] * (atom->E_old.y - atom->E_ind.y));
	      dp_sum +=
		dsquare(dp_alpha[typ1] * (atom->E_old.z - atom->E_ind.z));
	    }
	  }

	  dp_sum /= 3 * inconf[h];
	  dp_sum = sqrt(dp_sum);

	  if (dp_it) {
	    if ((dp_sum > max_diff) || (dp_it > 50)) {
	      dp_converged = 1;
	      for (i = 0; i < inconf[h]; i++) {	/* atoms */
		atom = conf_atoms + i + cnfstart[h] - firstatom;
		typ1 = atom->typ;
		if (dp_alpha[typ1]) {
		  atom->p_ind.x =
		    dp_alpha[typ1] * atom->E_stat.x + atom->p_sr.x;
		  atom->p_ind.y =
		    dp_alpha[typ1] * atom->E_stat.y + atom->p_sr.y;
		  atom->p_ind.z =
		    dp_alpha[typ1] * atom->E_stat.z + atom->p_sr.z;
		  atom->E_ind.x = atom->E_stat.x;
		  atom->E_ind.y = atom->E_stat.y;
		  atom->E_ind.z = atom->E_stat.z;
		}
	      }
	    }
	  }

	  if (dp_sum < dp_tol) {
	    dp_converged = 1;
	  }

	  dp_it++;
	}			/* end T H I R D loop over atoms */


	/* F O U R T H  loop: calculate monopole-dipole and dipole-dipole forces */
	double  rp_i, rp_j, pp_ij, tmp_1, tmp_2;
	double  grad_1, grad_2, srval, srgrad, srval_tail, srgrad_tail, fnval_sum,
	  grad_sum;
	for (i = 0; i < inconf[h]; i++) {	/* atoms */
	  atom = conf_atoms + i + cnfstart[h] - firstatom;
	  typ1 = atom->typ;
	  k = 3 * (cnfstart[h] + i);
	  for (j = 0; j < atom->n_neigh; j++) {	/* neighbors */
	    neigh = atom->neigh + j;
	    typ2 = neigh->typ;
	    col = neigh->col[0];

	    /* In small cells, an atom might interact with itself */
	    self = (neigh->nr == i + cnfstart[h]) ? 1 : 0;
	    if (neigh->r < dp_cut && (dp_alpha[typ1] || dp_alpha[typ2])) {

	      fnval_tail = -neigh->grad_el;
	      grad_tail = -neigh->ggrad_el;

	      if (dp_b[col] && dp_c[col]) {
		shortrange_term(neigh->r, dp_b[col], dp_c[col], &srval_tail,
		  &srgrad_tail);
		srval = fnval_tail * srval_tail;
		srgrad = fnval_tail * srgrad_tail + grad_tail * srval_tail;
	      }

	      if (self) {
		fnval_tail *= 0.5;
		grad_tail *= 0.5;
	      }

	      /* monopole-dipole contributions */
	      if (charge[typ1] && dp_alpha[typ2]) {

		if (dp_b[col] && dp_c[col]) {
		  fnval_sum = fnval_tail + srval;
		  grad_sum = grad_tail + srgrad;
		} else {
		  fnval_sum = fnval_tail;
		  grad_sum = grad_tail;
		}

		rp_j =
		  SPROD(conf_atoms[neigh->nr - firstatom].p_ind, neigh->dist);
		fnval = charge[typ1] * rp_j * fnval_sum * neigh->r;
		grad_1 = charge[typ1] * rp_j * grad_sum * neigh->r2;
		grad_2 = charge[typ1] * fnval_sum;

		forces[energy_p + h] -= fnval;

		if (uf) {
		  tmp_force.x =
		    neigh->dist.x * grad_1 + conf_atoms[neigh->nr -
		    firstatom].p_ind.x * grad_2;
		  tmp_force.y =
		    neigh->dist.y * grad_1 + conf_atoms[neigh->nr -
		    firstatom].p_ind.y * grad_2;
		  tmp_force.z =
		    neigh->dist.z * grad_1 + conf_atoms[neigh->nr -
		    firstatom].p_ind.z * grad_2;
		  forces[k] -= tmp_force.x;
		  forces[k + 1] -= tmp_force.y;
		  forces[k + 2] -= tmp_force.z;
		  /* actio = reactio */
		  l = 3 * neigh->nr;
		  forces[l] += tmp_force.x;
		  forces[l + 1] += tmp_force.y;
		  forces[l + 2] += tmp_force.z;

#ifdef STRESS
		  /* calculate stresses */
		  if (us) {
		    tmp_force.x *= neigh->r;
		    tmp_force.y *= neigh->r;
		    tmp_force.z *= neigh->r;
		    stresses = stress_p + 6 * h;
		    forces[stresses] += neigh->dist.x * tmp_force.x;
		    forces[stresses + 1] += neigh->dist.y * tmp_force.y;
		    forces[stresses + 2] += neigh->dist.z * tmp_force.z;
		    forces[stresses + 3] += neigh->dist.x * tmp_force.y;
		    forces[stresses + 4] += neigh->dist.y * tmp_force.z;
		    forces[stresses + 5] += neigh->dist.z * tmp_force.x;
		  }
#endif /* STRESS */
		}
	      }


	      /* dipole-monopole contributions */
	      if (dp_alpha[typ2] && charge[typ2]) {

		if (dp_b[col] && dp_c[col]) {
		  fnval_sum = fnval_tail + srval;
		  grad_sum = grad_tail + srgrad;
		} else {
		  fnval_sum = fnval_tail;
		  grad_sum = grad_tail;
		}

		rp_i = SPROD(atom->p_ind, neigh->dist);
		fnval = charge[typ2] * rp_i * fnval_sum * neigh->r;
		grad_1 = charge[typ2] * rp_i * grad_sum * neigh->r2;
		grad_2 = charge[typ2] * fnval_sum;

		forces[energy_p + h] += fnval;

		if (uf) {
		  tmp_force.x = neigh->dist.x * grad_1 + atom->p_ind.x * grad_2;
		  tmp_force.y = neigh->dist.y * grad_1 + atom->p_ind.y * grad_2;
		  tmp_force.z = neigh->dist.z * grad_1 + atom->p_ind.z * grad_2;
		  forces[k] += tmp_force.x;
		  forces[k + 1] += tmp_force.y;
		  forces[k + 2] += tmp_force.z;
		  /* actio = reactio */
		  l = 3 * neigh->nr;
		  forces[l] -= tmp_force.x;
		  forces[l + 1] -= tmp_force.y;
		  forces[l + 2] -= tmp_force.z;

#ifdef STRESS
		  /* calculate stresses */
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


	      /* dipole-dipole contributions */
	      if (dp_alpha[typ1] && dp_alpha[typ2]) {

		pp_ij =
		  SPROD(atom->p_ind, conf_atoms[neigh->nr - firstatom].p_ind);
		tmp_1 = 3 * rp_i * rp_j;
		tmp_2 = 3 * fnval_tail / neigh->r2;

		fnval = -(tmp_1 - pp_ij) * fnval_tail;
		grad_1 = (tmp_1 - pp_ij) * grad_tail;
		grad_2 = 2 * rp_i * rp_j;

		forces[energy_p + h] += fnval;

		if (uf) {
		  tmp_force.x =
		    grad_1 * neigh->r * neigh->dist.x -
		    tmp_2 * (grad_2 * neigh->r * neigh->dist.x -
		    rp_i * neigh->r * conf_atoms[neigh->nr -
		      firstatom].p_ind.x - rp_j * neigh->r * atom->p_ind.x);
		  tmp_force.y =
		    grad_1 * neigh->r * neigh->dist.y -
		    tmp_2 * (grad_2 * neigh->r * neigh->dist.y -
		    rp_i * neigh->r * conf_atoms[neigh->nr -
		      firstatom].p_ind.y - rp_j * neigh->r * atom->p_ind.y);
		  tmp_force.z =
		    grad_1 * neigh->r * neigh->dist.z -
		    tmp_2 * (grad_2 * neigh->r * neigh->dist.z -
		    rp_i * neigh->r * conf_atoms[neigh->nr -
		      firstatom].p_ind.z - rp_j * neigh->r * atom->p_ind.z);
		  forces[k] -= tmp_force.x;
		  forces[k + 1] -= tmp_force.y;
		  forces[k + 2] -= tmp_force.z;
		  /* actio = reactio */
		  l = 3 * neigh->nr;
		  forces[l] += tmp_force.x;
		  forces[l + 1] += tmp_force.y;
		  forces[l + 2] += tmp_force.z;

#ifdef STRESS
		  /* calculate stresses */
		  if (us) {
		    tmp_force.x *= neigh->r;
		    tmp_force.y *= neigh->r;
		    tmp_force.z *= neigh->r;
		    stresses = stress_p + 6 * h;
		    forces[stresses] += neigh->dist.x * tmp_force.x;
		    forces[stresses + 1] += neigh->dist.y * tmp_force.y;
		    forces[stresses + 2] += neigh->dist.z * tmp_force.z;
		    forces[stresses + 3] += neigh->dist.x * tmp_force.y;
		    forces[stresses + 4] += neigh->dist.y * tmp_force.z;
		    forces[stresses + 5] += neigh->dist.z * tmp_force.x;
		  }
#endif /* STRESS */
		}
	      }

	    }
	  }			/* loop over neighbours */
	}			/* end F O U R T H loop over atoms */
#endif /* DIPOLE */


	/* F I F T H  loop: self energy contributions and sum-up force contributions */
	double  qq, pp;
	for (i = 0; i < inconf[h]; i++) {	/* atoms */
	  atom = conf_atoms + i + cnfstart[h] - firstatom;
	  typ1 = atom->typ;
	  k = 3 * (cnfstart[h] + i);

	  /* self energy contributions */
	  if (charge[typ1]) {
	    qq = charge[typ1] * charge[typ1];
	    fnval = dp_eps * dp_kappa * qq / sqrt(M_PI);
	    forces[energy_p + h] -= fnval;
	  }
#ifdef DIPOLE
	  if (dp_alpha[typ1]) {
	    pp = SPROD(atom->p_ind, atom->p_ind);
	    fnval = pp / (2 * dp_alpha[typ1]);
	    forces[energy_p + h] += fnval;
	  }
	  /* alternative dipole self energy including kappa-dependence */
	  //if (dp_alpha[typ1]) {
	  // pp = SPROD(atom->p_ind, atom->p_ind);
	  // fnval = kkk * pp / sqrt(M_PI);
	  // forces[energy_p + h] += fnval;
	  //}
#endif /* DIPOLE */


	  /* sum-up: whole force contributions flow into tmpsum */
	  if (uf) {
#ifdef FWEIGHT
	    /* Weigh by absolute value of force */
	    forces[k] /= FORCE_EPS + atom->absforce;
	    forces[k + 1] /= FORCE_EPS + atom->absforce;
	    forces[k + 2] /= FORCE_EPS + atom->absforce;
#endif /* FWEIGHT */
#ifdef CONTRIB
	    if (atom->contrib)
#endif /* CONTRIB */
	    tmpsum +=
	      conf_weight[h] * (dsquare(forces[k]) + dsquare(forces[k + 1]) +
	      dsquare(forces[k + 2]));
	  }

	}			/* end F I F T H loop over atoms */


	/* S I X T H  loop: EAM force */
	if (uf) {		/* only required if we calc forces */
	  for (i = 0; i < inconf[h]; i++) {
	    atom = conf_atoms + i + cnfstart[h] - firstatom;
	    k = 3 * (cnfstart[h] + i);
	    for (j = 0; j < atom->n_neigh; j++) {
	      /* loop over neighbors */
	      neigh = atom->neigh + j;
	      /* In small cells, an atom might interact with itself */
	      self = (neigh->nr == i + cnfstart[h]) ? 1 : 0;
	      col_F = paircol + ntypes + atom->typ;	/* column of F */
	      r = neigh->r;
	      /* are we within reach? */
	      if ((r < calc_pot.end[neigh->col[1]])
		|| (r < calc_pot.end[col_F - ntypes])) {
		rho_grad =
		  (r < calc_pot.end[neigh->col[1]]) ? splint_grad_dir(&calc_pot,
		  xi, neigh->slot[1], neigh->shift[1], neigh->step[1]) : 0.;
		if (atom->typ == neigh->typ)	/* use actio = reactio */
		  rho_grad_j = rho_grad;
		else
		  rho_grad_j =
		    (r < calc_pot.end[col_F - ntypes]) ? splint_grad(&calc_pot,
		    xi, col_F - ntypes, r) : 0.;
		/* now we know everything - calculate forces */
		eam_force =
		  (rho_grad * atom->gradF +
		  rho_grad_j * conf_atoms[(neigh->nr) - firstatom].gradF);
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
		/* actio = reactio */
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
	      }			/* within reach */
	    }			/* loop over neighbours */
#ifdef FWEIGHT
	    /* Weigh by absolute value of force */
	    forces[k] /= FORCE_EPS + atom->absforce;
	    forces[k + 1] /= FORCE_EPS + atom->absforce;
	    forces[k + 2] /= FORCE_EPS + atom->absforce;
#endif /* FWEIGHT */
	    /* sum up forces  */
	    tmpsum +=
	      conf_weight[h] * (dsquare(forces[k]) + dsquare(forces[k + 1]) +
	      dsquare(forces[k + 2]));
	   }
	  }        	/* end S I X T H loop over atoms */	

	/* whole energy contributions flow into tmpsum */
	forces[energy_p + h] /= (double)inconf[h];
	forces[energy_p + h] -= force_0[energy_p + h];
#ifdef COMPAT
	tmpsum += conf_weight[h] * dsquare(eweight * forces[energy_p + h]);
#else
	tmpsum += conf_weight[h] * eweight * dsquare(forces[energy_p + h]);
#endif /* COMPAT */

#ifdef STRESS
	/* whole stress contributions flow into tmpsum */
	if (uf && us) {
	  for (i = 0; i < 6; i++) {
	    forces[stress_p + 6 * h + i] /= conf_vol[h - firstconf];
	    forces[stress_p + 6 * h + i] -= force_0[stress_p + 6 * h + i];
	    tmpsum +=
#ifdef COMPAT
	      conf_weight[h] * dsquare(sweight * forces[stress_p + 6 * h + i]);
#else
	      conf_weight[h] * sweight * dsquare(forces[stress_p + 6 * h + i]);
#endif /* COMPAT */
	  }
	}
#endif /* STRESS */
	/* limiting constraints per configuration */
	tmpsum += conf_weight[h] * dsquare(forces[limit_p + h]);
      }				/* end M A I N loop over configurations */
    }				/* parallel region */
#ifdef MPI
    /* Reduce rho_sum */
    MPI_Reduce(&rho_sum_loc, &rho_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else /* MPI */
    rho_sum = rho_sum_loc;
#endif /* MPI */

    /* dummy constraints (global) */
#ifdef APOT
    /* add punishment for out of bounds (mostly for powell_lsq) */
    if (myid == 0) {
      tmpsum += apot_punish(xi_opt, forces);
    }
#endif /* APOT */

    if (myid == 0) {
      int   g;
      for (g = 0; g < ntypes; g++) {
#ifdef NORESCALE
	/* clear field */
	forces[dummy_p + ntypes + g] = 0.;	/* Free end... */
	/* NEW: Constraint on U': U'(1.)=0; */
	forces[dummy_p + g] =
	  DUMMY_WEIGHT * splint_grad(&calc_pot, xi, paircol + ntypes + g, 1.);
#else /* NOTHING */
	forces[dummy_p + ntypes + g] = 0.;	/* Free end... */
/* constraints on U`(n) */
	forces[dummy_p + g] =
	  DUMMY_WEIGHT * splint_grad(&calc_pot, xi, paircol + ntypes + g,
	  .5 * (calc_pot.begin[paircol + ntypes + g] + calc_pot.end[paircol +
	      ntypes + g]))
	  - force_0[dummy_p + g];
#endif /* NORESCALE */
	tmpsum += dsquare(forces[dummy_p + ntypes + g]);
	tmpsum += dsquare(forces[dummy_p + g]);
      }				/* loop over types */
#ifdef NORESCALE
      /* NEW: Constraint on n: <n>=1. ONE CONSTRAINT ONLY */
      /* Calculate averages */
      rho_sum /= (double)natoms;
      /* ATTN: if there are invariant potentials, things might be problematic */
      forces[dummy_p + ntypes] = DUMMY_WEIGHT * (rho_sum - 1.);
      tmpsum += dsquare(forces[dummy_p + ntypes]);
#endif /* NORESCALE */
    }				/* only root process */


    sum = tmpsum;		/* global sum = local sum  */

#ifdef MPI
    /* reduce global sum */
    sum = 0.;
    MPI_Reduce(&tmpsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    /* gather forces, energies, stresses */
    /* forces */
    MPI_Gatherv(forces + firstatom * 3, myatoms, MPI_VECTOR,
      forces, atom_len, atom_dist, MPI_VECTOR, 0, MPI_COMM_WORLD);
    /* energies */
    MPI_Gatherv(forces + natoms * 3 + firstconf, myconf, MPI_DOUBLE,	
      forces + natoms * 3, conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /* stresses */
    MPI_Gatherv(forces + natoms * 3 + nconf + 6 * firstconf, myconf, MPI_STENS,
      forces + natoms * 3 + nconf, conf_len, conf_dist, MPI_STENS, 0,
      MPI_COMM_WORLD);
    /* punishment constraints */
    MPI_Gatherv(forces + natoms * 3 + 7 * nconf + firstconf, myconf, MPI_DOUBLE,
      forces + natoms * 3 + 7 * nconf, conf_len, conf_dist, MPI_DOUBLE, 0,
      MPI_COMM_WORLD);
    /* no need to pick up dummy constraints - are already @ root */
#endif /* MPI */

    /* root process exits this function now */
    if (myid == 0) {
      fcalls++;			/* Increase function call counter */
      if (isnan(sum)) {
#ifdef DEBUG
	printf("\n--> Force is nan! <--\n\n");
#endif /* DEBUG */
	return 10e10;
      } else
	return sum;
    }

  }

  /* once a non-root process arrives here, all is done. */
  return -1.;
}


#endif /* COULOMB && EAM */
