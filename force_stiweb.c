/****************************************************************
 *
 * force_stiweb.c: Routines used for calculating Stillinger-Weber
 * 	forces/energies
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

#ifdef STIWEB

#ifndef APOT
#error Stillinger-Weber needs the APOT flag to work
#endif /* APOT */

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

double calc_forces_stiweb(double *xi_opt, double *forces, int flag)
{
  int   col, i;
  double tmpsum = 0.0, sum = 0.0;
  const sw_t *sw = &apot_table.sw;

#ifndef MPI
  myconf = nconf;
#endif /* !MPI */

  /* This is the start of an infinite loop */
  while (1) {
    tmpsum = 0.0;		/* sum of squares of local process */

#ifndef MPI
  apot_check_params(xi_opt);
#endif /* !MPI */

#ifdef MPI
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (1 == flag)
      break;			/* Exception: flag 1 means clean up */

    if (0 == myid)
      apot_check_params(xi_opt);
    MPI_Bcast(xi_opt, ndimtot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif /* MPI */

    update_stiweb_pointers(xi_opt);

    /* region containing loop over configurations */
    {
      atom_t *atom;
      int   h, j, k, l;
      int   self, uf;
#ifdef STRESS
      int   us, stresses;
#endif /* STRESS */

      neigh_t *neigh_j;

      /* pair variables */
      double phi_r, phi_a, inv_c, f_cut;
      double power[2], x[2], y[2];
      double tmp, tmp_r;
      double v2_val, v2_grad;
      vector tmp_force;

      /* threebody variables */
      int jj, kk, m, ijk;
      neigh_t *neigh_k;
      angl *n_angl;
      double v3_val, tmp_grad1, tmp_grad2;
      double tmp_jj, tmp_jk, tmp_kk;
      double tmp_1, tmp_2;
      vector force_j, force_k;

      /* loop over configurations */
      for (h = firstconf; h < firstconf + myconf; h++) {
	uf = conf_uf[h - firstconf];
#ifdef STRESS
	us = conf_us[h - firstconf];
#endif /* STRESS */
	/* reset energies and stresses */
	forces[energy_p + h] = 0.;
#ifdef STRESS
	for (i = 0; i < 6; i++)
	  forces[stress_p + 6 * h + i] = 0.;
#endif /* STRESS */

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
	}
	/* end first loop */

	/* 2nd loop: calculate pair forces and energies */
	for (i = 0; i < inconf[h]; i++) {
	  atom = conf_atoms + i + cnfstart[h] - firstatom;
	  k = 3 * (cnfstart[h] + i);
	  /* loop over neighbors */
	  for (j = 0; j < atom->n_neigh; j++) {
	    neigh_j = atom->neigh + j;
	    /* In small cells, an atom might interact with itself */
	    self = (neigh_j->nr == i + cnfstart[h]) ? 1 : 0;

	    /* pair potential part */
	    col = neigh_j->col[0];
	    if (neigh_j->r < calc_pot.end[col]) {
	      /* fn value and grad are calculated in the same step */
		x[0] = neigh_j->r;
		x[1] = x[0];
		y[0] = - *(sw->p[col]);
		y[1] = - *(sw->q[col]);
		power_m(2, power, x, y);
		phi_r = *(sw->A[col]) * power[0];
		phi_a = *(sw->B[col]) * power[1];
		inv_c = 1.0 / (neigh_j->r - *(sw->a1[col]));
		f_cut = exp(*(sw->delta[col]) * inv_c);
		v2_val = (phi_r + phi_a) * f_cut;
	      if (uf) {
		v2_grad = - v2_val * *(sw->delta[col]) * inv_c * inv_c
			- f_cut * neigh_j->inv_r * (*(sw->p[col]) * phi_r + *(sw->q[col]) * phi_a);
	      }
	      /* avoid double counting if atom is interacting with a copy of itself */
	      if (self) {
		v2_val *= 0.5;
		v2_grad *= 0.5;
	      }
	      /* not double force: cohesive energy */
	      forces[energy_p + h] += v2_val;

	      if (uf) {
		tmp_force.x = neigh_j->dist.x * v2_grad;
		tmp_force.y = neigh_j->dist.y * v2_grad;
		tmp_force.z = neigh_j->dist.z * v2_grad;
		forces[k] += tmp_force.x;
		forces[k + 1] += tmp_force.y;
		forces[k + 2] += tmp_force.z;
		/* actio = reactio */
		l = 3 * neigh_j->nr;
		forces[l] -= tmp_force.x;
		forces[l + 1] -= tmp_force.y;
		forces[l + 2] -= tmp_force.z;
#ifdef STRESS
		/* also calculate pair stresses */
		if (us) {
		  tmp_force.x *= neigh_j->r;
		  tmp_force.y *= neigh_j->r;
		  tmp_force.z *= neigh_j->r;
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
	      /* calculate for later */
	    }

	    col = neigh_j->col[0];
	    if (neigh_j->r < calc_pot.end[neigh_j->col[1]]) {
	      tmp_r = neigh_j->r - *(sw->a2[col]);
	      if (tmp_r < -0.01 * *(sw->gamma[col])) {
		      tmp_r = 1.0 / tmp_r;
		      neigh_j->f = exp(*(sw->gamma[col])*tmp_r);
		      neigh_j->df = -neigh_j->f * *(sw->gamma[col]) * tmp_r * tmp_r;
	      } else {
		      neigh_j->f = 0.0;
		      neigh_j->df = 0.0;
	      }
	    }
	  }			/* loop over neighbours j */

	  /* count number of angles */
	  ijk = 0;
	  /* loop over all neighbors */
	  for (jj = 0; jj < atom->n_neigh - 1; jj++) {
            /* Get pointer to neighbor j */
            neigh_j = atom->neigh + jj;
	    /* Force location for atom j */
	    l = 3 * neigh_j->nr;
	    /* loop over remaining neighbors */
	    for (kk = jj + 1; kk < atom->n_neigh; kk++) {
	      /* Store pointer to angular part (g) */
	      n_angl = atom->angl_part + ijk;
	      /* Get pointer to neighbor k */
    	      neigh_k = atom->neigh + kk;
  	      /* Force location for atom k */
	      m = 3 * neigh_k->nr;
	      /* shortcut for types without threebody interaction */
	      if (0.0 == *(sw->lambda[atom->typ][neigh_j->typ][neigh_k->typ])) continue;
	      	/* potential term */
		tmp = n_angl->cos + 1.0 / 3.0;
		v3_val = *(sw->lambda[atom->typ][neigh_j->typ][neigh_k->typ]) * neigh_j->f * neigh_k->f * tmp * tmp;

		/* total potential */
		forces[energy_p + h] += v3_val;

		/* forces */
		tmp_grad1  = *(sw->lambda[atom->typ][neigh_j->typ][neigh_k->typ])
			* neigh_j->f * neigh_k->f * 2 * tmp;
		tmp_grad2  = *(sw->lambda[atom->typ][neigh_j->typ][neigh_k->typ]) * tmp * tmp;

		tmp_jj = 1.0 / ( neigh_j->r * neigh_j->r );
		tmp_jk = 1.0 / ( neigh_j->r * neigh_k->r );
		tmp_kk = 1.0 / ( neigh_k->r * neigh_k->r );
		tmp_1 = tmp_grad2 * neigh_j->df * neigh_k->f - tmp_grad1 * n_angl->cos * tmp_jj;
		tmp_2 = tmp_grad1 * tmp_jk;

		force_j.x = tmp_1 * neigh_j->dist.x + tmp_2 * neigh_k->dist.x;
		force_j.y = tmp_1 * neigh_j->dist.y + tmp_2 * neigh_k->dist.y;
		force_j.z = tmp_1 * neigh_j->dist.z + tmp_2 * neigh_k->dist.z;

		tmp_1 = tmp_grad2 * neigh_k->df * neigh_j->f - tmp_grad1 * n_angl->cos * tmp_kk;
		force_k.x = tmp_1 * neigh_k->dist.x + tmp_2 * neigh_j->dist.x;
		force_k.y = tmp_1 * neigh_k->dist.y + tmp_2 * neigh_j->dist.y;
		force_k.z = tmp_1 * neigh_k->dist.z + tmp_2 * neigh_j->dist.z;

		/* update force on particle i */
		forces[k] += force_j.x + force_k.x;
		forces[k + 1] += force_j.y + force_k.y;
		forces[k + 2] += force_j.z + force_k.z;

		/* update force on particle j */
		forces[l] -= force_j.x;
		forces[l + 1] -= force_j.y;
		forces[l + 2] -= force_j.z;

		/* update force on particle k */
		forces[m] -= force_k.x;
		forces[m + 1] -= force_k.y;
		forces[m + 2] -= force_k.z;

#ifdef STRESS /* Distribute stress among atoms */
		if (us) {
		  stresses = stress_p + 6 * h;
		  tmp = force_j.x * neigh_j->dist.x + force_k.x * neigh_k->dist.x;
		  forces[stresses] -= tmp;
		  tmp = force_j.y * neigh_j->dist.y + force_k.y * neigh_k->dist.y;
		  forces[stresses + 1] -= tmp;
		  tmp = force_j.z * neigh_j->dist.z + force_k.z * neigh_k->dist.z;
		  forces[stresses + 2] -= tmp;
		  tmp = 0.5 * ( force_j.y * neigh_j->dist.z + force_k.y * neigh_k->dist.z
                        + force_j.z * neigh_j->dist.y + force_k.z * neigh_k->dist.y );
		  forces[stresses + 3] -= tmp;
		  tmp =  0.5 * ( force_j.z * neigh_j->dist.x + force_k.z * neigh_k->dist.x
                         + force_j.x * neigh_j->dist.z + force_k.x * neigh_k->dist.z );
		  forces[stresses + 4] -= tmp;
		  tmp =  0.5 * ( force_j.x * neigh_j->dist.y + force_k.x * neigh_k->dist.y
                         + force_j.y * neigh_j->dist.x + force_k.y * neigh_k->dist.x );
		  forces[stresses + 5] -= tmp;
		}
#endif /* STRESS */

	    }  /* kk */
	  }  /* jj */

/*then we can calculate contribution of forces right away */
	  if (uf) {
#ifdef FWEIGHT
	    /* Weigh by absolute value of force */
	    forces[k] /= FORCE_EPS + atom->absforce;
	    forces[k + 1] /= FORCE_EPS + atom->absforce;
	    forces[k + 2] /= FORCE_EPS + atom->absforce;
#endif /* FWEIGHT */
	    /* sum up forces */
#ifdef CONTRIB
	    if (atom->contrib)
#endif /* CONTRIB */
	      tmpsum +=
		conf_weight[h] * (dsquare(forces[k]) + dsquare(forces[k + 1]) + dsquare(forces[k + 2]));
	  }			/* second loop over atoms */
	}

	/* energy contributions */
	forces[energy_p + h] /= (double)inconf[h];
	forces[energy_p + h] -= force_0[energy_p + h];
	tmpsum += conf_weight[h] * eweight * dsquare(forces[energy_p + h]);
#ifdef STRESS
	/* stress contributions */
	if (uf && us) {
	  for (i = 0; i < 6; i++) {
	    forces[stress_p + 6 * h + i] /= conf_vol[h - firstconf];
	    forces[stress_p + 6 * h + i] -= force_0[stress_p + 6 * h + i];
	    tmpsum += conf_weight[h] * sweight * dsquare(forces[stress_p + 6 * h + i]);
	  }
	}
#endif /* STRESS */
	/* limiting constraints per configuration */
      }				/* loop over configurations */
    }				/* parallel region */

    /* dummy constraints (global) */
    /* add punishment for out of bounds (mostly for powell_lsq) */
    if (myid == 0) {
      tmpsum += apot_punish(xi_opt, forces);
    }

    sum = tmpsum;		/* global sum = local sum  */
#ifdef MPI
    /* reduce global sum */
    sum = 0.;
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
    }
#endif /* MPI */

    /* root process exits this function now */
    if (myid == 0) {
      fcalls++;			/* Increase function call counter */
      if (isnan(sum)) {
#ifdef DEBUG
	printf("\n--> Force is nan! <--\n\n");
#endif /* DEBUG */
	return 10e30;
      } else
	return sum;
    }
  }  /* infinite while loop */

  /* once a non-root process arrives here, all is done. */
  return -1.0;
}

void update_stiweb_pointers(double *xi) {
  int i,j,k;
  int index = 2;
  sw_t *sw = &apot_table.sw;

  /* allocate if this has not been done */
  if (0 == sw->init) {
    sw->A = (double **)malloc(paircol * sizeof(double *));
    sw->B = (double **)malloc(paircol * sizeof(double *));
    sw->p = (double **)malloc(paircol * sizeof(double *));
    sw->q = (double **)malloc(paircol * sizeof(double *));
    sw->delta = (double **)malloc(paircol * sizeof(double *));
    sw->a1 = (double **)malloc(paircol * sizeof(double *));
    sw->gamma = (double **)malloc(paircol * sizeof(double *));
    sw->a2 = (double **)malloc(paircol * sizeof(double *));
    for (i=0;i<paircol;i++) {
      sw->A[i] = NULL;
      sw->B[i] = NULL;
      sw->p[i] = NULL;
      sw->q[i] = NULL;
      sw->delta[i] = NULL;
      sw->a1[i] = NULL;
      sw->gamma[i] = NULL;
      sw->a2[i] = NULL;
    }
    sw->lambda = (double ****)malloc(paircol * sizeof(double ***));
    for (i=0;i<ntypes;i++) {
      sw->lambda[i] = (double ***)malloc(ntypes * sizeof(double **));
      for (j=0;j<ntypes;j++) {
        sw->lambda[i][j] = (double **)malloc(ntypes * sizeof(double *));
      }
    }
   sw->init = 1;
  }

  /* update only if the address has changed */
  if (sw->A[0] != xi + index) {
    /* set the pair parameters */
    for (i = 0;i<paircol;i++) {
      sw->A[i] = xi + index++;
      sw->B[i] = xi + index++;
      sw->p[i] = xi + index++;
      sw->q[i] = xi + index++;
      sw->delta[i] = xi + index++;
      sw->a1[i] = xi + index++;
      index += 2;
    }
    for (i = 0;i<paircol;i++) {
      sw->gamma[i] = xi + index++;
      sw->a2[i] = xi + index++;
      index += 2;
    }
    /* set the lambda pointer */
    for (i=0;i<ntypes;i++)
      for (j=0;j<ntypes;j++)
	for (k=0;k<ntypes;k++)
          sw->lambda[i][j][k] = xi + index++;
  }

  return;
}

#endif /* STIWEB */
