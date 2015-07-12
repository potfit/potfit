/****************************************************************
 *
 * forces.c: General settings for force routines
 *
 ****************************************************************
 *
 * Copyright 2002-2015
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

#include "forces.h"
#include "functions.h"
#include "memory.h"
#include "splines.h"
#include "utils.h"

double (*g_calc_forces)(double* xi_opt, double* forces, int shutdown_flag);
double (*g_splint)(pot_table_t*, double*, int, double);
double (*g_splint_grad)(pot_table_t*, double*, int, double);
double (*g_splint_comb)(pot_table_t*, double*, int, double, double*);

double calc_forces_pair(double* xi_opt, double* forces, int shutdown_flag);
double calc_forces_eam(double* xi_opt, double* forces, int shutdown_flag);
double calc_forces_adp(double* xi_opt, double* forces, int shutdown_flag);
double calc_forces_meam(double* xi_opt, double* forces, int shutdown_flag);
double calc_forces_elstat(double* xi_opt, double* forces, int shutdown_flag);
double calc_forces_eam_elstat(double* xi_opt, double* forces, int shutdown_flag);
double calc_forces_stiweb(double* xi_opt, double* forces, int shutdown_flag);
double calc_forces_tersoff(double* xi_opt, double* forces, int shutdown_flag);
double calc_forces_tersoffmod(double* xi_opt, double* forces, int shutdown_flag);

/****************************************************************
 *
 *  init_forces
 *      called after all parameters and potentials are read
 *      additional assignments and initializations can be made here
 *
 ****************************************************************/

void init_forces(int is_worker)
{
// set proper force routine
#if defined(PAIR)
  g_calc_forces = &calc_forces_pair;
#endif

#if defined(EAM) && !defined(COULOMB) && !defined(DIPOLE)
  g_calc_forces = &calc_forces_eam;
#endif

#if defined(ADP)
  g_calc_forces = &calc_forces_adp;
#endif

#if defined(COULOMB) || defined(DIPOLE)
#if defined(EAM)
  g_calc_forces = &calc_forces_eam_elstat;
#else
  g_calc_forces = &calc_forces_elstat;
#endif
#endif

#if defined(MEAM)
  g_calc_forces = &calc_forces_meam;
#endif

#if defined(STIWEB)
  g_calc_forces = &calc_forces_stiweb;
#endif

#if defined(TERSOFF)
#if defined(TERSOFFMOD)
  g_calc_forces = &calc_forces_tersoffmod;
#else
  g_calc_forces = &calc_forces_tersoff;
#endif
#endif

/* Select correct spline interpolation and other functions */
#if defined(APOT)
  if (g_pot.format == 0)
  {
    g_splint = splint_ed;
    g_splint_comb = splint_comb_ed;
    g_splint_grad = splint_grad_ed;
  }
#else
  if (g_pot.format == 3)
  {
    g_splint = splint_ed;
    g_splint_comb = splint_comb_ed;
    g_splint_grad = splint_grad_ed;
  }
  else if (g_pot.format >= 4)
  {
    g_splint = splint_ne;
    g_splint_comb = splint_comb_ne;
    g_splint_grad = splint_grad_ne;
  }
#endif /* APOT */

#ifdef COULOMB
  if (g_pot.apot_table.sw_kappa)
    init_tails(g_pot.apot_table.dp_kappa[0]);
#endif /* COULOMB */

/* set spline density corrections to 0 */
#if defined(EAM) || defined(ADP) || defined(MEAM)
  g_todo.lambda = (double*)Malloc(g_param.ntypes * sizeof(double));
#endif /* EAM || ADP || MEAM */
}

/****************************************************************
 *
 *  set_force_vector_pointers() assigns the correct positions to the
 *  	pointers for the force vector
 *
 *  	energy_p 	... 	position of the first energy value
 *  	stress_p 	... 	position of the first stress value
 *  	limit_p 	... 	position of the first limiting constraint
 *  	dummy_p 	... 	position of the first dummy constraint
 *  	punish_par_p 	... 	position of the first parameter
 *punishment (APOT only)
 *  	punish_pot_p 	... 	position of the first potential
 *punishment (APOT only)
 *
 * 	limiting constraints: (EAM and similar potentials only)
 * 		text
 *
 * 	dummy constraints:
 * 		test
 *
 * 	parameter punishments: (APOT only)
 * 		These punishments are used if a certain parameter is
 * 		out of its predefined bounds.
 *
 * 	potential punishments: (APOT only)
 * 		These punishments are used if an analytic potential
 * 		function exhibits a bad behavior.
 * 		E.g. parameter A needs to be smaller than parameter B
 *
 ****************************************************************/

void set_force_vector_pointers()
{
  g_calc.energy_p = 3 * g_config.natoms;

#ifdef STRESS

  g_calc.stress_p = g_calc.energy_p + g_config.nconf;
#if defined EAM || defined ADP || defined MEAM
  g_calc.limit_p = g_calc.stress_p + 6 * g_config.nconf;
  g_calc.dummy_p = g_calc.limit_p + g_config.nconf;
#ifdef APOT
  g_calc.punish_par_p = g_calc.dummy_p + 2 * g_param.ntypes;
  g_calc.punish_pot_p =
      g_calc.punish_par_p + g_pot.apot_table.total_par - g_pot.apot_table.invar_pots;
#endif /* APOT */
#else  /* EAM || ADP || MEAM */
#ifdef APOT
  g_calc.punish_par_p = g_calc.stress_p + 6 * g_config.nconf;
  g_calc.punish_pot_p =
      g_calc.punish_par_p + g_pot.apot_table.total_par - g_pot.apot_table.invar_pots;
#endif /* APOT */
#endif /* EAM || ADP || MEAM */

#else /* STRESS */

#if defined EAM || defined ADP || defined MEAM
  g_calc.limit_p = g_calc.energy_p + g_config.nconf;
  g_calc.dummy_p = g_calc.limit_p + g_config.nconf;
#ifdef APOT
  g_calc.punish_par_p = g_calc.dummy_p + 2 * g_param.ntypes;
  g_calc.punish_pot_p =
      g_calc.punish_par_p + g_pot.apot_table.total_par - g_pot.apot_table.invar_pots;
#endif /* APOT */
#else  /* EAM || ADP || MEAM */
#ifdef APOT
  g_calc.punish_par_p = g_calc.energy_p + g_config.nconf;
  g_calc.punish_pot_p =
      g_calc.punish_par_p + g_pot.apot_table.total_par - g_pot.apot_table.invar_pots;
#endif /* APOT */
#endif /* EAM || ADP || MEAM */

#endif /* STRESS */
}

/****************************************************************
 *
 *  gather_forces
 *      called after all parameters and potentials are read
 *      additional assignments and initializations can be made here
 *
 ****************************************************************/

void gather_forces(double* error_sum, double* forces)
{
#if defined(MPI)
  double tmpsum = 0.0;

  MPI_Reduce(error_sum, &tmpsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* gather forces, energies, stresses */
  if (g_mpi.myid == 0)
  {
    /* root node already has data in place */
    /* forces */
    MPI_Gatherv(MPI_IN_PLACE, g_mpi.myatoms, g_mpi.MPI_VECTOR, forces, g_mpi.atom_len,
                g_mpi.atom_dist, g_mpi.MPI_VECTOR, 0, MPI_COMM_WORLD);
    /* energies */
    MPI_Gatherv(MPI_IN_PLACE, g_mpi.myconf, MPI_DOUBLE, forces + g_calc.energy_p,
                g_mpi.conf_len, g_mpi.conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if defined(STRESS)
    /* stresses */
    MPI_Gatherv(MPI_IN_PLACE, g_mpi.myconf, g_mpi.MPI_STENS, forces + g_calc.stress_p,
                g_mpi.conf_len, g_mpi.conf_dist, g_mpi.MPI_STENS, 0, MPI_COMM_WORLD);
#endif  // STRESS
  }
  else
  {
    /* forces */
    MPI_Gatherv(forces + g_mpi.firstatom * 3, g_mpi.myatoms, g_mpi.MPI_VECTOR, forces,
                g_mpi.atom_len, g_mpi.atom_dist, g_mpi.MPI_VECTOR, 0, MPI_COMM_WORLD);
    /* energies */
    MPI_Gatherv(forces + g_calc.energy_p + g_mpi.firstconf, g_mpi.myconf, MPI_DOUBLE,
                forces + g_calc.energy_p, g_mpi.conf_len, g_mpi.conf_dist, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
#if defined(STRESS)
    /* stresses */
    MPI_Gatherv(forces + g_calc.stress_p + 6 * g_mpi.firstconf, g_mpi.myconf,
                g_mpi.MPI_STENS, forces + g_calc.stress_p, g_mpi.conf_len,
                g_mpi.conf_dist, g_mpi.MPI_STENS, 0, MPI_COMM_WORLD);
#endif  // STRESS
  }

  *error_sum = tmpsum;

#else

  return;

#endif  // MPI
}
