/****************************************************************
 *
 * forces.c: General settings for force routines
 *
 ****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
 *
 * http://potfit.sourceforge.net/
 *
 ****************************************************************
 *
 * This file is part of potfit.
 *
 * potfit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * potfit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

#include "potfit.h"

#include "force.h"
#include "memory.h"
#include "splines.h"
#include "utils.h"

double (*g_splint)(pot_table_t*, double*, int, double);
double (*g_splint_grad)(pot_table_t*, double*, int, double);
double (*g_splint_comb)(pot_table_t*, double*, int, double, double*);

/****************************************************************
  init_force_common
    called after all parameters and potentials are read
    additional assignments and initializations can be made here
****************************************************************/

void init_force_common(int is_worker)
{
  // Select correct spline interpolation and other functions
  switch (g_pot.format_type) {
    case POTENTIAL_FORMAT_UNKNOWN:
      error(1, "Unknown potential format detected! (%s:%d)\n", __FILE__,
            __LINE__);
    case POTENTIAL_FORMAT_ANALYTIC:
    case POTENTIAL_FORMAT_TABULATED_EQ_DIST:
      g_splint = splint_ed;
      g_splint_comb = splint_comb_ed;
      g_splint_grad = splint_grad_ed;
      break;
    case POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST:
      g_splint = splint_ne;
      g_splint_comb = splint_comb_ne;
      g_splint_grad = splint_grad_ne;
      break;
  }
  // for force-specific initializations please use the
  // init_force function local to the force_xxx.c file
}

/****************************************************************
  set_force_vector_pointers() assigns the correct positions to the
    pointers for the force vector

    energy_p        position of the first energy value
    stress_p        position of the first stress value
    limit_p         position of the first limiting constraint
    dummy_p         position of the first dummy constraint
    punish_par_p    position of the first parameter *punishment (APOT only)
    punish_pot_p    position of the first potential *punishment (APOT only)

    limiting constraints: (EAM and similar potentials only)
      text

    dummy constraints:
      text

    parameter punishments: (APOT only)
      These punishments are used if a certain parameter is
      out of its predefined bounds.

    potential punishments: (APOT only)
      These punishments are used if an analytic potential
      function exhibits a bad behavior.
      E.g. parameter A needs to be smaller than parameter B
****************************************************************/

void set_force_vector_pointers()
{
  g_calc.energy_p = 3 * g_config.natoms;

#if defined(STRESS)

  g_calc.stress_p = g_calc.energy_p + g_config.nconf;
#if defined(EAM) || defined(ADP) || defined(MEAM)
  g_calc.limit_p = g_calc.stress_p + 6 * g_config.nconf;
  g_calc.dummy_p = g_calc.limit_p + g_config.nconf;
#if defined(APOT)
  g_calc.punish_par_p = g_calc.dummy_p + 2 * g_param.ntypes;
  g_calc.punish_pot_p = g_calc.punish_par_p + g_pot.apot_table.total_par -
                        g_pot.apot_table.invar_pots;
#endif  // APOT
#else   // EAM || ADP || MEAM
#if defined(APOT)
  g_calc.punish_par_p = g_calc.stress_p + 6 * g_config.nconf;
  g_calc.punish_pot_p = g_calc.punish_par_p + g_pot.apot_table.total_par -
                        g_pot.apot_table.invar_pots;
#endif  // APOT
#endif  // EAM || ADP || MEAM

#else  // STRESS

#if defined(EAM) || defined(ADP) || defined(MEAM)
  g_calc.limit_p = g_calc.energy_p + g_config.nconf;
  g_calc.dummy_p = g_calc.limit_p + g_config.nconf;
#if defined(APOT)
  g_calc.punish_par_p = g_calc.dummy_p + 2 * g_param.ntypes;
  g_calc.punish_pot_p = g_calc.punish_par_p + g_pot.apot_table.total_par -
                        g_pot.apot_table.invar_pots;
#endif  // APOT
#else   // EAM || ADP || MEAM
#if defined(APOT)
  g_calc.punish_par_p = g_calc.energy_p + g_config.nconf;
  g_calc.punish_pot_p = g_calc.punish_par_p + g_pot.apot_table.total_par -
                        g_pot.apot_table.invar_pots;
#endif  // APOT
#endif  // EAM || ADP || MEAM

#endif  // STRESS
}

/****************************************************************
  gather_variable
****************************************************************/

void gather_variable(double* var)
{
#if defined(MPI)
  // Reduce variable
  double tmpvar = 0.0;
  MPI_Reduce(&var, &tmpvar, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (g_mpi.myid == 0)
    *var = tmpvar;
#endif  // MPI
}

/****************************************************************
  gather_forces
    called after all parameters and potentials are read
    additional assignments and initializations can be made here
****************************************************************/

void gather_forces(double* error_sum, double* forces)
{
#if defined(MPI)
  double tmpsum = 0.0;

  MPI_Reduce(error_sum, &tmpsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  // gather forces, energies, stresses
  if (g_mpi.myid == 0) {
    // root node already has data in place
    // forces
    MPI_Gatherv(MPI_IN_PLACE, g_mpi.myatoms, g_mpi.MPI_VECTOR, forces,
                g_mpi.atom_len, g_mpi.atom_dist, g_mpi.MPI_VECTOR, 0,
                MPI_COMM_WORLD);
    // energies
    MPI_Gatherv(MPI_IN_PLACE, g_mpi.myconf, MPI_DOUBLE,
                forces + g_calc.energy_p, g_mpi.conf_len, g_mpi.conf_dist,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if defined(STRESS)
    // stresses
    MPI_Gatherv(MPI_IN_PLACE, g_mpi.myconf, g_mpi.MPI_STENS,
                forces + g_calc.stress_p, g_mpi.conf_len, g_mpi.conf_dist,
                g_mpi.MPI_STENS, 0, MPI_COMM_WORLD);
#endif  // STRESS
#if defined(RESCALE) && (defined(EAM) || defined(ADP) || defined(MEAM))
    // punishment constraints
    MPI_Gatherv(MPI_IN_PLACE, g_mpi.myconf, MPI_DOUBLE, forces + g_calc.limit_p,
                g_mpi.conf_len, g_mpi.conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif  // RESCALE && (EAM || ADP || MEAM)
  } else {
    // forces
    MPI_Gatherv(forces + g_mpi.firstatom * 3, g_mpi.myatoms, g_mpi.MPI_VECTOR,
                forces, g_mpi.atom_len, g_mpi.atom_dist, g_mpi.MPI_VECTOR, 0,
                MPI_COMM_WORLD);
    // energies
    MPI_Gatherv(forces + g_calc.energy_p + g_mpi.firstconf, g_mpi.myconf,
                MPI_DOUBLE, forces + g_calc.energy_p, g_mpi.conf_len,
                g_mpi.conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if defined(STRESS)
    // stresses
    MPI_Gatherv(forces + g_calc.stress_p + 6 * g_mpi.firstconf, g_mpi.myconf,
                g_mpi.MPI_STENS, forces + g_calc.stress_p, g_mpi.conf_len,
                g_mpi.conf_dist, g_mpi.MPI_STENS, 0, MPI_COMM_WORLD);
#endif  // STRESS
#if defined(RESCALE) && (defined(EAM) || defined(ADP) || defined(MEAM))
    // punishment constraints
    MPI_Gatherv(forces + g_calc.limit_p + g_mpi.firstconf, g_mpi.myconf,
                MPI_DOUBLE, forces + g_calc.limit_p, g_mpi.conf_len,
                g_mpi.conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif  // RESCALE && (EAM || ADP || MEAM)
  }

  *error_sum = tmpsum;

#endif  // MPI
}

/****************************************************************
  update_splines
****************************************************************/

void update_splines(double* xi, int start_col, int num_col, int grad_flag)
{
  for (int col = start_col; col < start_col + num_col; col++) {
    int first = g_pot.calc_pot.first[col];
    double grad_left = (grad_flag & 1) ? *(xi + first - 2) : 0.0;
    double grad_right = (grad_flag & 2) ? *(xi + first - 1) : 0.0;
    switch (g_pot.format_type) {
      case POTENTIAL_FORMAT_UNKNOWN:
        error(1, "Unknown potential format detected! (%s:%d)", __FILE__, __LINE__);
      case POTENTIAL_FORMAT_ANALYTIC:
      case POTENTIAL_FORMAT_TABULATED_EQ_DIST: {
        spline_ed(g_pot.calc_pot.step[col], xi + first,
                  g_pot.calc_pot.last[col] - first + 1,
                  grad_left, grad_right, g_pot.calc_pot.d2tab + first);
        break;
      }
      case POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST: {
        spline_ne(g_pot.calc_pot.xcoord + first, xi + first,
                  g_pot.calc_pot.last[col] - first + 1,
                  grad_left, grad_right, g_pot.calc_pot.d2tab + first);
        break;
      }
    }
  }
}
