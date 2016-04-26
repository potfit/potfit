/****************************************************************
 *
 * functions.h: potfit analytic functions header file
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

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#if defined(APOT)

/****************************************************************
  actual functions for different potentials
****************************************************************/

#define FUNCTION(name, npar) \
  void name##_value(const double r, const double* params, double* fvalue)

#include "functions.itm"

#undef FUNCTION

// functions for analytic potential initialization
void initialize_analytic_potentials(void);
void add_potential(const char* name, int npar, fvalue_pointer function);
int apot_get_num_parameters(const char* potential_name);
int apot_assign_function_pointers(apot_table_t* apot_table);
void check_correct_apot_functions(void);

// functions for analytic potential evaluation
int apot_check_params(double* params);
double apot_cutoff(const double r, const double r0, const double h);
double apot_gradient(const double r, const double* params, fvalue_pointer func);
double apot_punish(double*, double*);

#if defined(DEBUG)
void debug_apot();
#endif  // DEBUG

// special functions for electrostatic calculations
#if defined(COULOMB)
void elstat_value(double, double, double*, double*, double*);
void elstat_shift(double, double, double*, double*, double*);
void elstat_dsf(double, double, double*, double*, double*);
void init_tails(double);
#endif  // COULOMB
#if defined(DIPOLE)
double shortrange_value(double, double, double, double);
void shortrange_term(double, double, double, double*, double*);
#endif  // DIPOLE

#endif  // APOT

#endif  // FUNCTIONS_H_INCLUDED
