/****************************************************************
 *
 * forces.h: General settings for force routines
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

#ifndef FORCE_H_INCLUDED
#define FORCE_H_INCLUDED

double calc_forces(double* xi_opt, double* forces, int shutdown_flag);
extern double (*g_splint)(pot_table_t*, double*, int, double);
extern double (*g_splint_grad)(pot_table_t*, double*, int, double);
extern double (*g_splint_comb)(pot_table_t*, double*, int, double, double*);

// common force initialization (force_common.cc)
void init_force_common(int is_worker);
// individual force initialization (force_XXX.cc)
void init_force(int is_worker);

void set_force_vector_pointers();
void gather_variable(double* var);
void gather_forces(double* error_sum, double* forces);

void update_splines(double* xi, int start_col, int num_col, int grad_flag);

#if defined(STIWEB)
void update_stiweb_pointers(double*);
#endif  // STIWEB

#if defined(TERSOFF)
void update_tersoff_pointers(double*);
#endif  // TERSOFF

#endif  // FORCE_H_INCLUDED
