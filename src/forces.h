/****************************************************************
 *
 * forces.h: General settings for force routines
 *
 ****************************************************************
 *
 * Copyright 2002-2015
 *      Institute for Theoretical and Applied Physics
 *      University of Stuttgart, D-70550 Stuttgart, Germany
 *      http://potfit.sourceforge.net/
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

#ifndef POTFIT_FORCES_H
#define POTFIT_FORCES_H

double (*g_calc_forces)(double* xi_opt, double* forces, int shutdown_flag);
double (*g_splint)(pot_table_t *, double *, int, double);
double (*g_splint_grad)(pot_table_t *, double *, int, double);
double (*g_splint_comb)(pot_table_t *, double *, int, double, double *);

void  init_forces(int is_worker);
void  set_force_vector_pointers();

#endif // POTFIT_FORCES_H
