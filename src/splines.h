/****************************************************************
 *
 * potential.h: potfit spline functions header file
 *
 ****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
 *
 * https://www.potfit.net/
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

#ifndef SPLINES_H_INCLUDED
#define SPLINES_H_INCLUDED

void spline_ed(double, double*, int, double, double, double*);
double splint_ed(pot_table_t*, double*, int, double);
double splint_grad_ed(pot_table_t*, double*, int, double);
double splint_comb_ed(pot_table_t*, double*, int, double, double*);
double splint_dir(pot_table_t*, double*, int, double, double);
double splint_comb_dir(pot_table_t*, double*, int, double, double, double*);
double splint_grad_dir(pot_table_t*, double*, int, double, double);
void spline_ne(double*, double*, int, double, double, double*);
double splint_ne(pot_table_t*, double*, int, double);
double splint_ne_lin(pot_table_t*, double*, int, double);
double splint_comb_ne(pot_table_t*, double*, int, double, double*);
double splint_grad_ne(pot_table_t*, double*, int, double);

#endif  // SPLINES_H_INCLUDED
