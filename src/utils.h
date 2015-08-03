/****************************************************************
 *
 * utils.h: potfit utilities header file
 *
 ****************************************************************
 *
 * Copyright 2002-2014
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
 *****************************************************************/

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

/* vector and matrix allocations */
// int* vect_int(long);
// double* vect_double(long);
double** mat_double(long, long);
void free_vect_int(int*);
void free_vect_double(double*);
void free_mat_double(double**);

/* vector procuct */
vector vec_prod(vector, vector);

/* different power functions */
static inline int isquare(int a) { return a * a; }
static inline double dsquare(double a) { return a * a; }
void power_1(double* result, double* base, double* exponent);
void power_m(int count, double* result, double* base, double* exponent);

#endif  // UTILS_H_INCLUDED
