/****************************************************************
 *
 * utils.h: potfit utilities header file
 *
 ****************************************************************
 *
 * Copyright 2002-2012
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
 *****************************************************************/

/* vector and matrix allocations */
int  *vect_int(long dim);
double *vect_double(long dim);
double **mat_double(long rowdim, long coldim);
void  free_vect_int(int *vect);
void  free_vect_double(double *vect);
void  free_mat_double(double **matrix);

/* memory management */
void  reg_for_free(void *p, char *name, ...);
void  free_all_pointers();

/* vector procuct */
vector vec_prod(vector, vector);

/* pRNG with equal or normal distribution */
double eqdist();
double normdist();

/* different power functions */
inline int isquare(int);
inline double dsquare(double);
void  power_1(double *, double *, double *);
void  power_m(int, double *, double *, double *);

#if defined APOT && defined EVO
/* quicksort for ODE */
void  quicksort(double *x, int low, int high, double **p);
int   partition(double *x, int low, int high, int index, double **p);
void  swap_population(double *a, double *b);
#endif /* APOT && EVO */
