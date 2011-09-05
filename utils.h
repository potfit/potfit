/****************************************************************
 *
 * utils.h: potfit utilities header file
 *
 ****************************************************************
 *
 * Copyright 2002-2011 Peter Brommer, Daniel Schopf
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://www.itap.physik.uni-stuttgart.de/
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
real *vect_real(long dim);
real **mat_real(long rowdim, long coldim);
void  free_vect_int(int *vect);
void  free_vect_real(real *vect);
void  free_mat_real(real **matrix);

/* memory management */
void  reg_for_free(void *p, char *name, ...);
void  free_all_pointers();

/* vector procuct */
vector vec_prod(vector, vector);

/* pRNG with equal or normal distribution */
real  eqdist();
real  normdist();

/* different power functions */
inline int isquare(int);
inline real dsquare(real);
void  power_1(real *, real *, real *);
void  power_m(int, real *, real *, real *);

#if defined APOT && defined EVO
/* quicksort for ODE */
void  quicksort(real *x, int low, int high, real **p);
int   partition(real *x, int low, int high, int index, real **p);
void  swap_population(real *a, real *b);
#endif /* APOT && EVO */
