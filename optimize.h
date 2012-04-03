/****************************************************************
 *
 * optimize.h: potfit optimization header file
 *
 ****************************************************************
 *
 * Copyright 2011-2012
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
 ****************************************************************/

#ifndef POTFIT_H
#include "potfit.h"
#endif

#ifdef EVO
/* differential evolution [diff_evo.c] */
void  init_population(double **, double *, double *);
#ifdef APOT
void  opposite_check(double **, double *, int);
#endif /* APOT */
void  diff_evo(double *);
#else /* EVO */
/* simulated annealing [simann.c] */
#ifdef APOT
void  randomize_parameter(int, double *, double *);
#else
void  makebump(double *, double, double, int);
#endif /* APOT */
void  anneal(double *);
#endif /* EVO */

/* powell least squares [powell_lsq.c] */
void  powell_lsq(double *);
int   gamma_init(double **, double **, double *, double *);
int   gamma_update(double **, double, double, double *, double *, double *, int,
  int, int, double);
void  lineqsys_init(double **, double **, double *, double *, int, int);
void  lineqsys_update(double **, double **, double *, double *, int, int, int);
void  copy_matrix(double **, double **, int, int);
void  copy_vector(double *, double *, int);
void  matdotvec(double **, double *, double *, int, int);
double normalize_vector(double *, int);
