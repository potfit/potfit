/****************************************************************
 *
 * optimize.h: potfit optimization header file
 *
 ****************************************************************
 *
 * Copyright 2011
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
void  init_population(real **, real *, real *);
#ifdef APOT
void  opposite_check(real **, real *, int);
#endif /* APOT */
void  diff_evo(real *);
#else /* EVO */
/* simulated annealing [simann.c] */
#ifdef APOT
void  randomize_parameter(int, real *, real *);
#else
void  makebump(real *, real, real, int);
#endif /* APOT */
void  anneal(real *);
#endif /* EVO */

/* powell least squares [powell_lsq.c] */
void  powell_lsq(real *);
int   gamma_init(real **, real **, real *, real *);
int   gamma_update(real **, real, real, real *, real *, real *, int, int, int,
  real);
void  lineqsys_init(real **, real **, real *, real *, int, int);
void  lineqsys_update(real **, real **, real *, real *, int, int, int);
void  copy_matrix(real **, real **, int, int);
void  copy_vector(real *, real *, int);
void  matdotvec(real **, real *, real *, int, int);
real  normalize_vector(real *, int);
