/****************************************************************
 *
 * potential.h: potfit spline functions header file
 *
 ****************************************************************
 *
 * Copyright 2011 Daniel Schopf
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
 ****************************************************************/

#ifndef POTFIT_H
#include "potfit.h"
#endif

void  spline_ed(real, real *, int, real, real, real *);
real  splint_ed(pot_table_t *, real *, int, real);
real  splint_grad_ed(pot_table_t *, real *, int, real);
real  splint_comb_ed(pot_table_t *, real *, int, real, real *);
real  splint_dir(pot_table_t *, real *, int, real, real);
real  splint_comb_dir(pot_table_t *, real *, int, real, real, real *);
real  splint_grad_dir(pot_table_t *, real *, int, real, real);
void  spline_ne(real *, real *, int, real, real, real *);
real  splint_ne(pot_table_t *, real *, int, real);
real  splint_comb_ne(pot_table_t *, real *, int, real, real *);
real  splint_grad_ne(pot_table_t *, real *, int, real);
