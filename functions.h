/****************************************************************
 *
 * functions.h: potfit analytic functions header file
 *
 ****************************************************************
 *
 * Copyright 2011 Daniel Schopf, Philipp Beck
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

#ifdef APOT

#ifndef POTFIT_H
#include "potfit.h"
#endif

/* actual functions for different potentials */

void  lj_value(real, real *, real *);
void  eopp_value(real, real *, real *);
void  morse_value(real, real *, real *);
void  ms_value(real, real *, real *);
void  buck_value(real, real *, real *);
void  softshell_value(real, real *, real *);
void  eopp_exp_value(real, real *, real *);
void  meopp_value(real, real *, real *);
void  power_decay_value(real, real *, real *);
void  exp_decay_value(real, real *, real *);
void  pohlong_value(real, real *, real *);
void  parabola_value(real, real *, real *);
void  csw_value(real, real *, real *);
void  universal_value(real, real *, real *);
void  const_value(real, real *, real *);
void  sqrt_value(real, real *, real *);
void  mexp_decay_value(real, real *, real *);
void  strmm_value(real, real *, real *);
void  double_morse_value(real, real *, real *);
void  double_exp_value(real, real *, real *);
void  poly_5_value(real, real *, real *);
void  cbb_value(real, real *, real *);
void  exp_plus_value(real, real *, real *);
void  mishin_value(real, real *, real *);
void  gen_lj_value(real, real *, real *);
void  gljm_value(real, real *, real *);
void  vas_value(real, real *, real *);
void  vpair_value(real, real *, real *);
void  csw2_value(real, real *, real *);

/* template for new potential function called newpot */

/* "newpot" potential */
void  newpot_value(real, real *, real *);
/* end of template */

/* functions for analytic potential initialization */
void  apot_init(void);
void  add_potential(char *, int, fvalue_pointer);
int   apot_assign_functions(apot_table_t *);
int   apot_check_params(real *);
int   apot_parameters(char *);
real  apot_grad(real, real *, void (*function) (real, real *, real *));
real  apot_punish(real *, real *);
real  cutoff(real, real, real);

#ifdef DEBUG
void  debug_apot();
#endif /* DEBUG */

#ifdef PAIR
/* chemical potential [chempot.c] */
int   swap_chem_pot(int, int);
int   sort_chem_pot_2d(void);
real  chemical_potential(int, int *, real *);
real  chemical_potential_1d(int *, real *);
real  chemical_potential_2d(int *, real *);
real  chemical_potential_3d(int *, real *, int);
void  init_chemical_potential(int);
#endif /* PAIR */

/* functions for electrostatic calculations  */
#ifdef COULOMB
void  init_tails(void);
void  write_coulomb_table(void);
void  write_coul2imd(void);
void  ms_init(real, real *, real *, real *);
void  buck_init(real, real *, real *, real *);
void  ms_shift(real, real *, real *);
void  buck_shift(real, real *, real *);
void  elstat_value(real, real *, real *, real *);
void  elstat_shift(real, real *, real *, real *);
#endif /* COULOMB */
#ifdef DIPOLE
real  shortrange_value(real, real, real, real);
void  shortrange_term(real, real, real, real *, real *);
#endif /* DIPOLE */

#endif /* APOT */
