/****************************************************************
 *
 * functions.h: potfit analytic functions header file
 *
 ****************************************************************
 *
 * Copyright 2002-2013
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
 ****************************************************************/

#ifdef APOT

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

/* actual functions for different potentials */

void  lj_value(double, double *, double *);
void  eopp_value(double, double *, double *);
void  morse_value(double, double *, double *);
void  ms_value(double, double *, double *);
void  buck_value(double, double *, double *);
void  softshell_value(double, double *, double *);
void  eopp_exp_value(double, double *, double *);
void  meopp_value(double, double *, double *);
void  power_value(double, double *, double *);
void  power_decay_value(double, double *, double *);
void  exp_decay_value(double, double *, double *);
void  bjs_value(double, double *, double *);
void  parabola_value(double, double *, double *);
void  csw_value(double, double *, double *);
void  universal_value(double, double *, double *);
void  const_value(double, double *, double *);
void  sqrt_value(double, double *, double *);
void  mexp_decay_value(double, double *, double *);
void  strmm_value(double, double *, double *);
void  double_morse_value(double, double *, double *);
void  double_exp_value(double, double *, double *);
void  poly_5_value(double, double *, double *);
void  kawamura_value(double, double *, double *);
void  kawamura_mix_value(double, double *, double *);
void  exp_plus_value(double, double *, double *);
void  mishin_value(double, double *, double *);
void  gen_lj_value(double, double *, double *);
void  gljm_value(double, double *, double *);
void  vas_value(double, double *, double *);
void  vpair_value(double, double *, double *);
void  csw2_value(double, double *, double *);
void  sheng_phi1_value(double, double *, double *);
void  sheng_phi2_value(double, double *, double *);
void  sheng_rho_value(double, double *, double *);
void  sheng_F_value(double, double *, double *);

#ifdef STIWEB
void  stiweb_2_value(double, double *, double *);
void  stiweb_3_value(double, double *, double *);
void  lambda_value(double, double *, double *);
#endif /* STIWEB */

#ifdef TERSOFF
void  tersoff_pot_value(double, double *, double *);
void  tersoff_mix_value(double, double *, double *);
#endif /* TERSOFF */

/* template for new potential function called newpot */

/* "newpot" potential */
void  newpot_value(double, double *, double *);

/* end of template */

/* functions for analytic potential initialization */
void  apot_init(void);
void  add_potential(char *, int, fvalue_pointer);
int   apot_assign_functions(apot_table_t *);
int   apot_check_params(double *);
int   apot_parameters(char *);
void  check_apot_functions(void);
double apot_grad(double, double *, void (*function) (double, double *, double *));
double apot_punish(double *, double *);
double cutoff(double, double, double);

#ifdef DEBUG
void  debug_apot();
#endif /* DEBUG */

/* chemical potential [chempot.c] */
#ifdef PAIR
int   swap_chem_pot(int, int);
int   sort_chem_pot_2d(void);
double chemical_potential(int, int *, double *);
double chemical_potential_1d(int *, double *);
double chemical_potential_2d(int *, double *);
double chemical_potential_3d(int *, double *, int);
void  init_chemical_potential(int);
#endif /* PAIR */

/* functions for electrostatic calculations  */
#ifdef COULOMB
void  ms_init(double, double *, double *, double *);
void  buck_init(double, double *, double *, double *);
void  ms_shift(double, double *, double *);
void  buck_shift(double, double *, double *);
void  elstat_value(double, double, double *, double *, double *);
void  elstat_shift(double, double, double *, double *, double *);
#endif /* COULOMB */
#ifdef DIPOLE
double shortrange_value(double, double, double, double);
void  shortrange_term(double, double, double, double *, double *);
#endif /* DIPOLE */

#endif /* FUNCTIONS_H */

#endif /* APOT */
