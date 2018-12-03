/****************************************************************
 *
 * uq.h: Potential Ensemble method for Uncertainty Quantification
 *
 ****************************************************************
 *
 * Copyright 2002-2018 - the potfit development team
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

void ensemble_generation(double);


double single_param_pert_cost(double, int);
void subsection_pert(double*, double*, int, double);
double hess_bracketing(double, int, int);
double** calc_hessian(double, int);
int calc_h0_eigenvectors(double**, double, double, double**, double*);
int calc_svd(double**, double**, double*);
double generate_mc_sample(double** const, double** const, double, double, double*, int*, FILE*);
int mc_moves(double**,double*, double*, double, FILE*);

double** mat_double(int, int); /* in powell_lsq.c */
