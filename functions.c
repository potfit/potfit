/****************************************************************
 *
 * functions.c: Routines and function calls used for analytic potentials
 *
 ****************************************************************
 *
 * Copyright 2002-2013
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

#include "potfit.h"

#include "functions.h"
#include "utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif /* M_PI */

/* macro for simplified addition of new potential functions */
#define str(s) #s
#define add_pot(a,b) add_potential(str(a),b,&a ## _value)

/****************************************************************
 *
 * initialize the function_table for analytic potentials
 *
 ****************************************************************/

void apot_init(void)
{
  int   i;

  add_pot(lj, 2);
  add_pot(eopp, 6);
  add_pot(morse, 3);
#ifdef COULOMB
  add_potential("ms", 3, &ms_shift);
  add_potential("buck", 3, &buck_shift);
#else
  add_potential("ms", 3, &ms_value);
  add_potential("buck", 3, &buck_value);
#endif /* COULOMB */
  add_pot(softshell, 2);
  add_pot(eopp_exp, 6);
  add_pot(meopp, 7);
  add_pot(power, 2);
  add_pot(power_decay, 2);
  add_pot(exp_decay, 2);
  add_pot(bjs, 3);
  add_pot(parabola, 3);
  add_pot(csw, 4);
  add_pot(universal, 4);
  add_pot(const, 1);
  add_pot(sqrt, 2);
  add_pot(mexp_decay, 3);
  add_pot(strmm, 5);
  add_pot(double_morse, 7);
  add_pot(double_exp, 5);
  add_pot(poly_5, 5);
  add_pot(kawamura, 9);
  add_pot(kawamura_mix, 12);
  add_pot(exp_plus, 3);
  add_pot(mishin, 6);
  add_pot(gen_lj, 5);
  add_pot(gljm, 12);
  add_pot(vas, 2);
  add_pot(vpair, 7);
  add_pot(csw2, 4);
  add_pot(sheng_phi1, 5);
  add_pot(sheng_phi2, 4);
  add_pot(sheng_rho, 5);
  add_pot(sheng_F, 4);

  reg_for_free(function_table.name, "function_table.name");
  reg_for_free(function_table.n_par, "function_table.n_par");
  reg_for_free(function_table.fvalue, "function_table.fvalue");
  for (i = 0; i < n_functions; i++)
    reg_for_free(function_table.name[i], "function_table.name[i]");

  return;
}

/****************************************************************
 *
 * add analytic function to function_table
 *
 ****************************************************************/

void add_potential(char *name, int parameter, fvalue_pointer fval)
{
  int   i;
  int   k = n_functions;

  /* only add potentials with unused names */
  for (i = 0; i < k; i++) {
    if (strcmp(function_table.name[i], name) == 0) {
      error(1, "There already is a potential with the name \"%s\".", name);
    }
  }

  /* allocate memory */
  function_table.name = (char **)realloc(function_table.name, (k + 1) * sizeof(char *));
  function_table.name[k] = (char *)malloc(255 * sizeof(char));
  function_table.n_par = (int *)realloc(function_table.n_par, (k + 1) * sizeof(int));
  function_table.fvalue = (fvalue_pointer *) realloc(function_table.fvalue, (k + 1) * sizeof(fvalue_pointer));
  if (function_table.name[k] == NULL || function_table.n_par == NULL || function_table.fvalue == NULL)
    error(1, "Could not allocate memory for function_table!");

  /* assign values */
  for (i = 0; i < 255; i++)
    function_table.name[k][i] = '\0';
  strncpy(function_table.name[k], name, strlen(name));
  function_table.n_par[k] = parameter;
  function_table.fvalue[k] = fval;

  n_functions++;
}

/****************************************************************
 *
 * return the number of parameters for a specific analytic potential
 *
 ****************************************************************/

int apot_parameters(char *name)
{
  int   i;

  for (i = 0; i < n_functions; i++) {
    if (strcmp(function_table.name[i], name) == 0) {
      return function_table.n_par[i];
    }
  }

  return -1;
}

/****************************************************************
 *
 * assign function pointers to corresponding functions
 *
 ****************************************************************/

int apot_assign_functions(apot_table_t *apt)
{
  int   i, j;

  for (i = 0; i < apt->number; i++) {
    for (j = 0; j < n_functions; j++) {
      if (strcmp(apt->names[i], function_table.name[j]) == 0) {
	apt->fvalue[i] = function_table.fvalue[j];
	break;
      }
      if (j == n_functions - 1)
	return -1;
    }
  }

  return 0;
}

/****************************************************************
 *
 * actual functions representing the analytic potentials
 *
 ****************************************************************/

/****************************************************************
 *
 * lennard-jones potential
 *
 * http://dx.doi.org/doi:10.1098/rspa.1924.0082
 *
 ****************************************************************/

void lj_value(double r, double *p, double *f)
{
  static double sig_d_rad6, sig_d_rad12;

  sig_d_rad6 = (p[1] * p[1]) / (r * r);
  sig_d_rad6 = sig_d_rad6 * sig_d_rad6 * sig_d_rad6;
  sig_d_rad12 = dsquare(sig_d_rad6);

  *f = 4. * p[0] * (sig_d_rad12 - sig_d_rad6);
}

/****************************************************************
 *
 * empirical oscillating pair potential (eopp)
 *
 * http://arxiv.org/abs/0802.2926v2
 *
 ****************************************************************/

void eopp_value(double r, double *p, double *f)
{
  static double x[2], y[2], power[2];

  x[0] = r;
  x[1] = r;
  y[0] = p[1];
  y[1] = p[3];

  power_m(2, power, x, y);

  *f = p[0] / power[0] + (p[2] / power[1]) * cos(p[4] * r + p[5]);
}

/****************************************************************
 *
 * morse potential
 *
 * http://dx.doi.org/doi:10.1103/PhysRev.34.57
 *
 ****************************************************************/

void morse_value(double r, double *p, double *f)
{
  *f = p[0] * (exp(-2 * p[1] * (r - p[2])) - 2. * exp(-p[1] * (r - p[2])));
}

/****************************************************************
 *
 * morse-stretch potential (without derivative!)
 *
 * http://dx.doi.org/doi:10.1063/1.1513312
 *
 ****************************************************************/

void ms_value(double r, double *p, double *f)
{
  static double x;

  x = 1. - r / p[2];

  *f = p[0] * (exp(p[1] * x) - 2. * exp((p[1] * x) / 2.));
}

/****************************************************************
 *
 * buckingham potential (without derivative!) - slightly modified
 *
 * http://dx.doi.org/doi:10.1098/rspa.1977.0049
 *
 ****************************************************************/

void buck_value(double r, double *p, double *f)
{
  static double x, y;

  x = (p[1] * p[1]) / (r * r);
  y = x * x * x;

  *f = p[0] * exp(-r / p[1]) - p[2] * y;
}

/****************************************************************
 *
 * softshell potential
 *
 *****************************************************************/

void softshell_value(double r, double *p, double *f)
{
  static double x, y;

  x = p[0] / r;
  y = p[1];

  power_1(f, &x, &y);
}

/****************************************************************
 *
 * eopp_exp potential
 *
 * http://arxiv.org/abs/0802.2926v2
 *
 ****************************************************************/

void eopp_exp_value(double r, double *p, double *f)
{
  static double power;

  power_1(&power, &r, &p[3]);

  *f = p[0] * exp(-p[1] * r) + (p[2] / power) * cos(p[4] * r + p[5]);
}

/****************************************************************
 *
 * meopp potential
 *
 * http://arxiv.org/abs/0802.2926v2
 *
 ****************************************************************/

void meopp_value(double r, double *p, double *f)
{
  static double x[2], y[2], power[2];

  x[0] = r - p[6];
  x[1] = r;
  y[0] = p[1];
  y[1] = p[3];

  power_m(2, power, x, y);

  *f = p[0] / power[0] + (p[2] / power[1]) * cos(p[4] * r + p[5]);
}

/****************************************************************
 *
 * power potential
 *
 ****************************************************************/

void power_value(double r, double *p, double *f)
{
  static double x, y, power;

  x = r;
  y = p[1];

  power_1(&power, &x, &y);

  *f = p[0] * power;
}

/****************************************************************
 *
 * power_decay potential
 *
 ****************************************************************/

void power_decay_value(double r, double *p, double *f)
{
  static double x, y, power;

  x = 1. / r;
  y = p[1];

  power_1(&power, &x, &y);

  *f = p[0] * power;
}

/****************************************************************
 *
 * exp_decay potential
 *
 ****************************************************************/

void exp_decay_value(double r, double *p, double *f)
{
  *f = p[0] * exp(-p[1] * r);
}

/****************************************************************
 *
 * bjs potential
 *
 * http://dx.doi.org/doi:10.1103/PhysRevB.37.6632
 *
 ****************************************************************/

void bjs_value(double r, double *p, double *f)
{
  static double power;

  if (r == 0)
    *f = 0;
  else {
    power_1(&power, &r, &p[1]);
    *f = p[0] * (1. - p[1] * log(r)) * power + p[2] * r;
  }
}

/****************************************************************
 *
 * parabola potential
 *
 ****************************************************************/

void parabola_value(double r, double *p, double *f)
{
  *f = (r * r) * p[0] + r * p[1] + p[2];
}

/****************************************************************
 *
 * chantasiriwan (csw) and milstein potential
 *
 * http://dx.doi.org/doi:10.1103/PhysRevB.53.14080
 *
 ****************************************************************/

void csw_value(double r, double *p, double *f)
{
  static double power;

  power_1(&power, &r, &p[3]);

  *f = (1. + p[0] * cos(p[2] * r) + p[1] * sin(p[2] * r)) / power;
}

/****************************************************************
 *
 * chantasiriwan (csw) and milstein potential - slightly modified
 *
 * http://dx.doi.org/doi:10.1103/PhysRevB.53.14080
 *
 ****************************************************************/

void csw2_value(double r, double *p, double *f)
{
  static double power;

  power_1(&power, &r, &p[3]);

  *f = (1. + p[0] * cos(p[1] * r + p[2])) / power;
}

/****************************************************************
 *
 * universal embedding function
 *
 * http://dx.doi.org/doi:10.1557/jmr.1989.1195
 *
 ****************************************************************/

void universal_value(double r, double *p, double *f)
{
  static double x[2], y[2], power[2];

  x[0] = r;
  x[1] = r;
  y[0] = p[1];
  y[1] = p[2];

  power_m(2, power, x, y);

  *f = p[0] * (p[2] / (p[2] - p[1]) * power[0] - p[1] / (p[2] - p[1]) * power[1]) + p[3] * r;
}

/****************************************************************
 *
 * constant function
 *
 ****************************************************************/

void const_value(double r, double *p, double *f)
{
  *f = *p;
}

/****************************************************************
 *
 * square root function
 *
 * http://dx.doi.org/doi:10.1080/01418618408244210
 *
 ****************************************************************/

void sqrt_value(double r, double *p, double *f)
{
  *f = p[0] * sqrt(r / p[1]);
}

/****************************************************************
 *
 * mexp_decay potential
 *
 ****************************************************************/

void mexp_decay_value(double r, double *p, double *f)
{
  *f = p[0] * exp(-p[1] * (r - p[2]));
}

/****************************************************************
 *
 * streitz-mintmire (strmm) potential
 *
 * http://dx.doi.org/doi:10.1103/PhysRevB.50.11996
 *
 ****************************************************************/

void strmm_value(double r, double *p, double *f)
{
  static double r_0;

  r_0 = r - p[4];

  *f = 2. * p[0] * exp(-p[1] / 2. * r_0) - p[2] * (1. + p[3] * r_0) * exp(-p[3] * r_0);
}

/****************************************************************
 *
 * double morse potential
 *
 * http://dx.doi.org/doi:10.1557/proc-538-535
 *
 ****************************************************************/

void double_morse_value(double r, double *p, double *f)
{
  *f =
    (p[0] * (exp(-2. * p[1] * (r - p[2])) - 2. * exp(-p[1] * (r - p[2]))) +
    p[3] * (exp(-2. * p[4] * (r - p[5])) - 2. * exp(-p[4] * (r - p[5])))) + p[6];
}

/****************************************************************
 *
 * double exp potential
 *
 * http://dx.doi.org/doi:10.1557/proc-538-535
 *
 ****************************************************************/

void double_exp_value(double r, double *p, double *f)
{
  *f = (p[0] * exp(-p[1] * dsquare(r - p[2])) + exp(-p[3] * (r - p[4])));
}

/****************************************************************
 *
 * poly 5 potential
 *
 * http://dx.doi.org/doi:10.1557/proc-538-535
 *
 ****************************************************************/

void poly_5_value(double r, double *p, double *f)
{
  static double dr;

  dr = (r - 1.) * (r - 1.);

  *f = p[0] + .5 * p[1] * dr + p[2] * (r - 1.) * dr + p[3] * (dr * dr) + p[4] * (dr * dr) * (r - 1.);
}

/****************************************************************
 *
 * kawamura potential
 *
 * http://dx.doi.org/10.1016/S0925-8388(00)00806-9
 *
 ****************************************************************/

void kawamura_value(double r, double *p, double *f)
{
  static double r6;

  r6 = r * r * r;
  r6 *= r6;

  *f = p[0] * p[1] / r + p[2] * (p[5] + p[6]) * exp((p[3] + p[4] - r) / (p[5] + p[6])) - p[7] * p[8] / r6;

  return;
}

void kawamura_mix_value(double r, double *p, double *f)
{
  static double r6;

  r6 = r * r * r;
  r6 *= r6;

  *f = p[0] * p[1] / r + p[2] * (p[5] + p[6]) * exp((p[3] + p[4] - r) / (p[5] + p[6])) - p[7] * p[8] / r6
    + p[2] * p[9] * (exp(-2 * p[10] * (r - p[11])) - 2. * exp(-p[10] * (r - p[11])));

  return;
}

/****************************************************************
 *
 * exp_plus potential
 *
 ****************************************************************/

void exp_plus_value(double r, double *p, double *f)
{
  *f = p[0] * exp(-p[1] * r) + p[2];
}

/****************************************************************
 *
 * mishin potential
 *
 * http://dx.doi.org/doi:10.1016/j.actamat.2005.05.001
 *
 ****************************************************************/

void mishin_value(double r, double *p, double *f)
{
  static double z;
  static double temp;
  static double power;

  z = r - p[3];
  temp = exp(-p[5] * r);

  power_1(&power, &z, &p[4]);

  *f = p[0] * power * temp * (1. + p[1] * temp) + p[2];
}

/****************************************************************
 *
 * gen_lj potential, generalized lennard-jones
 *
 * http://dx.doi.org/doi:10.1016/j.actamat.2005.05.001
 *
 ****************************************************************/

void gen_lj_value(double r, double *p, double *f)
{
  static double x[2], y[2], power[2];

  x[0] = r / p[3];
  x[1] = x[0];
  y[0] = p[1];
  y[1] = p[2];

  power_m(2, power, x, y);

  *f = p[0] / (p[2] - p[1]) * (p[2] / power[0] - p[1] / power[1]) + p[4];
}

/****************************************************************
 *
 * gljm potential, generalized lennard-jones + mishin potential
 *
 * http://dx.doi.org/doi:10.1016/j.actamat.2005.05.001
 *
 ****************************************************************/

void gljm_value(double r, double *p, double *f)
{
  static double x[3], y[3], power[3];

  x[0] = r / p[3];
  x[1] = x[0];
  x[2] = r - p[9];
  y[0] = p[1];
  y[1] = p[2];
  y[2] = p[10];

  power_m(3, power, x, y);

  double temp = exp(-p[11] * power[2]);

  *f =
    p[0] / (p[2] - p[1]) * (p[2] / power[0] - p[1] / power[1]) + p[4] +
    p[5] * (p[6] * power[2] * temp * (1. + p[7] * temp) + p[8]);
}

/****************************************************************
 *
 * bond-stretching function of vashishta potential (f_c)
 *
 * http://dx.doi.org/doi:10.1016/0022-3093(94)90351-4
 *
 ****************************************************************/

void vas_value(double r, double *p, double *f)
{
  *f = exp(p[0] / (r - p[1]));
}

/****************************************************************
 *
 * original pair contributions of vashishta potential
 * (V_2 without second "Coulomb"-term)
 *
 * http://dx.doi.org/doi:10.1016/0022-3093(94)90351-4
 *
 ****************************************************************/

void vpair_value(double r, double *p, double *f)
{
  static double x[7], y, z;

  y = r;
  z = p[1];

  power_1(&x[0], &y, &z);
  x[1] = r * r;
  x[2] = x[1] * x[1];
  x[3] = p[2] * p[2];
  x[4] = p[3] * p[3];
  x[5] = p[4] * x[4] + p[5] * x[3];
  x[6] = exp(-r / p[6]);

  *f = 14.4 * (p[0] / x[0] - 0.5 * (x[5] / x[2]) * x[6]);
}

/****************************************************************
 *
 * analytical fits to sheng-aluminum-EAM-potential
 *
 ****************************************************************/

void sheng_phi1_value(double r, double *p, double *f)
{
  static double x, y, z;

  x = -p[1] * r * r;
  y = r - p[4];
  z = -p[3] * y * y;

  *f = p[0] * exp(x) + p[2] * exp(z);
}

void sheng_phi2_value(double r, double *p, double *f)
{
  static double x, y, z;

  x = -p[1] * r * r;
  y = r - p[3];
  z = p[2] * p[2] + y * y;

  *f = p[0] * exp(x) + p[2] / z;
}

void sheng_rho_value(double r, double *p, double *f)
{
  static double sig_d_rad6, sig_d_rad12, x, y, power;
  static int h, k;

  h = (r > 1.45) ? 1 : 0;
  k = (r <= 1.45) ? 1 : 0;

  x = r;
  y = p[1];
  power_1(&power, &x, &y);

  sig_d_rad6 = (p[4] * p[4]) / (r * r);
  sig_d_rad6 = sig_d_rad6 * sig_d_rad6 * sig_d_rad6;
  sig_d_rad12 = dsquare(sig_d_rad6);

  *f = (p[0] * power + p[2]) * k + (4. * p[3] * (sig_d_rad12 - sig_d_rad6)) * h;
}

void sheng_F_value(double r, double *p, double *f)
{
  static double x, y, power;

  x = r;
  y = p[1];
  power_1(&power, &x, &y);

  *f = p[0] * power + p[2] * r + p[3];
}

/****************************************************************
 *
 * template for new potential function called mypotential
 * for further information plase have a look at the online documentation
 *
 * http://potfit.itap.physik.uni-stuttgart.de/
 *
 ****************************************************************/

/* template for new function */

/****************************************************************
 *
 * newpot potential
 *
 ****************************************************************/

void newpot_value(double r, double *p, double *f)
{
  *f = r * p[0] + p[1];
}

/* end of template */

/****************************************************************
 *
 * end of analytic potentials
 *
 ****************************************************************/

/****************************************************************
 *
 * function for smooth cutoff radius
 *
 ****************************************************************/

double cutoff(double r, double r0, double h)
{
  if ((r - r0) > 0)
    return 0;

  static double val;

  val = (r - r0) / h;
  val *= val;
  val *= val;

  return val / (1. + val);
}

/****************************************************************
 *
 * check analytic parameters for special conditions
 *
 ****************************************************************/

int apot_check_params(double *params)
{
  int   i, j = 2, k;

  for (i = 0; i < apot_table.number; i++) {

    /* last parameter of eopp potential is 2 pi periodic */
    if (strcmp(apot_table.names[i], "eopp") == 0) {
      k = j + 5;
      if (params[k] > 2 * M_PI)
	do {
	  params[k] -= 2 * M_PI;
	} while (params[k] > 2 * M_PI);
      if (params[k] < 0)
	do {
	  params[k] += 2 * M_PI;
	} while (params[k] < 0);
    }
    /* the third parameter of csw2 potential is 2 pi periodic */
    if (strcmp(apot_table.names[i], "csw2") == 0) {
      k = j + 2;
      if (params[k] > 2 * M_PI)
	do {
	  params[k] -= 2 * M_PI;
	} while (params[k] > 2 * M_PI);
      if (params[k] < 0)
	do {
	  params[k] += 2 * M_PI;
	} while (params[k] < 0);
    }

    /* jump to next potential */
    j += 2 + apot_parameters(apot_table.names[i]) + smooth_pot[i];
  }
  return 0;
}

/****************************************************************
 *
 * punish analytic potential for bad habits
 *
 ****************************************************************/

double apot_punish(double *params, double *forces)
{
  int   i, j;
  double x, tmpsum = 0., min, max;

  /* loop over parameters */
  for (i = 0; i < ndim; i++) {
    min = apot_table.pmin[apot_table.idxpot[i]][apot_table.idxparam[i]];
    max = apot_table.pmax[apot_table.idxpot[i]][apot_table.idxparam[i]];
    /* punishment for out of bounds */
    if (x = params[idx[i]] - min, x < 0) {
      tmpsum += APOT_PUNISH * x * x;
      forces[punish_par_p + i] = APOT_PUNISH * x * x;
    } else if (x = params[idx[i]] - max, x > 0) {
      tmpsum += APOT_PUNISH * x * x;
      forces[punish_par_p + i] = APOT_PUNISH * x * x;
    }
  }

  j = 2;
  /* loop over potentials */
  for (i = 0; i < apot_table.number; i++) {

    /* punish eta_1 < eta_2 for eopp function */
    if (strcmp(apot_table.names[i], "eopp") == 0) {
      x = params[j + 1] - params[j + 3];
      if (x < 0) {
	forces[punish_pot_p + i] = apot_punish_value * (1 + x) * (1 + x);
	tmpsum += apot_punish_value * (1 + x) * (1 + x);
      }
    }
#ifdef EAM
    /* punish m=n for universal embedding function */
    if (strcmp(apot_table.names[i], "universal") == 0) {
      x = params[j + 2] - params[j + 1];
      if (fabs(x) < 1e-6) {
	forces[punish_pot_p + i] = apot_punish_value / (x * x);
	tmpsum += apot_punish_value / (x * x);
      }
    }
#endif /* EAM */

    /* jump to next potential */
    j += 2 + apot_table.n_par[i];
  }

  return tmpsum;
}

/****************************************************************
 *
 * calculate gradient for analytic potential
 *
 ****************************************************************/

double apot_grad(double r, double *p, void (*function) (double, double *, double *))
{
  double a, b, h = 0.0001;

  function(r + h, p, &a);
  function(r - h, p, &b);

  return (a - b) / (2. * h);
}

#ifdef COULOMB

/****************************************************************
 *
 * ms potential + first derivative
 *
 ****************************************************************/

void ms_init(double r, double *pot, double *grad, double *p)
{
  static double x[4];

  x[0] = 1 - r / p[2];
  x[1] = exp(p[1] * x[0]);
  x[2] = exp(p[1] * x[0] / 2);
  x[3] = p[0] * p[1] / (r * p[2]);

  *pot = p[0] * (x[1] - 2 * x[2]);
  *grad = x[3] * (-x[1] + x[2]);
}

/****************************************************************
 *
 * buckingham potential + first derivative
 *
 ****************************************************************/

void buck_init(double r, double *pot, double *grad, double *p)
{
  static double x[3];

  x[0] = dsquare(p[1]) / dsquare(r);
  x[1] = p[2] * x[0] * x[0] * x[0];
  x[2] = p[0] * exp(-r / p[1]);

  *pot = x[2] - x[1];
  *grad = -x[2] / p[1] + 6 * x[1] / r;
}

/****************************************************************
 *
 * shifted ms potential
 *
 ****************************************************************/

void ms_shift(double r, double *p, double *f)
{
  static double pot, grad, pot_cut, grad_cut;

  ms_init(r, &pot, &grad, p);
  ms_init(dp_cut, &pot_cut, &grad_cut, p);

  *f = pot - pot_cut - (r - dp_cut) * grad_cut;
}

/****************************************************************
 *
 * shifted buckingham potential
 *
 ****************************************************************/

void buck_shift(double r, double *p, double *f)
{
  static double pot, grad, pot_cut, grad_cut;

  buck_init(r, &pot, &grad, p);
  buck_init(dp_cut, &pot_cut, &grad_cut, p);

  *f = pot - pot_cut - r * (r - dp_cut) * grad_cut;
}

/****************************************************************
 *
 * tail of electrostatic potential and first two derivatives
 *
 ****************************************************************/

void elstat_value(double r, double dp_kappa, double *ftail, double *gtail, double *ggtail)
{
  static double x[4];

  x[0] = r * r;
  x[1] = dp_kappa * dp_kappa;
  x[2] = 2 * dp_eps * dp_kappa / sqrt(M_PI);
  x[3] = exp(-x[0] * x[1]);

  *ftail = dp_eps * erfc(dp_kappa * r) / r;
  *gtail = -(*ftail + x[2] * x[3]) / x[0];
  *ggtail = (2 * x[1] * x[2] * x[3] - *gtail * 3) / x[0];
}

/****************************************************************
 *
 * shifted tail of coloumb potential
 *
 ****************************************************************/

void elstat_shift(double r, double dp_kappa, double *fnval_tail, double *grad_tail, double *ggrad_tail)
{
  static double ftail, gtail, ggtail, ftail_cut, gtail_cut, ggtail_cut;
  static double x[3];

  x[0] = r * r;
  x[1] = dp_cut * dp_cut;
  x[2] = x[0] - x[1];

  elstat_value(r, dp_kappa, &ftail, &gtail, &ggtail);
  elstat_value(dp_cut, dp_kappa, &ftail_cut, &gtail_cut, &ggtail_cut);

  *fnval_tail = ftail - ftail_cut - x[2] * gtail_cut / 2;
  *grad_tail = gtail - gtail_cut;
  *ggrad_tail = 0.;
#ifdef DIPOLE
  *fnval_tail -= x[2] * x[2] * ggtail_cut / 8;
  *grad_tail -= x[2] * ggtail_cut / 2;
  *ggrad_tail = ggtail - ggtail_cut;
#endif /* DIPOLE */
}

#endif /* COULOMB */

#ifdef DIPOLE

/****************************************************************
 *
 * short-range part of dipole moments
 *
 ****************************************************************/

double shortrange_value(double r, double a, double b, double c)
{
  static double x[5];

  x[0] = b * r;
  x[1] = x[0] * x[0];
  x[2] = x[1] * x[0];
  x[3] = x[1] * x[1];
  x[4] = 1 + x[0] + x[1] / 2 + x[2] / 6 + x[3] / 24;

  return a * c * x[4] * exp(-x[0]) / dp_eps;
}

/****************************************************************
 *
 * tail of additional short-range contribution to energy and forces
 *
 ****************************************************************/

void shortrange_term(double r, double b, double c, double *srval_tail, double *srgrad_tail)
{
  static double x[6];

  x[0] = b * r;
  x[1] = x[0] * x[0];
  x[2] = x[1] * x[0];
  x[3] = x[1] * x[1];
  x[4] = 1 + x[0] + x[1] / 2 + x[2] / 6 + x[3] / 24;
  x[5] = exp(-x[0]);

  *srval_tail = c * x[4] * x[5] / dp_eps;
  *srgrad_tail = -c * b * x[3] * x[5] / (24 * dp_eps * r);
}

#endif /* DIPOLE */

/****************************************************************
 *
 *  debug function to print the potentials
 *
 ****************************************************************/

#ifdef DEBUG

void debug_apot()
{
  int   i, j;

  fflush(stdout);
  fprintf(stderr, "\n##############################################\n");
  fprintf(stderr, "###########      DEBUG OUTPUT      ###########\n");
  fprintf(stderr, "##############################################\n");
  fprintf(stderr, "\nThere are %d potentials with a total of %d parameters.\n",
    apot_table.number, apot_table.total_par);
  for (i = 0; i < apot_table.number; i++) {
    fprintf(stderr, "\npotential #%d (type=%s, smooth=%d)\n", i + 1, apot_table.names[i], smooth_pot[i]);
    fprintf(stderr, "begin=%f end=%f\n", apot_table.begin[i], apot_table.end[i]);
    for (j = 0; j < apot_table.n_par[i]; j++) {
      fprintf(stderr, "parameter %d: name=%s value=%f min=%f max=%f\n", j + 1,
	apot_table.param_name[i][j], apot_table.values[i][j], apot_table.pmin[i][j], apot_table.pmax[i][j]);
    }
  }
#ifdef PAIR
  if (!enable_cp) {
    fprintf(stderr, "\nchemical potentials are DISABLED!\n");
  } else {
    fprintf(stderr, "\nchemical potentials:\n");
    for (i = 0; i < ntypes; i++)
      fprintf(stderr, "cp_%d=%f min=%f max=%f\n", i, apot_table.chempot[i],
	apot_table.pmin[apot_table.number][i], apot_table.pmax[apot_table.number][i]);
    if (compnodes > 0) {
      if (ntypes == 2) {
	fprintf(stderr, "composition nodes:\n");
	for (j = 0; j < compnodes; j++)
	  fprintf(stderr, "composition=%f value=%f min=%f max=%f\n",
	    compnodelist[j], apot_table.chempot[ntypes + j],
	    apot_table.pmin[apot_table.number][ntypes + j], apot_table.pmax[apot_table.number][ntypes + j]);
      }
    }
  }
#endif /* PAIR */
  exit(EXIT_FAILURE);
}

#endif /* DEBUG */
