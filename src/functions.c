/****************************************************************
 *
 * functions.c: Routines and function calls used for analytic potentials
 *
 ****************************************************************
 *
 * Copyright 2002-2015
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

#include "potfit.h"

#include "functions.h"
#include "memory.h"
#include "utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif  // M_PI

/* macro for simplified addition of new potential functions */
#define ADD_POT(a, b) add_potential(#a, b, &a##_value)

/****************************************************************
 *
 *  function table: stores all available analytic potentials
 *
 ****************************************************************/

struct
{
  char** name;            /* identifier of the potential */
  int* num_params;        /* number of parameters */
  fvalue_pointer* fvalue; /* function pointer */
  int num_functions;      /* number of analytic function prototypes */
} function_table;

/****************************************************************
 *
 * initialize the function_table for analytic potentials
 *
 ****************************************************************/

void apot_init(void)
{
  ADD_POT(lj, 2);
  ADD_POT(eopp, 6);
  ADD_POT(morse, 3);
#ifdef COULOMB
  add_potential("ms", 3, &ms_shift);
  add_potential("buck", 3, &buck_shift);
#else
  add_potential("ms", 3, &ms_value);
  add_potential("buck", 3, &buck_value);
#endif /* COULOMB */
  ADD_POT(softshell, 2);
  ADD_POT(eopp_exp, 6);
  ADD_POT(meopp, 7);
  ADD_POT(power, 2);
  ADD_POT(power_decay, 2);
  ADD_POT(exp_decay, 2);
  ADD_POT(bjs, 3);
  ADD_POT(parabola, 3);
  ADD_POT(csw, 4);
  ADD_POT(universal, 4);
  ADD_POT(const, 1);
  ADD_POT(sqrt, 2);
  ADD_POT(mexp_decay, 3);
  ADD_POT(strmm, 5);
  ADD_POT(double_morse, 7);
  ADD_POT(double_exp, 5);
  ADD_POT(poly_5, 5);
  ADD_POT(kawamura, 9);
  ADD_POT(kawamura_mix, 12);
  ADD_POT(exp_plus, 3);
  ADD_POT(mishin, 6);
  ADD_POT(gen_lj, 5);
  ADD_POT(gljm, 12);
  ADD_POT(vas, 2);
  ADD_POT(vpair, 7);
  ADD_POT(csw2, 4);
  ADD_POT(sheng_phi1, 5);
  ADD_POT(sheng_phi2, 4);
  ADD_POT(sheng_rho, 5);
  ADD_POT(sheng_F, 4);

#if defined(STIWEB)
  ADD_POT(stiweb_2, 6);
  ADD_POT(stiweb_3, 2);
  ADD_POT(lambda, (int)(0.5 * g_param.ntypes * g_param.ntypes * (g_param.ntypes + 1)));
#endif  // STIWEB

#if defined(TERSOFF)
#if !defined(TERSOFFMOD)
  ADD_POT(tersoff_pot, 11);
  ADD_POT(tersoff_mix, 2);
#else
  ADD_POT(tersoff_mod_pot, 16);
#endif  // !TERSOFFMOD
#endif  // TERSOFF
}

/****************************************************************
 *
 * add analytic function to function_table
 *
 ****************************************************************/

void add_potential(const char* name, int parameter, fvalue_pointer fval)
{
  const int k = function_table.num_functions;

  /* only add potentials with unused names */
  for (int i = 0; i < k; i++)
  {
    if (strcmp(function_table.name[i], name) == 0)
    {
      error(1, "There already is a potential with the name \"%s\".", name);
    }
  }

  /* allocate memory */
  function_table.name = (char**)Realloc(function_table.name, (k + 1) * sizeof(char*));
  function_table.name[k] = (char*)Malloc((strlen(name) + 1) * sizeof(char));
  function_table.num_params =
      (int*)Realloc(function_table.num_params, (k + 1) * sizeof(int));
  function_table.fvalue =
      (fvalue_pointer*)Realloc(function_table.fvalue, (k + 1) * sizeof(fvalue_pointer));

  /* assign values */
  strncpy(function_table.name[k], name, strlen(name));
  function_table.num_params[k] = parameter;
  function_table.fvalue[k] = fval;

  function_table.num_functions++;
}

/****************************************************************
 *
 * return the number of parameters for a specific analytic potential
 *
 ****************************************************************/

int apot_parameters(char* name)
{
  int i;

  for (i = 0; i < function_table.num_functions; i++)
  {
    if (strcmp(function_table.name[i], name) == 0)
    {
      return function_table.num_params[i];
    }
  }

  return -1;
}

/****************************************************************
 *
 * assign function pointers to corresponding functions
 *
 ****************************************************************/

int apot_assign_functions(apot_table_t* apt)
{
  for (int i = 0; i < apt->number; i++)
  {
    for (int j = 0; j < function_table.num_functions; j++)
    {
      if (strcmp(apt->names[i], function_table.name[j]) == 0)
      {
        apt->fvalue[i] = function_table.fvalue[j];
        break;
      }
      if (j == function_table.num_functions - 1)
        return -1;
    }
  }

  return 0;
}

/****************************************************************
 *
 * check for special functions needed by certain potential models
 *
 ****************************************************************/

void check_apot_functions(void)
{
#if defined(STIWEB)
  /* paircol is not yet defined at this point */
  int pcol = (g_param.ntypes * (g_param.ntypes + 1)) / 2;

  /* check for the correct function types for SW potential */
  for (int i = 0; i < pcol; i++)
  {
    if (strcmp(g_pot.apot_table.names[i], "stiweb_2") != 0)
      error(1, "Only stiweb_2 potential is allowed for the %d. potential!\n", i + 1);
    if (strcmp(g_pot.apot_table.names[pcol + i], "stiweb_3") != 0)
      error(1, "Only stiweb_3 potential is allowed for the %d. potential!\n",
            pcol + i + 1);
  }

  if (strcmp(g_pot.apot_table.names[2 * pcol], "lambda") != 0)
    error(1,
          "The last potential for Stillinger-Weber has to be of the \"lambda\" "
          "type!\n");

  /* make sure the cutoff parameters (a1,a2) can be optimized correctly */
  for (int i = 0; i < pcol; i++)
  {
    if (g_pot.apot_table.pmax[i][5] > g_pot.apot_table.end[i])
    {
      error(0,
            "The upper bound for the parameter a1 exceeds the cutoff radius in "
            "potential %d.\n",
            i + 1);
      error(1, "a1 needs to be less or equal to the potential cutoff.\n");
    }
    if (g_pot.apot_table.pmax[pcol + i][1] > g_pot.apot_table.end[pcol + i])
    {
      error(0,
            "The upper bound for the parameter a2 exceeds the cutoff radius in "
            "potential %d.\n",
            i + 1);
      error(1, "a1 needs to be less or equal to the potential cutoff.\n");
    }
  }
#endif // STIWEB

#if defined(TERSOFF)
  /* paircol is not yet defined at this point */
  int pcol = (g_param.ntypes * (g_param.ntypes + 1)) / 2;

#if !defined(TERSOFFMOD)
  /* check for the correct function types for TERSOFF potential */
  for (int i = 0; i < pcol; i++)
  {
    if (strcmp(g_pot.apot_table.names[i], "tersoff_pot") != 0)
      error(1, "Only tersoff_pot potential is allowed for the %d. potential!\n", i + 1);
  }
  for (int i = 0; i < g_param.ntypes * (g_param.ntypes - 1) / 2.0; i++)
  {
    if (strcmp(g_pot.apot_table.names[pcol + i], "tersoff_mix") != 0)
      error(1, "Only tersoff_mix potential is allowed for the %d. potential!\n", i + 1);
  }

  /* make sure the cutoff parameters (R, S) can be optimized correctly */
  for (int i = 0; i < pcol; i++)
  {
    if (g_pot.apot_table.pmax[i][9] < g_pot.apot_table.pmin[i][10])
    {
      error(0,
            "The upper bound for the parameter S is smaller than the lower "
            "bound\n");
      error(0, "for the parameter R in potential %d.\n", i + 1);
      error(1, "Please change it, that the condition R < S can be fulfilled.\n");
    }
  }
#else
  /* check for the correct function types for TERSOFFMOD potential */
  for (int i = 0; i < pcol; i++)
  {
    if (strcmp(g_pot.apot_table.names[i], "tersoff_mod_pot") != 0)
      error(1, "Only tersoff_pot potential is allowed for the %d. potential!\n", i + 1);
  }

  /* make sure the cutoff parameters (R, S) can be optimized correctly */
  for (int i = 0; i < pcol; i++)
  {
    if (g_pot.apot_table.pmax[i][15] < g_pot.apot_table.pmin[i][14])
    {
      error(0,
            "The upper bound for the parameter R2 is smaller than the lower "
            "bound\n");
      error(0, "for the parameter R1 in potential %d.\n", i + 1);
      error(1, "Please change it, that the condition R1 < R2 can be fulfilled.\n");
    }
  }

#endif // !TERSOFFMOD
#endif // TERSOFF
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

void lj_value(double r, double* p, double* f)
{
  double x = (p[1] * p[1]) / (r * r);
  x = x * x * x;

  *f = 4.0 * p[0] * x * (x - 1.0);
}

/****************************************************************
 *
 * empirical oscillating pair potential (eopp)
 *
 * http://arxiv.org/abs/0802.2926v2
 *
 ****************************************************************/

void eopp_value(double r, double* p, double* f)
{
  double x[2] = {r, r};
  double y[2] = {p[1], p[3]};
  double power[2] = {0, 0};

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

void morse_value(double r, double* p, double* f)
{
  *f = p[0] * (exp(-2.0 * p[1] * (r - p[2])) - 2.0 * exp(-p[1] * (r - p[2])));
}

/****************************************************************
 *
 * morse-stretch potential (without derivative!)
 *
 * http://dx.doi.org/doi:10.1063/1.1513312
 *
 ****************************************************************/

void ms_value(double r, double* p, double* f)
{
  double x = 1.0 - r / p[2];

  *f = p[0] * (exp(p[1] * x) - 2.0 * exp((p[1] * x) / 2.0));
}

/****************************************************************
 *
 * buckingham potential (without derivative!) - slightly modified
 *
 * http://dx.doi.org/doi:10.1098/rspa.1977.0049
 *
 ****************************************************************/

void buck_value(double r, double* p, double* f)
{
  double x = (p[1] * p[1]) / (r * r);
  double y = x * x * x;

  *f = p[0] * exp(-r / p[1]) - p[2] * y;
}

/****************************************************************
 *
 * softshell potential
 *
 *****************************************************************/

void softshell_value(double r, double* p, double* f)
{
  double x = p[0] / r;
  double y = p[1];

  power_1(f, &x, &y);
}

/****************************************************************
 *
 * eopp_exp potential
 *
 * http://arxiv.org/abs/0802.2926v2
 *
 ****************************************************************/

void eopp_exp_value(double r, double* p, double* f)
{
  double power = 0;

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

void meopp_value(double r, double* p, double* f)
{
  double x[2] = {r-p[6], r};
  double y[2] = {p[1], p[3]};
  double power[2] = {0, 0};

  power_m(2, power, x, y);

  *f = p[0] / power[0] + (p[2] / power[1]) * cos(p[4] * r + p[5]);
}

/****************************************************************
 *
 * power potential
 *
 ****************************************************************/

void power_value(double r, double* p, double* f)
{
  double power = 0;

  power_1(&power, &r, &p[1]);

  *f = p[0] * power;
}

/****************************************************************
 *
 * power_decay potential
 *
 ****************************************************************/

void power_decay_value(double r, double* p, double* f)
{
  double power = 0;

  power_1(&power, &r, &p[1]);

  *f = p[0] / power;
}

/****************************************************************
 *
 * exp_decay potential
 *
 ****************************************************************/

void exp_decay_value(double r, double* p, double* f) { *f = p[0] * exp(-p[1] * r); }

/****************************************************************
 *
 * bjs potential
 *
 * http://dx.doi.org/doi:10.1103/PhysRevB.37.6632
 *
 ****************************************************************/

void bjs_value(double r, double* p, double* f)
{
  if (r == 0.0)
    *f = 0.0;
  else
  {
    double power = 0;

    power_1(&power, &r, &p[1]);

    *f = p[0] * (1.0 - p[1] * log(r)) * power + p[2] * r;
  }
}

/****************************************************************
 *
 * parabola potential
 *
 ****************************************************************/

void parabola_value(double r, double* p, double* f)
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

void csw_value(double r, double* p, double* f)
{
  static double power;

  power_1(&power, &r, &p[3]);

  *f = (1.0 + p[0] * cos(p[2] * r) + p[1] * sin(p[2] * r)) / power;
}

/****************************************************************
 *
 * chantasiriwan (csw) and milstein potential - slightly modified
 *
 * http://dx.doi.org/doi:10.1103/PhysRevB.53.14080
 *
 ****************************************************************/

void csw2_value(double r, double* p, double* f)
{
  static double power;

  power_1(&power, &r, &p[3]);

  *f = (1.0 + p[0] * cos(p[1] * r + p[2])) / power;
}

/****************************************************************
 *
 * universal embedding function
 *
 * http://dx.doi.org/doi:10.1557/jmr.1989.1195
 *
 ****************************************************************/

void universal_value(double r, double* p, double* f)
{
  double x[2] = {r, r};
  double y[2] = {p[1], p[2]};
  double power[2] = {0, 0};

  power_m(2, power, x, y);

  *f = p[0] * (p[2] / (p[2] - p[1]) * power[0] - p[1] / (p[2] - p[1]) * power[1]) + p[3] * r;
}

/****************************************************************
 *
 * constant function
 *
 ****************************************************************/

void const_value(double r, double* p, double* f) { *f = *p; }

/****************************************************************
 *
 * square root function
 *
 * http://dx.doi.org/doi:10.1080/01418618408244210
 *
 ****************************************************************/

void sqrt_value(double r, double* p, double* f) { *f = p[0] * sqrt(r / p[1]); }

/****************************************************************
 *
 * mexp_decay potential
 *
 ****************************************************************/

void mexp_decay_value(double r, double* p, double* f)
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

void strmm_value(double r, double* p, double* f)
{
  double r_0 = r - p[4];

  *f = 2.0 * p[0] * exp(-p[1] / 2.0 * r_0) - p[2] * (1.0 + p[3] * r_0) * exp(-p[3] * r_0);
}

/****************************************************************
 *
 * double morse potential
 *
 * http://dx.doi.org/doi:10.1557/proc-538-535
 *
 ****************************************************************/

void double_morse_value(double r, double* p, double* f)
{
  *f = (p[0] * (exp(-2.0 * p[1] * (r - p[2])) - 2.0 * exp(-p[1] * (r - p[2]))) +
        p[3] * (exp(-2.0 * p[4] * (r - p[5])) - 2.0 * exp(-p[4] * (r - p[5])))) +
       p[6];
}

/****************************************************************
 *
 * double exp potential
 *
 * http://dx.doi.org/doi:10.1557/proc-538-535
 *
 ****************************************************************/

void double_exp_value(double r, double* p, double* f)
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

void poly_5_value(double r, double* p, double* f)
{
  double dr = (r - 1.0) * (r - 1.0);

  *f = p[0] + 0.5 * p[1] * dr + p[2] * (r - 1.0) * dr + p[3] * (dr * dr) +
       p[4] * (dr * dr) * (r - 1.0);
}

/****************************************************************
 *
 * kawamura potential
 *
 * http://dx.doi.org/10.1016/S0925-8388(00)00806-9
 *
 ****************************************************************/

void kawamura_value(double r, double* p, double* f)
{
  double r6 = r * r * r;

  r6 = r6 * r6;

  *f = p[0] * p[1] / r + p[2] * (p[5] + p[6]) * exp((p[3] + p[4] - r) / (p[5] + p[6])) -
       p[7] * p[8] / r6;
}

void kawamura_mix_value(double r, double* p, double* f)
{
  double r6 = r * r * r;

  r6 = r6 * r6;

  *f = p[0] * p[1] / r + p[2] * (p[5] + p[6]) * exp((p[3] + p[4] - r) / (p[5] + p[6])) -
       p[7] * p[8] / r6 +
       p[2] * p[9] * (exp(-2 * p[10] * (r - p[11])) - 2.0 * exp(-p[10] * (r - p[11])));
}

/****************************************************************
 *
 * exp_plus potential
 *
 ****************************************************************/

void exp_plus_value(double r, double* p, double* f) { *f = p[0] * exp(-p[1] * r) + p[2]; }

/****************************************************************
 *
 * mishin potential
 *
 * http://dx.doi.org/doi:10.1016/j.actamat.2005.05.001
 *
 ****************************************************************/

void mishin_value(double r, double* p, double* f)
{
  double z = r - p[3];
  double temp = exp(-p[5] * r);
  double power = 0;

  power_1(&power, &z, &p[4]);

  *f = p[0] * power * temp * (1.0 + p[1] * temp) + p[2];
}

/****************************************************************
 *
 * gen_lj potential, generalized lennard-jones
 *
 * http://dx.doi.org/doi:10.1016/j.actamat.2005.05.001
 *
 ****************************************************************/

void gen_lj_value(double r, double* p, double* f)
{
  double x[2] = {r / p[3], r / p[3]};
  double y[2] = {p[1], p[2]};
  double power[2] = {0, 0};

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

void gljm_value(double r, double* p, double* f)
{
  double x[3] = {r/p[3], r/p[3], r-p[9]};
  double y[3] = {p[1], p[2], p[10]};
  double power[3] = {0, 0, 0};

  power_m(3, power, x, y);

  double temp = exp(-p[11] * power[2]);

  *f = p[0] / (p[2] - p[1]) * (p[2] / power[0] - p[1] / power[1]) + p[4] +
       p[5] * (p[6] * power[2] * temp * (1.0 + p[7] * temp) + p[8]);
}

/****************************************************************
 *
 * bond-stretching function of vashishta potential (f_c)
 *
 * http://dx.doi.org/doi:10.1016/0022-3093(94)90351-4
 *
 ****************************************************************/

void vas_value(double r, double* p, double* f) { *f = exp(p[0] / (r - p[1])); }

/****************************************************************
 *
 * original pair contributions of vashishta potential
 * (V_2 without second "Coulomb"-term)
 *
 * http://dx.doi.org/doi:10.1016/0022-3093(94)90351-4
 *
 ****************************************************************/

void vpair_value(double r, double* p, double* f)
{
  double power = 0;
  double x = r * r;

  x = x * x;

  power_1(&power, &r, &p[1]);

  *f = 14.4 * (p[0] / power - 0.5 * (p[4] * p[3] * p[3] + p[5] * p[2] * p[2] / x) * exp(-r / p[6]));
}

/****************************************************************
 *
 * analytical fits to sheng-aluminum-EAM-potential
 *
 ****************************************************************/

void sheng_phi1_value(double r, double* p, double* f)
{
  double y = r - p[4];
  double z = -p[3] * y * y;

  *f = p[0] * exp(-p[1] * r * r) + p[2] * exp(z);
}

void sheng_phi2_value(double r, double* p, double* f)
{
  double y = r - p[3];
  double z = p[2] * p[2] + y * y;

  *f = p[0] * exp(-p[1] * r * r) + p[2] / z;
}

void sheng_rho_value(double r, double* p, double* f)
{
  double power = 0;

  power_1(&power, &r, &p[1]);

  double x = (p[4] * p[4]) / (r * r);
  x = x * x * x;

  if (r > 1.45)
    *f = 4.0 * p[3] * x * (x - 1.0);
  else
    *f = p[0] * power + p[2];
}

void sheng_F_value(double r, double* p, double* f)
{
  double power = 0;

  power_1(&power, &r, &p[1]);

  *f = p[0] * power + p[2] * r + p[3];
}

#if defined(STIWEB)

/****************************************************************
 *
 * Stillinger-Weber pair potential
 *
 ****************************************************************/

void stiweb_2_value(double r, double* p, double* f)
{
  double x[2] = {r, r};
  double y[2] = {-p[2], -p[3]};
  double power[2] = {0, 0};

  power_m(2, power, x, y);

  *f = (p[0] * power[0] - p[1] * power[1]) * exp(p[4] / (r - p[5]));
}

/****************************************************************
 *
 * Stillinger-Weber exp functions for threebody potential
 *
 ****************************************************************/

void stiweb_3_value(double r, double* p, double* f) { *f = exp(p[0] / (r - p[1])); }

/****************************************************************
 *
 * pseudo Stillinger-Weber potential function to store lamda values
 *
 ****************************************************************/

void lambda_value(double r, double* p, double* f) { *f = 0.0 * r* p[0]; }

#endif // STIWEB

#if defined(TERSOFF)
#if !defined(TERSOFFMOD)

/****************************************************************
 *
 * pseudo Tersoff potential function to store potential parameters
 *
 ****************************************************************/

void tersoff_pot_value(double r, double* p, double* f) { *f = 0.0 * r* p[0]; }

/****************************************************************
 *
 * pseudo Tersoff potential function to store mixing potential values
 *
 ****************************************************************/

void tersoff_mix_value(double r, double* p, double* f) { *f = 0.0 * r* p[0]; }

#else

/****************************************************************
 *
 * pseudo modified Tersoff potential function to store potential parameters
 *
 ****************************************************************/

void tersoff_mod_pot_value(double r, double* p, double* f) { *f = 0.0 * r * p[0]; }

#endif // !TERSOFFMOD
#endif // TERSOFF

/****************************************************************
 *
 * template for new potential function called mypotential
 *
 * for further information plase have a look at the online
 * documentation available at
 * 	http://potfit.sourceforge.net/
 *
 ****************************************************************/

/* template for new function */

/****************************************************************
 *
 * newpot potential
 *
 ****************************************************************/

void newpot_value(double r, double* p, double* f) { *f = r * p[0] + p[1]; }

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

  return val / (1.0 + val);
}

/****************************************************************
 *
 * check analytic parameters for special conditions
 *
 ****************************************************************/

int apot_check_params(double* params)
{
  int i, j = 2, k;
#ifdef TERSOFF
  double temp;
#endif /* TERSOFF */

  for (i = 0; i < g_pot.apot_table.number; i++)
  {
    /* last parameter of eopp potential is 2 pi periodic */
    if (strcmp(g_pot.apot_table.names[i], "eopp") == 0)
    {
      k = j + 5;
      if (params[k] > 2 * M_PI)
        do
        {
          params[k] -= 2 * M_PI;
        } while (params[k] > 2 * M_PI);
      if (params[k] < 0)
        do
        {
          params[k] += 2 * M_PI;
        } while (params[k] < 0);
    }

    /* the third parameter of csw2 potential is 2 pi periodic */
    if (strcmp(g_pot.apot_table.names[i], "csw2") == 0)
    {
      k = j + 2;
      if (params[k] > 2 * M_PI)
        do
        {
          params[k] -= 2 * M_PI;
        } while (params[k] > 2 * M_PI);
      if (params[k] < 0)
        do
        {
          params[k] += 2 * M_PI;
        } while (params[k] < 0);
    }
#ifdef TERSOFF
    /* the parameter S has to be greater than the parameter R */
    /* switch them if this is not the case */
    if (strcmp(g_pot.apot_table.names[i], "tersoff_pot") == 0)
    {
      k = j + 9;
      if (params[k] < params[k + 1])
      {
        temp = params[k];
        params[k] = params[k + 1];
        params[k + 1] = temp;
      }
    }
#endif /* TERSOFF */

    /* jump to next potential */
    j += 2 + apot_parameters(g_pot.apot_table.names[i]) + g_pot.smooth_pot[i];
  }

  return 0;
}

/****************************************************************
 *
 * punish analytic potential for bad habits
 *
 ****************************************************************/

double apot_punish(double* params, double* forces)
{
  int i, j;
  double x, tmpsum = 0.0, min, max;

  /* loop over individual parameters */
  for (i = 0; i < g_calc.ndim; i++)
  {
    min = g_pot.apot_table.pmin[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]];
    max = g_pot.apot_table.pmax[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]];
    /* punishment for out of bounds */
    if (x = params[g_todo.idx[i]] - min, x < 0)
    {
      tmpsum += APOT_PUNISH * x * x;
      forces[g_calc.punish_par_p + i] = APOT_PUNISH * x * x;
    }
    else if (x = params[g_todo.idx[i]] - max, x > 0)
    {
      tmpsum += APOT_PUNISH * x * x;
      forces[g_calc.punish_par_p + i] = APOT_PUNISH * x * x;
    }
  }

  j = 2;
  /* loop over potentials */
  for (i = 0; i < g_pot.apot_table.number; i++)
  {
    /* punish eta_1 < eta_2 for eopp function */
    if (strcmp(g_pot.apot_table.names[i], "eopp") == 0)
    {
      x = params[j + 1] - params[j + 3];
      if (x < 0)
      {
        forces[g_calc.punish_pot_p + i] = g_param.apot_punish_value * (1 + x) * (1 + x);
        tmpsum += g_param.apot_punish_value * (1 + x) * (1 + x);
      }
    }
#if defined EAM || defined ADP || defined MEAM
    /* punish m=n for universal embedding function */
    if (strcmp(g_pot.apot_table.names[i], "universal") == 0)
    {
      x = params[j + 2] - params[j + 1];
      if (fabs(x) < 1e-6)
      {
        forces[g_calc.punish_pot_p + i] = g_param.apot_punish_value / (x * x);
        tmpsum += g_param.apot_punish_value / (x * x);
      }
    }
#endif /* EAM */

    /* jump to next potential */
    j += 2 + g_pot.apot_table.n_par[i];
  }

  return tmpsum;
}

/****************************************************************
 *
 * calculate gradient for analytic potential
 *
 ****************************************************************/

double apot_grad(double r, double* p, void (*function)(double, double*, double*))
{
  double a, b, h = 0.0001;

  function(r + h, p, &a);
  function(r - h, p, &b);

  return (a - b) / (2.0 * h);
}

#if defined(COULOMB)

/****************************************************************
 *
 * ms potential + first derivative
 *
 ****************************************************************/

void ms_init(double r, double* pot, double* grad, double* p)
{
  double x[4];

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

void buck_init(double r, double* pot, double* grad, double* p)
{
  double x[3];

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

void ms_shift(double r, double* p, double* f)
{
  double pot, grad, pot_cut, grad_cut;

  ms_init(r, &pot, &grad, p);
  ms_init(g_todo.dp_cut, &pot_cut, &grad_cut, p);

  *f = pot - pot_cut - (r - g_todo.dp_cut) * grad_cut;
}

/****************************************************************
 *
 * shifted buckingham potential
 *
 ****************************************************************/

void buck_shift(double r, double* p, double* f)
{
  double pot, grad, pot_cut, grad_cut;

  buck_init(r, &pot, &grad, p);
  buck_init(g_todo.dp_cut, &pot_cut, &grad_cut, p);

  *f = pot - pot_cut - r*(r - g_todo.dp_cut) * grad_cut;
}

/****************************************************************
 *
 * tail of electrostatic potential and first two derivatives
 *
 ****************************************************************/

void elstat_value(double r, double dp_kappa, double* ftail, double* gtail, double* ggtail)
{
  double x[4];

  x[0] = r * r;
  x[1] = dp_kappa * dp_kappa;
  x[2] = 2 * g_todo.dp_eps * dp_kappa / sqrt(M_PI);
  x[3] = exp(-x[0] * x[1]);

  *ftail = g_todo.dp_eps* erfc(dp_kappa * r) / r;
  *gtail = -(*ftail + x[2] * x[3]) / x[0];
  *ggtail = (2 * x[1] * x[2] * x[3] - *gtail * 3) / x[0];
}

/****************************************************************
 *
 * shifted tail of coloumb potential
 *
 ****************************************************************/

void elstat_shift(double r, double dp_kappa, double* fnval_tail, double* grad_tail,
                  double* ggrad_tail)
{
  double ftail, gtail, ggtail, ftail_cut, gtail_cut, ggtail_cut;
  double x[3];

  x[0] = r * r;
  x[1] = g_todo.dp_cut * g_todo.dp_cut;
  x[2] = x[0] - x[1];

  elstat_value(r, dp_kappa, &ftail, &gtail, &ggtail);
  elstat_value(g_todo.dp_cut, dp_kappa, &ftail_cut, &gtail_cut, &ggtail_cut);

  *fnval_tail = ftail - ftail_cut - x[2] * gtail_cut / 2;
  *grad_tail = gtail - gtail_cut;
  *ggrad_tail = 0.0;
#if defined(DIPOLE)
  *fnval_tail -= x[2] * x[2] * ggtail_cut / 8;
  *grad_tail -= x[2] * ggtail_cut / 2;
  *ggrad_tail = ggtail - ggtail_cut;
#endif /* DIPOLE */
}

/****************************************************************
 *
 *  calculate tail of coulomb-potential and its first derivative
 *
 ****************************************************************/

void init_tails(double dp_kappa)
{
  for (int i = 0; i < g_config.natoms; i++)
  {
    for (int j = 0; j < g_config.atoms[i].num_neigh; j++)
    {
      elstat_shift(
          g_config.atoms[i].neigh[j].r, dp_kappa, &g_config.atoms[i].neigh[j].fnval_el,
          &g_config.atoms[i].neigh[j].grad_el, &g_config.atoms[i].neigh[j].ggrad_el);
    }
  }
}

#endif // COULOMB

#if defined(DIPOLE)

/****************************************************************
 *
 * short-range part of dipole moments
 *
 ****************************************************************/

double shortrange_value(double r, double a, double b, double c)
{
  double x[5];

  x[0] = b * r;
  x[1] = x[0] * x[0];
  x[2] = x[1] * x[0];
  x[3] = x[1] * x[1];
  x[4] = 1 + x[0] + x[1] / 2 + x[2] / 6 + x[3] / 24;

  return a * c * x[4] * exp(-x[0]) / g_todo.dp_eps;
}

/****************************************************************
 *
 * tail of additional short-range contribution to energy and forces
 *
 ****************************************************************/

void shortrange_term(double r, double b, double c, double* srval_tail,
                     double* srgrad_tail)
{
  double x[6];

  x[0] = b * r;
  x[1] = x[0] * x[0];
  x[2] = x[1] * x[0];
  x[3] = x[1] * x[1];
  x[4] = 1 + x[0] + x[1] / 2 + x[2] / 6 + x[3] / 24;
  x[5] = exp(-x[0]);

  *srval_tail = c* x[4] * x[5] / g_todo.dp_eps;
  *srgrad_tail = -c* b* x[3] * x[5] / (24 * g_todo.dp_eps * r);
}

#endif // DIPOLE

/****************************************************************
 *
 *  debug function to print the potentials
 *
 ****************************************************************/

#if defined(DEBUG)

void debug_apot()
{
  fflush(stdout);
  fprintf(stderr, "\n##############################################\n");
  fprintf(stderr, "###########      DEBUG OUTPUT      ###########\n");
  fprintf(stderr, "##############################################\n");
  fprintf(stderr, "\nThere are %d potentials with a total of %d parameters.\n",
          g_pot.apot_table.number, g_pot.apot_table.total_par);
  for (int i = 0; i < g_pot.apot_table.number; i++)
  {
    fprintf(stderr, "\npotential #%d (type=%s, smooth=%d)\n", i + 1,
            g_pot.apot_table.names[i], g_pot.smooth_pot[i]);
    fprintf(stderr, "begin=%f end=%f\n", g_pot.apot_table.begin[i],
            g_pot.apot_table.end[i]);
    for (int j = 0; j < g_pot.apot_table.n_par[i]; j++)
    {
      fprintf(stderr, "parameter %d: name=%s value=%f min=%f max=%f\n", j + 1,
              g_pot.apot_table.param_name[i][j], g_pot.apot_table.values[i][j],
              g_pot.apot_table.pmin[i][j], g_pot.apot_table.pmax[i][j]);
    }
  }

#if defined(PAIR)
  if (!g_param.enable_cp)
  {
    fprintf(stderr, "\nchemical potentials are DISABLED!\n");
  }
  else
  {
    fprintf(stderr, "\nchemical potentials:\n");
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(stderr, "cp_%d=%f min=%f max=%f\n", i, g_pot.apot_table.chempot[i],
              g_pot.apot_table.pmin[g_pot.apot_table.number][i],
              g_pot.apot_table.pmax[g_pot.apot_table.number][i]);
    if (g_param.compnodes > 0 && g_param.ntypes == 2)
    {
      fprintf(stderr, "composition nodes:\n");
      for (int j = 0; j < g_param.compnodes; j++)
        fprintf(stderr, "composition=%f value=%f min=%f max=%f\n",
                g_pot.compnodelist[j], g_pot.apot_table.chempot[g_param.ntypes + j],
                g_pot.apot_table.pmin[g_pot.apot_table.number][g_param.ntypes + j],
                g_pot.apot_table.pmax[g_pot.apot_table.number][g_param.ntypes + j]);
    }
  }
#endif // PAIR

  exit(EXIT_FAILURE);
}

#endif // DEBUG
