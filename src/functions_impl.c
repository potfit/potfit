/****************************************************************
 *
 * functions_impl.c: Routines and function calls used for analytic potentials
 *
 ****************************************************************
 *
 * Copyright 2002-2017 - the potfit development team
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

#include "potfit.h"

#include "utils.h"

/****************************************************************

  actual functions representing the analytic potentials

****************************************************************/

/****************************************************************
  lennard-jones potential
    http://dx.doi.org/doi:10.1098/rspa.1924.0082
****************************************************************/

void lj_value(const double r, const double* p, double* f)
{
  double x = (p[1] * p[1]) / (r * r);
  x = x * x * x;

  *f = 4.0 * p[0] * x * (x - 1.0);
}

/****************************************************************
  empirical oscillating pair potential (eopp)
    http://arxiv.org/abs/0802.2926v2
****************************************************************/

void eopp_value(const double r, const double* p, double* f)
{
  double x[2] = {r, r};
  double y[2] = {p[1], p[3]};
  double power[2] = {0, 0};

  power_m(2, power, x, y);

  *f = p[0] / power[0] + (p[2] / power[1]) * cos(p[4] * r + p[5]);
}

/****************************************************************
  morse potential
    http://dx.doi.org/doi:10.1103/PhysRev.34.57
****************************************************************/

void morse_value(const double r, const double* p, double* f)
{
  *f = p[0] * (exp(-2.0 * p[1] * (r - p[2])) - 2.0 * exp(-p[1] * (r - p[2])));
}

/****************************************************************
  morse-stretch potential (without derivative!)
    http://dx.doi.org/doi:10.1063/1.1513312
****************************************************************/

#if defined(COULOMB)

void ms_init(const double r, const double* p, double* pot, double* grad)
{
  double x[4];

  x[0] = 1 - r / p[2];
  x[1] = exp(p[1] * x[0]);
  x[2] = exp(p[1] * x[0] / 2);
  x[3] = p[0] * p[1] / (r * p[2]);

  *pot = p[0] * (x[1] - 2 * x[2]);
  *grad = x[3] * (-x[1] + x[2]);
}

void ms_value(const double r, const double* p, double* f)
{
  double pot = 0.0;
  double grad = 0.0;
  double pot_cut = 0.0;
  double grad_cut = 0.0;

  ms_init(r, p, &pot, &grad);
  ms_init(g_config.dp_cut, p, &pot_cut, &grad_cut);

  *f = pot - pot_cut - (r - g_config.dp_cut) * grad_cut;
}

#else

void ms_value(const double r, const double* p, double* f)
{
  double x = 1.0 - r / p[2];

  *f = p[0] * (exp(p[1] * x) - 2.0 * exp((p[1] * x) / 2.0));
}

#endif  // COULOMB

/****************************************************************
  buckingham potential (without derivative!) - slightly modified
    http://dx.doi.org/doi:10.1098/rspa.1977.0049
****************************************************************/

#if defined(COULOMB)

void buck_init(const double r, const double* p, double* pot, double* grad)
{
  double x[3];

  x[0] = dsquare(p[1]) / dsquare(r);
  x[1] = p[2] * x[0] * x[0] * x[0];
  x[2] = p[0] * exp(-r / p[1]);

  *pot = x[2] - x[1];
  *grad = -x[2] / p[1] + 6 * x[1] / r;
}

void buck_value(const double r, const double* p, double* f)
{
  double pot = 0.0;
  double grad = 0.0;
  double pot_cut = 0.0;
  double grad_cut = 0.0;

  buck_init(r, p, &pot, &grad);
  buck_init(g_config.dp_cut, p, &pot_cut, &grad_cut);

  *f = pot - pot_cut - r * (r - g_config.dp_cut) * grad_cut;
}

#else

void buck_value(const double r, const double* p, double* f)
{
  double x = (p[1] * p[1]) / (r * r);
  double y = x * x * x;

  *f = p[0] * exp(-r / p[1]) - p[2] * y;
}

#endif  // COULOMB

/****************************************************************
  softshell potential
    unknown reference
*****************************************************************/

void softshell_value(const double r, const double* p, double* f)
{
  double x = p[0] / r;

  power_1(f, &x, &p[1]);
}

/****************************************************************
  eopp_exp potential
    http://arxiv.org/abs/0802.2926v2
****************************************************************/

void eopp_exp_value(const double r, const double* p, double* f)
{
  double power = 0.0;

  power_1(&power, &r, &p[3]);

  *f = p[0] * exp(-p[1] * r) + (p[2] / power) * cos(p[4] * r + p[5]);
}

/****************************************************************
  meopp potential
    http://arxiv.org/abs/0802.2926v2
****************************************************************/

void meopp_value(const double r, const double* p, double* f)
{
  double x[2] = {r - p[6], r};
  double y[2] = {p[1], p[3]};
  double power[2] = {0, 0};

  power_m(2, power, x, y);

  *f = p[0] / power[0] + (p[2] / power[1]) * cos(p[4] * r + p[5]);
}

/****************************************************************
  power potential
    unknown reference
****************************************************************/

void power_value(const double r, const double* p, double* f)
{
  double power = 0;

  power_1(&power, &r, &p[1]);

  *f = p[0] * power;
}

/****************************************************************
  power_decay potential
    unknown reference
****************************************************************/

void power_decay_value(const double r, const double* p, double* f)
{
  double power = 0;

  power_1(&power, &r, &p[1]);

  *f = p[0] / power;
}

/****************************************************************
  exp_decay potential
    unknown reference
****************************************************************/

void exp_decay_value(const double r, const double* p, double* f)
{
  *f = p[0] * exp(-p[1] * r);
}

/****************************************************************
  bjs potential
    http://dx.doi.org/doi:10.1103/PhysRevB.37.6632
****************************************************************/

void bjs_value(const double r, const double* p, double* f)
{
  if (r == 0.0)
    *f = 0.0;
  else {
    double power = 0.0;

    power_1(&power, &r, &p[1]);

    *f = p[0] * (1.0 - p[1] * log(r)) * power + p[2] * r;
  }
}

/****************************************************************
  parabola potential
    unknown reference
****************************************************************/

void parabola_value(const double r, const double* p, double* f)
{
  *f = (r * r) * p[0] + r * p[1] + p[2];
}

/****************************************************************
  harmonic potential
    unknown reference
****************************************************************/

void harmonic_value(const double r, const double* p, double* f)
{
  *f = p[0] * (r - p[1]) * (r - p[1]);
}

/****************************************************************
  angular harmonic potential
    unknown reference
****************************************************************/

void angharmonic_value(const double r, const double* p, double* f)
{
  double theta = acos(r);

  *f = p[0] * (theta - p[1]) * (theta - p[1]);
}

/****************************************************************
  chantasiriwan (csw) and milstein potential
    http://dx.doi.org/doi:10.1103/PhysRevB.53.14080
****************************************************************/

void csw_value(const double r, const double* p, double* f)
{
  static double power;

  power_1(&power, &r, &p[3]);

  *f = (1.0 + p[0] * cos(p[2] * r) + p[1] * sin(p[2] * r)) / power;
}

/****************************************************************
  chantasiriwan (csw) and milstein potential - slightly modified
    http://dx.doi.org/doi:10.1103/PhysRevB.53.14080
****************************************************************/

void csw2_value(const double r, const double* p, double* f)
{
  static double power;

  power_1(&power, &r, &p[3]);

  *f = (1.0 + p[0] * cos(p[1] * r + p[2])) / power;
}

/****************************************************************
  universal embedding function
    http://dx.doi.org/doi:10.1557/jmr.1989.1195
****************************************************************/

void universal_value(const double r, const double* p, double* f)
{
  double x[2] = {r, r};
  double y[2] = {p[1], p[2]};
  double power[2] = {0, 0};

  power_m(2, power, x, y);

  *f = p[0] *
           (p[2] / (p[2] - p[1]) * power[0] - p[1] / (p[2] - p[1]) * power[1]) +
       p[3] * r;
}

/****************************************************************
  constant function
    unknown reference
****************************************************************/

void const_value(const double r, const double* p, double* f) { *f = *p; }
/****************************************************************
  square root function
    http://dx.doi.org/doi:10.1080/01418618408244210
****************************************************************/

void sqrt_value(const double r, const double* p, double* f)
{
  *f = p[0] * sqrt(r / p[1]);
}

/****************************************************************
  mexp_decay potential
    unknown reference
****************************************************************/

void mexp_decay_value(const double r, const double* p, double* f)
{
  *f = p[0] * exp(-p[1] * (r - p[2]));
}

/****************************************************************
  streitz-mintmire (strmm) potential
    http://dx.doi.org/doi:10.1103/PhysRevB.50.11996
****************************************************************/

void strmm_value(const double r, const double* p, double* f)
{
  double r_0 = r - p[4];

  *f = 2.0 * p[0] * exp(-p[1] / 2.0 * r_0) -
       p[2] * (1.0 + p[3] * r_0) * exp(-p[3] * r_0);
}

/****************************************************************
  double morse potential
    http://dx.doi.org/doi:10.1557/proc-538-535
****************************************************************/

void double_morse_value(const double r, const double* p, double* f)
{
  *f =
      (p[0] * (exp(-2.0 * p[1] * (r - p[2])) - 2.0 * exp(-p[1] * (r - p[2]))) +
       p[3] * (exp(-2.0 * p[4] * (r - p[5])) - 2.0 * exp(-p[4] * (r - p[5])))) +
      p[6];
}

/****************************************************************
  double exp potential
    http://dx.doi.org/doi:10.1557/proc-538-535
****************************************************************/

void double_exp_value(const double r, const double* p, double* f)
{
  *f = (p[0] * exp(-p[1] * dsquare(r - p[2])) + exp(-p[3] * (r - p[4])));
}

/****************************************************************
  poly 5 potential
    http://dx.doi.org/doi:10.1557/proc-538-535
****************************************************************/

void poly_5_value(const double r, const double* p, double* f)
{
  double dr = (r - 1.0) * (r - 1.0);

  *f = p[0] + 0.5 * p[1] * dr + p[2] * (r - 1.0) * dr + p[3] * (dr * dr) +
       p[4] * (dr * dr) * (r - 1.0);
}

/****************************************************************
  kawamura potential
    http://dx.doi.org/10.1016/S0925-8388(00)00806-9
****************************************************************/

void kawamura_value(const double r, const double* p, double* f)
{
  double r6 = r * r * r;

  r6 = r6 * r6;

  *f = p[0] * p[1] / r +
       p[2] * (p[5] + p[6]) * exp((p[3] + p[4] - r) / (p[5] + p[6])) -
       p[7] * p[8] / r6;
}

/****************************************************************
  kawamura mixing potential
    http://dx.doi.org/10.1016/S0925-8388(00)00806-9
****************************************************************/

void kawamura_mix_value(const double r, const double* p, double* f)
{
  double r6 = r * r * r;

  r6 = r6 * r6;

  *f = p[0] * p[1] / r +
       p[2] * (p[5] + p[6]) * exp((p[3] + p[4] - r) / (p[5] + p[6])) -
       p[7] * p[8] / r6 +
       p[2] * p[9] *
           (exp(-2 * p[10] * (r - p[11])) - 2.0 * exp(-p[10] * (r - p[11])));
}

/****************************************************************
  exp_plus potential
    unknown reference
****************************************************************/

void exp_plus_value(const double r, const double* p, double* f)
{
  *f = p[0] * exp(-p[1] * r) + p[2];
}

/****************************************************************
  mishin potential
    http://dx.doi.org/doi:10.1016/j.actamat.2005.05.001
****************************************************************/

void mishin_value(const double r, const double* p, double* f)
{
  double z = r - p[3];
  double temp = exp(-p[5] * r);
  double power = 0;

  power_1(&power, &z, &p[4]);

  *f = p[0] * power * temp * (1.0 + p[1] * temp) + p[2];
}

/****************************************************************
  gen_lj potential, generalized lennard-jones
    http://dx.doi.org/doi:10.1016/j.actamat.2005.05.001
****************************************************************/

void gen_lj_value(const double r, const double* p, double* f)
{
  double x[2] = {r / p[3], r / p[3]};
  double y[2] = {p[1], p[2]};
  double power[2] = {0, 0};

  power_m(2, power, x, y);

  *f = p[0] / (p[2] - p[1]) * (p[2] / power[0] - p[1] / power[1]) + p[4];
}

/****************************************************************
  gljm potential, generalized lennard-jones + mishin potential
    http://dx.doi.org/doi:10.1016/j.actamat.2005.05.001
****************************************************************/

void gljm_value(const double r, const double* p, double* f)
{
  double x[3] = {r / p[3], r / p[3], r - p[9]};
  double y[3] = {p[1], p[2], p[10]};
  double power[3] = {0.0, 0.0, 0.0};

  power_m(3, power, x, y);

  double temp = exp(-p[11] * power[2]);

  *f = p[0] / (p[2] - p[1]) * (p[2] / power[0] - p[1] / power[1]) + p[4] +
       p[5] * (p[6] * power[2] * temp * (1.0 + p[7] * temp) + p[8]);
}

/****************************************************************
  bond-stretching function of vashishta potential (f_c)
    http://dx.doi.org/doi:10.1016/0022-3093(94)90351-4
****************************************************************/

void vas_value(const double r, const double* p, double* f)
{
  *f = exp(p[0] / (r - p[1]));
}

/****************************************************************
  original pair contributions of vashishta potential
  (V_2 without second "Coulomb"-term)
    http://dx.doi.org/doi:10.1016/0022-3093(94)90351-4
****************************************************************/

void vpair_value(const double r, const double* p, double* f)
{
  double power = 0.0;
  double x = r * r;

  x *= x;

  power_1(&power, &r, &p[1]);

  *f = 14.4 *
       (p[0] / power -
        0.5 * (p[4] * p[3] * p[3] + p[5] * p[2] * p[2] / x) * exp(-r / p[6]));
}

/****************************************************************
  analytical fits to sheng-aluminum-EAM-potential
    unknown reference
****************************************************************/

void sheng_phi1_value(const double r, const double* p, double* f)
{
  double y = r - p[4];
  double z = -p[3] * y * y;

  *f = p[0] * exp(-p[1] * r * r) + p[2] * exp(z);
}

/****************************************************************
  analytical fits to sheng-aluminum-EAM-potential
    unknown reference
****************************************************************/

void sheng_phi2_value(const double r, const double* p, double* f)
{
  double y = r - p[3];
  double z = p[2] * p[2] + y * y;

  *f = p[0] * exp(-p[1] * r * r) + p[2] / z;
}

/****************************************************************
  analytical fits to sheng-aluminum-EAM-potential
    unknown reference
****************************************************************/

void sheng_rho_value(const double r, const double* p, double* f)
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

/****************************************************************
  analytical fits to sheng-aluminum-EAM-potential
    unknown reference
****************************************************************/

void sheng_F_value(const double r, const double* p, double* f)
{
  double power = 0;

  power_1(&power, &r, &p[1]);

  *f = p[0] * power + p[2] * r + p[3];
}

#if defined(STIWEB)

/****************************************************************
  Stillinger-Weber pair potential
****************************************************************/

void stiweb_2_value(const double r, const double* p, double* f)
{
  double x[2] = {r, r};
  double y[2] = {-p[2], -p[3]};
  double power[2] = {0, 0};

  power_m(2, power, x, y);

  *f = (p[0] * power[0] - p[1] * power[1]) * exp(p[4] / (r - p[5]));
}

/****************************************************************
  Stillinger-Weber exp functions for threebody potential
 ****************************************************************/

void stiweb_3_value(const double r, const double* p, double* f)
{
  *f = exp(p[0] / (r - p[1]));
}

/****************************************************************
  pseudo Stillinger-Weber potential function to store lamda values
****************************************************************/

void lambda_value(const double r, const double* p, double* f)
{
  (void)r;
  (void)p;
  (void)f;
}

#endif  // STIWEB

#if defined(TERSOFF)
#if !defined(TERSOFFMOD)

/****************************************************************
 *
 * pseudo Tersoff potential function to store potential parameters
 *
 ****************************************************************/

void tersoff_pot_value(const double r, const double* p, double* f)
{
  (void)r;
  (void)p;
  (void)f;
}

/****************************************************************
 *
 * pseudo Tersoff potential function to store mixing potential values
 *
 ****************************************************************/

void tersoff_mix_value(const double r, const double* p, double* f)
{
  (void)r;
  (void)p;
  (void)f;
}

#else

/****************************************************************
 *
 * pseudo modified Tersoff potential function to store potential parameters
 *
 ****************************************************************/

void tersoff_mod_pot_value(const double r, const double* p, double* f)
{
  (void)r;
  (void)p;
  (void)f;
}

#endif  // !TERSOFFMOD
#endif  // TERSOFF

/****************************************************************

  end of analytic potentials

****************************************************************/
