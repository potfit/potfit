/****************************************************************
 *
 * functions.c: Routines and function calls used for analytic potentials
 *
 ****************************************************************
 *
 * Copyright 2008-2010 Daniel Schopf, Philipp Beck
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

#include "potfit.h"
#ifndef ACML
#include <mkl_vml.h>
#else
#include <acml_mv.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

/****************************************************************
 *
 * return the number of parameters for a specific analytic potential
 *
 ****************************************************************/

int apot_parameters(char *name)
{
  if (strcmp(name, "lj") == 0) {
    return 2;
  } else if (strcmp(name, "eopp") == 0) {
    return 6;
  } else if (strcmp(name, "morse") == 0) {
    return 3;
  } else if (strcmp(name, "ms") == 0) {
    return 3;
  } else if (strcmp(name, "buck") == 0) {
    return 3;
  } else if (strcmp(name, "softshell") == 0) {
    return 2;
  } else if (strcmp(name, "eopp_exp") == 0) {
    return 6;
  } else if (strcmp(name, "meopp") == 0) {
    return 7;
  } else if (strcmp(name, "power_decay") == 0) {
    return 2;
  } else if (strcmp(name, "exp_decay") == 0) {
    return 2;
  } else if (strcmp(name, "pohlong") == 0) {
    return 3;
  } else if (strcmp(name, "parabola") == 0) {
    return 3;
  } else if (strcmp(name, "csw") == 0) {
    return 4;
  } else if (strcmp(name, "universal") == 0) {
    return 4;
  } else if (strcmp(name, "const") == 0) {
    return 1;
  } else if (strcmp(name, "sqrt") == 0) {
    return 2;
  } else if (strcmp(name, "mexp_decay") == 0) {
    return 3;
  } else if (strcmp(name, "strmm") == 0) {
    return 5;
  } else if (strcmp(name, "double_morse") == 0) {
    return 7;
  } else if (strcmp(name, "double_exp") == 0) {
    return 5;
  } else if (strcmp(name, "poly_5") == 0) {
    return 5;
  } else if (strcmp(name, "cbb") == 0) {
    return 12;
  } else if (strcmp(name, "exp_plus") == 0) {
    return 3;
  } else if (strcmp(name, "mishin") == 0) {
    return 6;
  } else if (strcmp(name, "gen_lj") == 0) {
    return 5;
  } else if (strcmp(name, "gljm") == 0) {
    return 12;
  } else if (strcmp(name, "vas") == 0) {
    return 2;
  } else if (strcmp(name, "vpair") == 0) {
    return 7;
  }

  /* template for new potential function called newpot */

  else if (strcmp(name, "newpot") == 0) {
    return 2;
  }

  /* end of template */

  return -1;
}

/****************************************************************
 *
 * assign function pointers to corresponding functions
 *
 ****************************************************************/

int apot_assign_functions(apot_table_t *apt)
{
  int   i;

  for (i = 0; i < apt->number; i++) {
    if (strcmp(apt->names[i], "lj") == 0) {
      apt->fvalue[i] = &lj_value;
    } else if (strcmp(apt->names[i], "eopp") == 0) {
      apt->fvalue[i] = &eopp_value;
    } else if (strcmp(apt->names[i], "morse") == 0) {
      apt->fvalue[i] = &morse_value;
    } else if (strcmp(apt->names[i], "ms") == 0) {
#ifdef COULOMB
      apt->fvalue[i] = &ms_shift;
#else
      apt->fvalue[i] = &ms_value;
#endif
    } else if (strcmp(apt->names[i], "buck") == 0) {
#ifdef COULOMB
      apt->fvalue[i] = &buck_shift;
#else
      apt->fvalue[i] = &buck_value;
#endif
    } else if (strcmp(apt->names[i], "softshell") == 0) {
      apt->fvalue[i] = &softshell_value;
    } else if (strcmp(apt->names[i], "eopp_exp") == 0) {
      apt->fvalue[i] = &eopp_exp_value;
    } else if (strcmp(apt->names[i], "meopp") == 0) {
      apt->fvalue[i] = &meopp_value;
    } else if (strcmp(apt->names[i], "power_decay") == 0) {
      apt->fvalue[i] = &power_decay_value;
    } else if (strcmp(apt->names[i], "exp_decay") == 0) {
      apt->fvalue[i] = &exp_decay_value;
    } else if (strcmp(apt->names[i], "pohlong") == 0) {
      apt->fvalue[i] = &pohlong_value;
    } else if (strcmp(apt->names[i], "parabola") == 0) {
      apt->fvalue[i] = &parabola_value;
    } else if (strcmp(apt->names[i], "csw") == 0) {
      apt->fvalue[i] = &csw_value;
    } else if (strcmp(apt->names[i], "universal") == 0) {
      apt->fvalue[i] = &universal_value;
    } else if (strcmp(apt->names[i], "const") == 0) {
      apt->fvalue[i] = &const_value;
    } else if (strcmp(apt->names[i], "sqrt") == 0) {
      apt->fvalue[i] = &sqrt_value;
    } else if (strcmp(apt->names[i], "mexp_decay") == 0) {
      apt->fvalue[i] = &mexp_decay_value;
    } else if (strcmp(apt->names[i], "strmm") == 0) {
      apt->fvalue[i] = &strmm_value;
    } else if (strcmp(apt->names[i], "double_morse") == 0) {
      apt->fvalue[i] = &double_morse_value;
    } else if (strcmp(apt->names[i], "double_exp") == 0) {
      apt->fvalue[i] = &double_exp_value;
    } else if (strcmp(apt->names[i], "poly_5") == 0) {
      apt->fvalue[i] = &poly_5_value;
    } else if (strcmp(apt->names[i], "cbb") == 0) {
      apt->fvalue[i] = &cbb_value;
    } else if (strcmp(apt->names[i], "exp_plus") == 0) {
      apt->fvalue[i] = &exp_plus_value;
    } else if (strcmp(apt->names[i], "mishin") == 0) {
      apt->fvalue[i] = &mishin_value;
    } else if (strcmp(apt->names[i], "gen_lj") == 0) {
      apt->fvalue[i] = &gen_lj_value;
    } else if (strcmp(apt->names[i], "gljm") == 0) {
      apt->fvalue[i] = &gljm_value;
    } else if (strcmp(apt->names[i], "vas") == 0) {
      apt->fvalue[i] = &vas_value;
    } else if (strcmp(apt->names[i], "vpair") == 0) {
      apt->fvalue[i] = &vpair_value;
    }

/* template for new potential function called newpot */

    else if (strcmp(apt->names[i], "newpot") == 0) {
      apt->fvalue[i] = &newpot_value;
    }

/* end of template */

    else
      return -1;
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
 ****************************************************************/

void lj_value(real r, real *p, real *f)
{
  real  sig_d_rad6, sig_d_rad12;

  sig_d_rad6 = SQR(p[1]) / SQR(r);
  sig_d_rad6 = sig_d_rad6 * sig_d_rad6 * sig_d_rad6;
  sig_d_rad12 = SQR(sig_d_rad6);

  *f = 4 * p[0] * (sig_d_rad12 - sig_d_rad6);
}

/****************************************************************
 *
 * empirical oscillating pair potential
 *
 ****************************************************************/

void eopp_value(real r, real *p, real *f)
{
  static real x[2], y[2], power[2];

  x[0] = r;
  x[1] = r;
  y[0] = p[1];
  y[1] = p[3];
#ifndef ACML
  vdPow(2, x, y, power);
#else
  power[0] = fastpow(x[0], y[0]);
  power[1] = fastpow(x[1], y[1]);
#endif

  *f = p[0] / power[0] + (p[2] / power[1]) * cos(p[4] * r + p[5]);
}

/****************************************************************
 *
 * morse potential
 *
 ****************************************************************/

void morse_value(real r, real *p, real *f)
{
  *f = p[0] * (exp(-2 * p[1] * (r - p[2])) - 2 * exp(-p[1] * (r - p[2])));
}

/****************************************************************
 *
 * morse-stretch potential (without derivative!)
 *
 ****************************************************************/

void ms_value(real r, real *p, real *f)
{
  static real x;

  x = 1 - r / p[2];

  *f = p[0] * (exp(p[1] * x) - 2 * exp((p[1] * x) / 2));
}

/****************************************************************
 *
 * buckingham potential (without derivative!)
 *
 ****************************************************************/

void buck_value(real r, real *p, real *f)
{
  static real x, y;
  
  x = SQR(p[1]) / SQR(r);
  y = x * x * x;
  
  *f = p[0] * exp(-r / p[1]) - p[2] * y;
}

/****************************************************************
 *
 * softshell potential
 *
 *****************************************************************/

void softshell_value(real r, real *p, real *f)
{
  static real x, y;

  x = p[0] / r;
  y = p[1];
#ifndef ACML
  vdPow(1, &x, &y, f);
#else
  *f = fastpow(x, y);
#endif
}

/****************************************************************
 *
 * eopp_exp potential
 *
 ****************************************************************/

void eopp_exp_value(real r, real *p, real *f)
{
  static real power;

#ifndef ACML
  vdPow(1, &r, &p[3], &power);
#else
  power = fastpow(r, p[3]);
#endif

  *f = p[0] * exp(-p[1] * r) + (p[2] / power) * cos(p[4] * r + p[5]);
}

/****************************************************************
 *
 * meopp potential
 *
 ****************************************************************/

void meopp_value(real r, real *p, real *f)
{
  static real x[2], y[2], power[2];

  x[0] = r - p[6];
  x[1] = r;
  y[0] = p[1];
  y[1] = p[3];
#ifndef ACML
  vdPow(2, x, y, power);
#else
  power[0] = fastpow(x[0], y[0]);
  power[1] = fastpow(x[1], y[1]);
#endif

  *f = p[0] / power[0] + (p[2] / power[1]) * cos(p[4] * r + p[5]);
}

/****************************************************************
 *
 * power_decay potential
 *
 ****************************************************************/

void power_decay_value(real r, real *p, real *f)
{
  static real x, y, power;

  x = 1. / r;
  y = p[1];
#ifndef ACML
  vdPow(1, &x, &y, &power);
#else
  power = fastpow(x, y);
#endif

  *f = p[0] * power;
}

/****************************************************************
 *
e pot * exp_decay potential
 *
 ****************************************************************/

void exp_decay_value(real r, real *p, real *f)
{
  *f = p[0] * exp(-p[1] * r);
}

/****************************************************************
 *
 * pohlong potential
 *
 ****************************************************************/

void pohlong_value(real r, real *p, real *f)
{
  real  power;

  if (r == 0)
    *f = 0;
  else {
#ifndef ACML
    vdPow(1, &r, &p[1], &power);
#else
    power = fastpow(r, p[1]);
#endif

    *f = p[0] * (1 - p[1] * log(r)) * power + p[2] * r;
  }
}

/****************************************************************
 *
 * parabola potential
 *
 ****************************************************************/

void parabola_value(real r, real *p, real *f)
{
  *f = SQR(r) * p[0] + r * p[1] + p[2];
}

/****************************************************************
 *
 * chantasiriwan (csw) and milstein potential
 *
 ****************************************************************/

void csw_value(real r, real *p, real *f)
{
  real  power;

#ifndef ACML
  vdPow(1, &r, &p[3], &power);
#else
  power = fastpow(r, p[3]);
#endif

  *f = (1 + p[0] * cos(p[2] * r) + p[1] * sin(p[2] * r)) / power;
}

/****************************************************************
 *
 * universal embedding function
 *
 ****************************************************************/

void universal_value(real r, real *p, real *f)
{
  static real x[2], y[2], power[2];

  x[0] = r;
  x[1] = r;
  y[0] = p[1];
  y[1] = p[2];
#ifndef ACML
  vdPow(2, x, y, power);
#else
  power[0] = fastpow(x[0], y[0]);
  power[1] = fastpow(x[1], y[1]);
#endif

  *f =
    p[0] * (p[2] / (p[2] - p[1]) * power[0] - p[1] / (p[2] - p[1]) * power[1]) +
    p[3] * r;
}

/****************************************************************
 *
 * constant function
 *
 ****************************************************************/

void const_value(real r, real *p, real *f)
{
  *f = *p;
}

/****************************************************************
 *
 * square root function
 *
 ****************************************************************/

void sqrt_value(real r, real *p, real *f)
{
  *f = p[0] * sqrt(r / p[1]);
}

/****************************************************************
 *
 * mexp_decay potential
 *
 ****************************************************************/

void mexp_decay_value(real r, real *p, real *f)
{
  *f = p[0] * exp(-p[1] * (r - p[2]));
}

/****************************************************************
 *
 * streitz-mintmire (strmm) potential
 *
 ****************************************************************/

void strmm_value(real r, real *p, real *f)
{
  real  r_0 = r - p[4];

  *f =
    2 * p[0] * exp(-p[1] / 2 * r_0) - p[2] * (1 +
    p[3] * r_0) * exp(-p[3] * r_0);
}

/****************************************************************
 *
 * double morse potential
 *
 ****************************************************************/

void double_morse_value(real r, real *p, real *f)
{
  *f =
    (p[0] * (exp(-2 * p[1] * (r - p[2])) - 2 * exp(-p[1] * (r - p[2]))) +
    p[3] * (exp(-2 * p[4] * (r - p[5])) - 2 * exp(-p[4] * (r - p[5])))) + p[6];
}

/****************************************************************
 *
 * double exp potential
 *
 ****************************************************************/

void double_exp_value(real r, real *p, real *f)
{
  *f = (p[0] * exp(-p[1] * SQR(r - p[2])) + exp(-p[3] * (r - p[4])));
}

/****************************************************************
 *
 * poly 5 potential
 *
 ****************************************************************/

void poly_5_value(real r, real *p, real *f)
{
  real  dr = (r - 1) * (r - 1);
  *f =
    p[0] + 0.5 * p[1] * dr + p[2] * (r - 1) * dr + p[3] * SQR(dr) +
    p[4] * SQR(dr) * (r - 1);
}

/****************************************************************
 *
 * cbb potential, from C. B. Basak
 *
 * see http://dx.doi.org/10.1016/S0925-8388(03)00350-5
 *
 ****************************************************************/

void cbb_value(real r, real *p, real *f)
{
  real  r6 = r * r * r;
  r6 *= r6;

  *f =
    p[0] * p[1] / r + p[2] * (p[5] + p[6]) * exp((p[3] + p[4] - r) / (p[5] +
      p[6])) - p[7] * p[8] / r6 + p[2] * p[9] * (exp(-2 * p[10] * (r - p[11])) -
    2 * exp(-p[10] * (r - p[11])));
}

/****************************************************************
 *
 * exp_plus potential
 *
 ****************************************************************/

void exp_plus_value(real r, real *p, real *f)
{
  *f = p[0] * exp(-p[1] * r) + p[2];
}

/****************************************************************
 *
 * mishin potential
 *
 ****************************************************************/

void mishin_value(real r, real *p, real *f)
{
  real  z = r - p[3];
  real  temp = exp(-p[5] * r);
  real  power;

#ifndef ACML
  vdPow(1, &z, &p[4], &power);
#else
  power = fastpow(z, p[4]);
#endif /* ACML */

  *f = p[0] * power * temp * (1 + p[1] * temp) + p[2];
}

/****************************************************************
 *
 * gen_lj potential, generalized lennard-jones
 *
 ****************************************************************/

void gen_lj_value(real r, real *p, real *f)
{
  static real x[2], y[2], power[2];

  x[0] = r / p[3];
  x[1] = x[0];
  y[0] = p[1];
  y[1] = p[2];
#ifndef ACML
  vdPow(2, x, y, power);
#else
  power[0] = fastpow(x[0], y[0]);
  power[1] = fastpow(x[1], y[1]);
#endif

  *f = p[0] / (p[2] - p[1]) * (p[2] / power[0] - p[1] / power[1]) + p[4];
}

/****************************************************************
 *
 * gljm potential, generalized lennard-jones + mishin potential
 *
 ****************************************************************/

void gljm_value(real r, real *p, real *f)
{
  real  x[3], y[3], power[3];

  x[0] = r / p[3];
  x[1] = x[0];
  x[2] = r - p[9];
  y[0] = p[1];
  y[1] = p[2];
  y[2] = p[10];
#ifndef ACML
  vdPow(3, x, y, power);
#else
  power[0] = fastpow(x[0], y[0]);
  power[1] = fastpow(x[1], y[1]);
  power[2] = fastpow(x[2], y[2]);
#endif /* ACML */

  real  temp = exp(-p[11] * power[2]);

  *f =
    p[0] / (p[2] - p[1]) * (p[2] / power[0] - p[1] / power[1]) + p[4] +
    p[5] * (p[6] * power[2] * temp * (1 + p[7] * temp) + p[8]);
}

/****************************************************************
 *
 * bond-stretching function of vashishta potential (f_c)
 *
 ****************************************************************/

void vas_value(real r, real *p, real *f)
{
  *f = exp(p[0] / (r - p[1]));
}

/****************************************************************
 *
 * original pair contributions of vashishta potential (V_2)
 *
 ****************************************************************/

void vpair_value(real r, real *p, real *f)
{
  real  x[7], y, z;

  y = r;
  z = p[1];
  vdPow(1, &y, &z, &x[0]);
  x[1] = r * r;
  x[2] = x[1] * x[1];
  x[3] = p[2] * p[2];
  x[4] = p[3] * p[3];
  x[5] = p[4] * x[4] + p[5] * x[3];
  x[6] = exp(-r / p[6]);

  *f = 14.4 * (p[0] / x[0] + p[2] * p[3] / r - 0.5 * x[5] / x[2] * x[6]);
}

/****************************************************************
 *
 * template for new potential function called mypotential
 * for further information plase have a look at the online documentation
 *
 * http://www.itap.physik.uni-stuttgart.de/~imd/potfit/potfit.html
 *
 ****************************************************************/

/* template for new function */

/****************************************************************
 *
 * newpot potential
 *
 ****************************************************************/

void newpot_value(real r, real *p, real *f)
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

real cutoff(real r, real r0, real h)
{
  if ((r - r0) > 0)
    return 0;

  real  val = 0;

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

int apot_check_params(real *params)
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

    /* jump to next potential */
    j += 1 + apot_table.n_par[i] + smooth_pot[i];
  }
  return 0;
}

/****************************************************************
 *
 * punish analytic potential for bad habits
 *
 ****************************************************************/

real apot_punish(real *params, real *forces)
{
  int   i, j;
  real  x, tmpsum = 0, min, max;

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
#endif

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

real apot_grad(real r, real *p, void (*function) (real, real *, real *))
{
  real  a, b, h = 0.0001;

  function(r + h, p, &a);
  function(r - h, p, &b);

  return (a - b) / (2 * h);
}

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
    fprintf(stderr, "\npotential #%d (type=%s, smooth=%d)\n", i + 1,
      apot_table.names[i], smooth_pot[i]);
    fprintf(stderr, "begin=%f end=%f\n", apot_table.begin[i],
      apot_table.end[i]);
    for (j = 0; j < apot_table.n_par[i]; j++) {
      fprintf(stderr, "parameter %d: name=%s value=%f min=%f max=%f\n", j + 1,
	apot_table.param_name[i][j], apot_table.values[i][j],
	apot_table.pmin[i][j], apot_table.pmax[i][j]);
    }
  }
#ifdef PAIR
  if (!enable_cp) {
    fprintf(stderr, "\nchemical potentials are DISABLED!\n");
  } else {
    fprintf(stderr, "\nchemical potentials:\n");
    for (i = 0; i < ntypes; i++)
      fprintf(stderr, "cp_%d=%f min=%f max=%f\n", i, apot_table.chempot[i],
	apot_table.pmin[apot_table.number][i],
	apot_table.pmax[apot_table.number][i]);
    if (compnodes > 0) {
      if (ntypes == 2) {
	fprintf(stderr, "composition nodes:\n");
	for (j = 0; j < compnodes; j++)
	  fprintf(stderr, "composition=%f value=%f min=%f max=%f\n",
	    compnodelist[j], apot_table.chempot[ntypes + j],
	    apot_table.pmin[apot_table.number][ntypes + j],
	    apot_table.pmax[apot_table.number][ntypes + j]);
      }
    }
  }
#endif /* PAIR */
  exit(EXIT_FAILURE);
}

#endif /* DEBUG */

#ifdef COULOMB

/******************************************************************************
*
* ms potential + first derivative
*
******************************************************************************/

void ms_init(real r, real *pot, real *grad, real *p)
{
  static real x[4];

  x[0] = 1 - r / p[2];
  x[1] = exp(p[1] * x[0]);
  x[2] = exp(p[1] * x[0] / 2);
  x[3] = p[0] * p[1] / (r * p[2]);

  *pot = p[0] * (x[1] - 2 * x[2]);
  *grad = x[3] * (-x[1] + x[2]);
}

/******************************************************************************
*
* buckingham potential + first derivative
*
******************************************************************************/

void buck_init(real r, real *pot, real *grad, real *p)
{
  static real x[3];

  x[0] = SQR(p[1]) / SQR(r);
  x[1] = p[2] * x[0] * x[0] * x[0];
  x[2] = p[0] * exp(-r / p[1]);

  *pot = x[2] - x[1];
  *grad = - x[2] / p[1] + 6 * p[1] * x[1] / r;
}

/******************************************************************************
*
* shifted ms potential
*
******************************************************************************/

void ms_shift(real r, real *p, real *f)
{
  static real pot, grad, pot_cut, grad_cut;

  ms_init(r, &pot, &grad, p);
  ms_init(dp_cut, &pot_cut, &grad_cut, p);

  *f = pot - pot_cut - r * (r - dp_cut) * grad_cut;
}

/******************************************************************************
*
* shifted buckingham potential
*
******************************************************************************/

void buck_shift(real r, real *p, real *f)
{
  static real pot, grad, pot_cut, grad_cut;

  buck_init(r, &pot, &grad, p);
  buck_init(dp_cut, &pot_cut, &grad_cut, p);

  *f = pot - pot_cut - r * (r - dp_cut) * grad_cut;
}

/******************************************************************************
*
* tail of electrostatic potential and first two derivatives
*
******************************************************************************/

void elstat_value(real r, real *ftail, real *gtail, real *ggtail)
{
  static real x[4];

  x[0] = r * r;
  x[1] = dp_kappa * dp_kappa;
  x[2] = 2 * dp_eps * dp_kappa / sqrt(M_PI);
  x[3] = exp(-x[0] * x[1]);

  *ftail = dp_eps * erfc(dp_kappa * r) / r;
  *gtail = -(*ftail + x[2] * x[3]) / x[0];
  *ggtail = (2 * x[1] * x[2] * x[3] - *gtail * 3) / x[0];
}

/******************************************************************************
*
* shifted tail of coloumb potential
*
******************************************************************************/

void elstat_shift(real r, real *fnval_tail, real *grad_tail, real *ggrad_tail)
{
  static real ftail, gtail, ggtail, ftail_cut, gtail_cut, ggtail_cut;
  static real x[3];

  x[0] = r * r;
  x[1] = dp_cut * dp_cut;
  x[2] = x[0] - x[1];

  elstat_value(r, &ftail, &gtail, &ggtail);
  elstat_value(dp_cut, &ftail_cut, &gtail_cut, &ggtail_cut);

  *fnval_tail = ftail - ftail_cut - x[2] * gtail_cut / 2;
  *grad_tail = gtail - gtail_cut;
  *ggrad_tail = 0.;
#ifdef DIPOLE
  *fnval_tail -= x[2] * x[2] * ggtail_cut / 8;
  *grad_tail -= x[2] * ggtail_cut / 2;
  *ggrad_tail = ggtail - ggtail_cut;	// ? richtig so? hier alles checken im Falle Dipole!
#endif
}

#endif /* COULOMB */
#ifdef DIPOLE

/******************************************************************************
*
* short-range part of dipole moments 
*
******************************************************************************/

real shortrange_value(real r, real a, real b, real c)
{
  static real x[5];

  x[0] = b * r;
  x[1] = x[0] * x[0];
  x[2] = x[1] * x[0];
  x[3] = x[1] * x[1];
  x[4] = 1 + x[0] + x[1] / 2 + x[2] / 6 + x[3] / 24;

  return a * c * x[4] * exp(-x[0]) / dp_eps;
}

/******************************************************************************
*
* tail of additional short-range contribution to energy and forces  
*
******************************************************************************/

void shortrange_term(real r, real b, real c, real *srval_tail,
  real *srgrad_tail)
{
  static real x[6];

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
#endif /* APOT */
