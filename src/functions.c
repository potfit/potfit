/****************************************************************
 *
 * functions.c: Routines and function calls used for analytic potentials
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

#include "potfit.h"

#include "functions.h"
#include "memory.h"
#include "utils.h"

// eopp = 0
// universal = 1
#define NUM_PUNISH_FUNCTIONS 2

/****************************************************************
  function table
    stores all available analytic potentials
****************************************************************/

struct {
  char** name;             // identifier of the potential
  int* num_params;         // number of parameters
  fvalue_pointer* fvalue;  // function pointer
  int num_functions;       // number of analytic function prototypes
  int** punish_index;      // array to index which functions may be punished
} function_table;

/****************************************************************
  initialize_analytic_potentials
    initialize the function_table for analytic potentials
****************************************************************/

void initialize_analytic_potentials(void)
{
#define FUNCTION(name, npar) add_potential(#name, npar, &name##_value)

#include "functions.itm"

#undef FUNCTION
  function_table.punish_index =
      (int**)Malloc(NUM_PUNISH_FUNCTIONS * sizeof(int*));
  for (int i = 0; i < NUM_PUNISH_FUNCTIONS; ++i)
    function_table.punish_index[i] = (int*)Malloc(sizeof(int));
}

/****************************************************************
  add_potential
    add analytic function to function_table
****************************************************************/

void add_potential(const char* name, int npar, fvalue_pointer function)
{
  const int k = function_table.num_functions;

  // only add potentials with unused names
  for (int i = 0; i < k; i++) {
    if (strcmp(function_table.name[i], name) == 0)
      error(1, "There already is a potential with the name \"%s\".\n", name);
  }

  // allocate memory
  function_table.name =
      (char**)Realloc(function_table.name, (k + 1) * sizeof(char*));
  function_table.name[k] = (char*)Malloc((strlen(name) + 1) * sizeof(char));
  function_table.num_params =
      (int*)Realloc(function_table.num_params, (k + 1) * sizeof(int));
  function_table.fvalue = (fvalue_pointer*)Realloc(
      function_table.fvalue, (k + 1) * sizeof(fvalue_pointer));

  // assign values
  sprintf(function_table.name[k], "%s", name);
  function_table.num_params[k] = npar;
  function_table.fvalue[k] = function;

  function_table.num_functions++;
}

/****************************************************************
  apot_get_num_parameters
    return the number of parameters for a specific analytic potential
****************************************************************/

int apot_get_num_parameters(const char* name)
{
  for (int i = 0; i < function_table.num_functions; i++) {
    if (strcmp(function_table.name[i], name) == 0)
      return function_table.num_params[i];
  }

  return -1;
}

/****************************************************************
  apot_assign_function_pointers
    assign function pointers to corresponding functions
****************************************************************/

int apot_assign_function_pointers(apot_table_t* apt)
{
  for (int i = 0; i < apt->number; i++) {
    for (int j = 0; j < function_table.num_functions; j++) {
      if (strcmp(apt->names[i], function_table.name[j]) == 0) {
        apt->fvalue[i] = function_table.fvalue[j];
        apot_assign_punish_functions(apt->names[i], i);
        break;
      }
      if (j == function_table.num_functions - 1)
        return -1;
    }
  }

  return 0;
}

/****************************************************************
  apot_assign_punish_functions
    assign function punishment indices
****************************************************************/

void apot_assign_punish_functions(char const* name, int index)
{
  // eopp is index 0
  if (strncmp(name, "eopp", 4) == 0) {
    int num_pot = ++function_table.punish_index[0][0];
    function_table.punish_index[0] =
        (int*)Realloc(function_table.punish_index[0], (num_pot + 1) * sizeof(int));
    function_table.punish_index[0][num_pot] = index;
  } else
      // universal is index 1
      if (strncmp(name, "universal", 9) == 0) {
    int num_pot = ++function_table.punish_index[1][0];
    function_table.punish_index[1] =
        (int*)Realloc(function_table.punish_index[1], (num_pot + 1) * sizeof(int));
    function_table.punish_index[1][num_pot] = index;
  }
}

/****************************************************************
  check_correct_apot_functions
    check for special functions needed by certain potential models
****************************************************************/

void check_correct_apot_functions(void)
{
#if defined(STIWEB)
  // paircol is not yet defined at this point
  int pcol = (g_param.ntypes * (g_param.ntypes + 1)) / 2;

  // check for the correct function types for SW potential
  for (int i = 0; i < pcol; i++) {
    if (strcmp(g_pot.apot_table.names[i], "stiweb_2") != 0)
      error(1, "Only stiweb_2 potential is allowed for the %d. potential!\n",
            i + 1);
    if (strcmp(g_pot.apot_table.names[pcol + i], "stiweb_3") != 0)
      error(1, "Only stiweb_3 potential is allowed for the %d. potential!\n",
            pcol + i + 1);
  }

  if (strcmp(g_pot.apot_table.names[2 * pcol], "lambda") != 0)
    error(1,
          "The last potential for Stillinger-Weber has to be of the \"lambda\" "
          "type!\n");

  // make sure the cutoff parameters (a1,a2) can be optimized correctly
  for (int i = 0; i < pcol; i++) {
    if (g_pot.apot_table.pmax[i][5] > g_pot.apot_table.end[i]) {
      error(0,
            "The upper bound for the parameter a1 exceeds the cutoff radius in "
            "potential %d.\n",
            i + 1);
      error(1, "a1 needs to be less or equal to the potential cutoff.\n");
    }
    if (g_pot.apot_table.pmax[pcol + i][1] > g_pot.apot_table.end[pcol + i]) {
      error(0,
            "The upper bound for the parameter a2 exceeds the cutoff radius in "
            "potential %d.\n",
            i + 1);
      error(1, "a1 needs to be less or equal to the potential cutoff.\n");
    }
  }
#endif  // STIWEB

#if defined(TERSOFF)
  // paircol is not yet defined at this point
  int pcol = (g_param.ntypes * (g_param.ntypes + 1)) / 2;

#if !defined(TERSOFFMOD)
  // check for the correct function types for TERSOFF potential
  for (int i = 0; i < pcol; i++) {
    if (strcmp(g_pot.apot_table.names[i], "tersoff_pot") != 0)
      error(1, "Only tersoff_pot potential is allowed for the %d. potential!\n",
            i + 1);
  }
  for (int i = 0; i < g_param.ntypes * (g_param.ntypes - 1) / 2.0; i++) {
    if (strcmp(g_pot.apot_table.names[pcol + i], "tersoff_mix") != 0)
      error(1, "Only tersoff_mix potential is allowed for the %d. potential!\n",
            i + 1);
  }

  // make sure the cutoff parameters (R, S) can be optimized correctly
  for (int i = 0; i < pcol; i++) {
    if (g_pot.apot_table.pmax[i][9] < g_pot.apot_table.pmin[i][10]) {
      error(0,
            "The upper bound for the parameter S is smaller than the lower "
            "bound\n");
      error(0, "for the parameter R in potential %d.\n", i + 1);
      error(1,
            "Please change it, that the condition R < S can be fulfilled.\n");
    }
  }
#else
  // check for the correct function types for TERSOFFMOD potential
  for (int i = 0; i < pcol; i++) {
    if (strcmp(g_pot.apot_table.names[i], "tersoff_mod_pot") != 0)
      error(1, "Only tersoff_pot potential is allowed for the %d. potential!\n",
            i + 1);
  }

  // make sure the cutoff parameters (R, S) can be optimized correctly
  for (int i = 0; i < pcol; i++) {
    if (g_pot.apot_table.pmax[i][15] < g_pot.apot_table.pmin[i][14]) {
      error(0,
            "The upper bound for the parameter R2 is smaller than the lower "
            "bound\n");
      error(0, "for the parameter R1 in potential %d.\n", i + 1);
      error(1,
            "Please change it, that the condition R1 < R2 can be fulfilled.\n");
    }
  }

#endif  // !TERSOFFMOD
#endif  // TERSOFF
}

/****************************************************************
  check analytic parameters for special conditions
****************************************************************/

int apot_check_params(double* params)
{
  int j = 2;
  int k = 0;

  for (int i = 0; i < g_pot.apot_table.number; i++) {
    // last parameter of eopp potential is 2 pi periodic
    if (strcmp(g_pot.apot_table.names[i], "eopp") == 0) {
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

    // the third parameter of csw2 potential is 2 pi periodic
    if (strcmp(g_pot.apot_table.names[i], "csw2") == 0) {
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

#if defined(TERSOFF)
    // the parameter S has to be greater than the parameter R
    // switch them if this is not the case
    if (strcmp(g_pot.apot_table.names[i], "tersoff_pot") == 0) {
      k = j + 9;
      if (params[k] < params[k + 1]) {
        double temp = params[k];
        params[k] = params[k + 1];
        params[k + 1] = temp;
      }
    }
#endif  // TERSOFF

    // jump to next potential
    j += 2 + apot_get_num_parameters(g_pot.apot_table.names[i]) +
         g_pot.smooth_pot[i];
  }

  return 0;
}

/****************************************************************
  apot_cutoff
    function for smooth cutoff radius
****************************************************************/

double apot_cutoff(const double r, const double r0, const double h)
{
  if ((r - r0) >= 0)
    return 0.0;

  double val = (r - r0) / h;
  val *= val;
  val *= val;

  return val / (1.0 + val);
}

/****************************************************************
  apot_gradient
    calculate gradient for analytic potential
****************************************************************/

double apot_gradient(const double r, const double* p, fvalue_pointer func)
{
  double a = 0.0;
  double b = 0.0;
  double h = 0.0001;

  func(r + h, p, &a);
  func(r - h, p, &b);

  return (a - b) / (2.0 * h);
}

/****************************************************************
  apot_punish
    punish analytic potential for bad habits
****************************************************************/

double apot_punish(double* params, double* forces)
{
  double tmpsum = 0.0;

  // loop over individual parameters
  for (int i = 0; i < g_calc.ndim; i++) {
    double min =
        g_pot.apot_table
            .pmin[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]];
    double max =
        g_pot.apot_table
            .pmax[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]];
    // punishment for out of bounds
    if (params[g_pot.opt_pot.idx[i]] < min) {
      double x = params[g_pot.opt_pot.idx[i]] - min;
      tmpsum += APOT_PUNISH * x * x;
      forces[g_calc.punish_par_p + i] = APOT_PUNISH * x * x;
    } else if (params[g_pot.opt_pot.idx[i]] > max) {
      double x = params[g_pot.opt_pot.idx[i]] - max;
      tmpsum += APOT_PUNISH * x * x;
      forces[g_calc.punish_par_p + i] = APOT_PUNISH * x * x;
    }
  }

  // eopp (index 0)
  // punish eta_1 < eta_2 for eopp function
  for (int i = 1; i <= function_table.punish_index[0][0]; ++i) {
    int idx = g_pot.opt_pot.first[function_table.punish_index[0][i]];
    double x = params[idx + 1] - params[idx + 3];
    if (x < 0) {
      forces[g_calc.punish_pot_p + function_table.punish_index[0][i]] =
          g_param.apot_punish_value * (1 + x) * (1 + x);
      tmpsum += g_param.apot_punish_value * (1 + x) * (1 + x);
    }
  }

  // universal (index 1)
  // punish m=n for universal embedding function
  for (int i = 1; i <= function_table.punish_index[1][0]; ++i) {
    int idx = g_pot.opt_pot.first[function_table.punish_index[1][i]];
    double x = fabs(params[idx + 2] - params[idx + 1]);
    if (x < 1e-6) {
      forces[g_calc.punish_pot_p + function_table.punish_index[1][i]] =
          g_param.apot_punish_value / (x * x);
      tmpsum += g_param.apot_punish_value / (x * x);
    }
  }

  return tmpsum;
}

#if defined(COULOMB)

/****************************************************************
  elstat_value
    tail of electrostatic potential and first two derivatives
****************************************************************/

void elstat_value(double r, double dp_kappa, double* ftail, double* gtail,
                  double* ggtail)
{
  double x[4];

  x[0] = r * r;
  x[1] = dp_kappa * dp_kappa;
  x[2] = 2 * DP_EPS * dp_kappa / sqrt(M_PI);
  x[3] = exp(-x[0] * x[1]);

  *ftail = DP_EPS * erfc(dp_kappa * r) / r;
  *gtail = -(*ftail + x[2] * x[3]) / x[0];                /* 1/r df/dr */
  *ggtail = (2 * x[1] * x[2] * x[3] - *gtail * 3) / x[0]; /* 1/r dg/dr */
}

/****************************************************************
  elstat_shift
    shifted tail of coulomb potential
****************************************************************/

void elstat_shift(double r, double dp_kappa, double* fnval_tail,
                  double* grad_tail, double* ggrad_tail)
{
  double ftail, gtail, ggtail, ftail_cut, gtail_cut, ggtail_cut;
  double x[3];

  x[0] = r * r;
  x[1] = g_config.dp_cut * g_config.dp_cut;
  x[2] = x[0] - x[1];

  elstat_value(r, dp_kappa, &ftail, &gtail, &ggtail);
  elstat_value(g_config.dp_cut, dp_kappa, &ftail_cut, &gtail_cut, &ggtail_cut);

  *fnval_tail = ftail - ftail_cut - x[2] * gtail_cut / 2;
  *grad_tail = gtail - gtail_cut;
  *ggrad_tail = 0.0;
#if defined(DIPOLE)
  *fnval_tail -= x[2] * x[2] * ggtail_cut / 8;
  *grad_tail -= x[2] * ggtail_cut / 2;
  *ggrad_tail = ggtail - ggtail_cut;
#endif  // DIPOLE
}

/***************************************************************
 *
 * damped shifted force Coulomb potential
 * http://dx.doi.org/10.1063/1.2206581
 *
 ****************************************************************/

#if defined(DSF)
void elstat_dsf(double r, double dp_kappa, double* fnval_tail,
                double* grad_tail, double* ggrad_tail)
{
  static double ftail, gtail, ggtail, ftail_cut, gtail_cut, ggtail_cut;
  static double x[3];

  x[0] = r * r;
  x[1] = g_config.dp_cut * g_config.dp_cut;
  x[2] = x[0] - x[1];

  elstat_value(r, dp_kappa, &ftail, &gtail, &ggtail);
  elstat_value(g_config.dp_cut, dp_kappa, &ftail_cut, &gtail_cut, &ggtail_cut);

  *fnval_tail =
      ftail - ftail_cut - (r - g_config.dp_cut) * gtail_cut * g_config.dp_cut;
  *grad_tail = gtail - gtail_cut * g_config.dp_cut / r; /*  1/r dV/r */
  *ggrad_tail = 0.0;
#ifdef DIPOLE
  *fnval_tail -= x[2] * x[2] * ggtail_cut / 8;
  *grad_tail -= x[2] * ggtail_cut / 2;
  *ggrad_tail = ggtail - ggtail_cut;
#endif /* DIPOLE */
}
#endif  // DSF

/****************************************************************
  init_tails
    calculate tail of coulomb-potential and its first derivative
****************************************************************/

void init_tails(double dp_kappa)
{
  for (int i = 0; i < g_config.natoms; i++) {
    for (int j = 0; j < g_config.atoms[i].num_neigh; j++) {
#if defined(DSF)
      elstat_dsf(g_config.atoms[i].neigh[j].r, dp_kappa,
                 &g_config.atoms[i].neigh[j].fnval_el,
                 &g_config.atoms[i].neigh[j].grad_el,
                 &g_config.atoms[i].neigh[j].ggrad_el);
#else
      elstat_shift(g_config.atoms[i].neigh[j].r, dp_kappa,
                   &g_config.atoms[i].neigh[j].fnval_el,
                   &g_config.atoms[i].neigh[j].grad_el,
                   &g_config.atoms[i].neigh[j].ggrad_el);
#endif  // DSF
    }
  }
}

#endif  // COULOMB

#if defined(DIPOLE)

/****************************************************************
  shortrange_value
    short-range part of dipole moments
****************************************************************/

double shortrange_value(double r, double a, double b, double c)
{
  double x[5];

  x[0] = b * r;
  x[1] = x[0] * x[0];
  x[2] = x[1] * x[0];
  x[3] = x[1] * x[1];
  x[4] = 1 + x[0] + x[1] / 2 + x[2] / 6 + x[3] / 24;

  return a * c * x[4] * exp(-x[0]) / DP_EPS;
}

/****************************************************************
  shortrange_term
    tail of additional short-range contribution to energy and forces
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

  *srval_tail = c * x[4] * x[5] / DP_EPS;
  *srgrad_tail = -c * b * x[3] * x[5] / (24 * DP_EPS * r);
}

#endif  // DIPOLE

#if defined(DEBUG)

/****************************************************************
  debug function to print the potentials
****************************************************************/

void debug_apot()
{
  fflush(stdout);
  fprintf(stderr, "\n##############################################\n");
  fprintf(stderr, "###########      DEBUG OUTPUT      ###########\n");
  fprintf(stderr, "##############################################\n");
  fprintf(stderr, "\nThere are %d potentials with a total of %d parameters.\n",
          g_pot.apot_table.number, g_pot.apot_table.total_par);
  for (int i = 0; i < g_pot.apot_table.number; i++) {
    fprintf(stderr, "\npotential #%d (type=%s, smooth=%d)\n", i + 1,
            g_pot.apot_table.names[i], g_pot.smooth_pot[i]);
    fprintf(stderr, "begin=%f end=%f\n", g_pot.apot_table.begin[i],
            g_pot.apot_table.end[i]);
    for (int j = 0; j < g_pot.apot_table.n_par[i]; j++) {
      fprintf(stderr, "parameter %d: name=%s value=%f min=%f max=%f\n", j + 1,
              g_pot.apot_table.param_name[i][j], g_pot.apot_table.values[i][j],
              g_pot.apot_table.pmin[i][j], g_pot.apot_table.pmax[i][j]);
    }
  }

#if defined(PAIR)
  if (!g_param.enable_cp) {
    fprintf(stderr, "\nchemical potentials are DISABLED!\n");
  } else {
    fprintf(stderr, "\nchemical potentials:\n");
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(stderr, "cp_%d=%f min=%f max=%f\n", i,
              g_pot.apot_table.chempot[i],
              g_pot.apot_table.pmin[g_pot.apot_table.number][i],
              g_pot.apot_table.pmax[g_pot.apot_table.number][i]);
    if (g_param.compnodes > 0 && g_param.ntypes == 2) {
      fprintf(stderr, "composition nodes:\n");
      for (int j = 0; j < g_param.compnodes; j++)
        fprintf(
            stderr, "composition=%f value=%f min=%f max=%f\n",
            g_pot.compnodelist[j], g_pot.apot_table.chempot[g_param.ntypes + j],
            g_pot.apot_table.pmin[g_pot.apot_table.number][g_param.ntypes + j],
            g_pot.apot_table.pmax[g_pot.apot_table.number][g_param.ntypes + j]);
    }
  }
#endif  // PAIR

  exit(EXIT_FAILURE);
}

#endif  // DEBUG
