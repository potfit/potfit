/****************************************************************
 *
 * potential_input.c: Routines for reading a potential table
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

#include "chempot.h"
#include "functions.h"
#include "kim.h"
#include "memory.h"
#include "potential_input.h"
#include "utils.h"

void read_pot_line_F(char const* pbuf, potential_state* pstate);
void read_pot_line_T(char const* pbuf, potential_state* pstate);
void read_pot_line_I(char* pbuf, potential_state* pstate);
void read_pot_line_G(char* pbuf, potential_state* pstate);
void read_pot_line_C(char* pbuf, potential_state* pstate);

void allocate_memory_for_potentials(potential_state* pstate);

void calculate_cutoffs();
void read_maxch_file();

#if defined(APOT)
#if !defined(KIM)
void read_pot_table0(char const* potential_filename, FILE* pfile);
#else
void read_pot_table5(char const* potential_filename, FILE* pfile);
#endif // KIM
#else
void read_pot_table3(char const* potential_filename, FILE* pfile,
                     potential_state* pstate);
void read_pot_table4(char const* potential_filename, FILE* pfile,
                     potential_state* pstate);
#endif // APOT

/****************************************************************
 *
 * read potential tables
 *
 ****************************************************************/

void read_pot_table(char const* potential_filename)
{
  char buffer[1024];
  char* res = NULL;

  potential_state pstate;

  memset(&pstate, 0, sizeof(pstate));

  pstate.filename = potential_filename;

  // set paircol to the number of pair potentials
  g_calc.paircol = (g_param.ntypes * (g_param.ntypes + 1)) / 2;

  // open file
  FILE* pfile = fopen(potential_filename, "r");

  if (NULL == pfile)
    error(1, "Could not open file %s\n", potential_filename);

  printf("Starting to read the potential file:\n");

  // read the header
  do {
    // read one line
    res = fgets(buffer, 1024, pfile);

    if (NULL == res)
      error(1, "Unexpected end of file in %s\n", potential_filename);

    // check if it is a header line
    if (buffer[0] != '#')
      error(1, "Header corrupt in file %s\n", potential_filename);

    switch (buffer[1]) {
      case 'C':
        read_pot_line_C(buffer, &pstate);
        break;
      case 'E':
        pstate.end_header = 1;
        break;
      case 'G':
        read_pot_line_G(buffer, &pstate);
        break;
      case 'I':
        read_pot_line_I(buffer, &pstate);
        break;
      case 'F':
        read_pot_line_F(buffer, &pstate);
        break;
      case 'T':
        read_pot_line_T(buffer, &pstate);
      case '#':
        continue;
      default:
        warning("Ignoring potential header line starting with #%c.\n", buffer[1]);
    }

  } while (!pstate.end_header);

  // do we have a format in the header?
  if (!pstate.have_format)
    error(1, "Format not specified in header of potential file %s\n",
          potential_filename);

  allocate_memory_for_potentials(&pstate);

  switch (g_pot.format_type) {
    case POTENTIAL_FORMAT_UNKNOWN:
      error(1, "Unknown potential format detected! (%s:%d)\n", __FILE__,
            __LINE__);
    case POTENTIAL_FORMAT_ANALYTIC:
#if defined(APOT) && !defined(KIM)
      read_pot_table0(potential_filename, pfile);
#endif  // APOT
      break;
    case POTENTIAL_FORMAT_TABULATED_EQ_DIST:
#if !defined(APOT)
      read_pot_table3(potential_filename, pfile, &pstate);
#endif  // !APOT
      break;
    case POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST:
#if !defined(APOT)
      read_pot_table4(potential_filename, pfile, &pstate);
#endif  // !APOT
      break;
    case POTENTIAL_FORMAT_KIM:
#if defined(APOT) && defined(KIM)
      read_pot_table5(potential_filename, pfile);
#endif  // KIM
      break;
  }

  fclose(pfile);

  printf("Reading potential file >> %s << ... done\n", potential_filename);

  calculate_cutoffs();

  read_maxch_file();
}

/****************************************************************
  read_pot_line_F
****************************************************************/

void read_pot_line_F(char const* pbuf, potential_state* pstate)
{
  int format = POTENTIAL_FORMAT_UNKNOWN;

  int rval = sscanf(pbuf + 2, "%d %d", &format, &pstate->num_pots);
  if (rval != 2)
    error(1, "Corrupt format header line in file %s\n", pstate->filename);

  switch (format) {
    case POTENTIAL_FORMAT_ANALYTIC:
      printf(" - Potential file format %d detected: analytic potentials\n", format);
      break;
    case POTENTIAL_FORMAT_TABULATED_EQ_DIST:
      printf(" - Potential file format %d detected: tabulated eqdist\n", format);
      break;
    case POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST:
      printf(" - Potential file format %d detected: tabulated non-eqdist\n", format);
      break;
    case POTENTIAL_FORMAT_KIM:
      printf(" - Potential file format %d detected: KIM\n", format);
      break;
    case POTENTIAL_FORMAT_UNKNOWN:
    default:
      break;
  }

#if defined(APOT)
  if (format != POTENTIAL_FORMAT_ANALYTIC)
    error(1, "This potfit binary only supports analytic potentials.\n");
#elif defined(KIM)
  if (format != POTENTIAL_FORMAT_KIM)
    error(1, "This potfit binary only supports KIM potentials.\n");
#else
  if (format != POTENTIAL_FORMAT_TABULATED_EQ_DIST && format != POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST)
    error(1, "This potfit binary only supports tabulated potentials.\n");
#endif

  g_pot.format_type = format;

  // only pair potentials for
  // - pair interactions
  // - coulomb interactions
  // - dipole interactions

  int npots = g_calc.paircol;

// more potential functions for other interactions
#if defined(EAM)
  npots += 2 * g_param.ntypes;
#if defined(TBEAM)
  npots += 2 * g_param.ntypes;
#endif  // TBEAM
#endif  // EAM

#if defined(ADP)
  npots = 3 * g_calc.paircol + 2 * g_param.ntypes;
#endif  // ADP

#if defined(MEAM)
  npots = 2 * g_calc.paircol + 3 * g_param.ntypes;
#endif  // MEAM

#if defined(STIWEB)
  npots = 2 * g_calc.paircol + 1;
#endif  // STIWEB

#if defined(TERSOFF) && !defined(TERSOFFMOD)
  npots = g_param.ntypes * g_param.ntypes;
#endif  // TERSOFF && !TERSOFFMOD

#if defined(KIM)
  npots = 1;
#endif // KIM

  if (pstate->num_pots == npots) {
    printf(" - Using %d %s potentials to calculate forces\n", npots,
           g_pot.interaction_name);
    fflush(stdout);
  } else {
#if defined(KIM)
    warning("The number of potentials should always be '1' when KIM is used.\n"
            "You specified %d in file '%s', and it is reset to %d.\n", pstate->num_pots, pstate->filename, npots);
    pstate->num_pots = npots;
#else
    error(0, "Wrong number of data columns in %s potential file \"%s\".\n",
          g_pot.interaction_name, pstate->filename);
    error(1, "For g_param.ntypes=%d there should be %d, but there are %d.\n",
          g_param.ntypes, npots, pstate->num_pots);
#endif // KIM
  }

  g_pot.gradient = (int*)Malloc(npots * sizeof(int));
  g_pot.invar_pot = (int*)Malloc(npots * sizeof(int));
#if defined(APOT)
  g_pot.smooth_pot = (int*)Malloc(npots * sizeof(int));
#endif  // APOT

  pstate->have_format = 1;
}

/****************************************************************
  read_pot_line_T
****************************************************************/

void read_pot_line_T(char const* pbuf, potential_state* pstate)
{
  char* pchar = strchr(pbuf + 3, '\n');

  if (pchar != NULL)
    *pchar = '\0';

  if (strcmp(pbuf + 3, g_pot.interaction_name) != 0) {
    error(0, "Wrong potential type found in potential file!\n");
    error(0, "This binary only supports %s potentials.\n",
          g_pot.interaction_name);
    error(1, "Your potential file contains a %s potential.\n", pbuf + 3);
  }
}

/****************************************************************
  read_pot_line_I
****************************************************************/

void read_pot_line_I(char* pbuf, potential_state* pstate)
{
  int have_invar = 0;
  char* str = NULL;

  if (pstate->have_format) {
#if defined(APOT)
    g_pot.apot_table.invar_pots = 0;
#endif  // APOT

    for (int i = 0; i < pstate->num_pots; i++) {
      str = strtok(((i == 0) ? pbuf + 2 : NULL), " \t\r\n");

      if (str == NULL) {
        error(1, "Not enough items in #I header line.\n");
      } else {
        g_pot.invar_pot[i] = atoi(str);
        if (g_pot.invar_pot[i] == 1) {
          have_invar = 1;
#if defined(APOT)
          g_pot.apot_table.invar_pots++;
#endif  // APOT
        }
      }
    }

    g_pot.have_invar = have_invar;
  } else
    error(1, "#I needs to be specified after #F in file %s\n", pstate->filename);
}

/****************************************************************
  read_pot_line_G
****************************************************************/

void read_pot_line_G(char* pbuf, potential_state* pstate)
{
#if defined(APOT)
  warning("#G header line is not supported for analytic potentials\n");
#else

  if (pstate->have_format) {
    // gradient complete
    for (int i = 0; i < pstate->num_pots; i++) {
      char* str = strtok(((i == 0) ? pbuf + 2 : NULL), " \t\r\n");

      if (str == NULL)
        error(1, "Not enough items in #G header line.\n");
      else
        g_pot.gradient[i] = atoi(str);
    }
    pstate->have_gradient = 1;
  } else
    error(1, "#G needs to be specified after #F in file %s\n", pstate->filename);
#endif  // !APOT
}

/****************************************************************
  read_pot_line_C
****************************************************************/

void read_pot_line_C(char* pbuf, potential_state* pstate)
{
  char names[g_param.ntypes][5];

  // check if there are enough items
  if (pstate->have_format) {
    for (int i = 0; i < g_param.ntypes; i++) {
      char* str = strtok(((i == 0) ? pbuf + 2 : NULL), " \t\r\n");
      if (str == NULL)
        error(1, "Not enough items in #C header line in file %s.\n",
              pstate->filename);
      int len = max(strlen(str), 4);
      strncpy(names[i], str, len);
      names[i][len] = '\0';
    }
  } else
    error(1, "#C needs to be specified after #F in file %s\n", pstate->filename);

  g_config.elements = (char const**)Malloc(g_param.ntypes * sizeof(char*));
  for (int i = 0; i < g_param.ntypes; i++) {
    g_config.elements[i] = (char*)Malloc((strlen(names[i]) + 1) * sizeof(char));
    strncpy((char*)g_config.elements[i], names[i], strlen(names[i]));
    *((char*)g_config.elements[i] + strlen(names[i])) = '\0';
  }
}

/****************************************************************
  allocate_memory_for_potentials
****************************************************************/

void allocate_memory_for_potentials(potential_state* pstate)
{
  pot_table_t* pt = &g_pot.opt_pot;
  const int size = pstate->num_pots;

  // allocate info block of function table
  pt->len = 0;
  pt->ncols = size;
  pt->begin = (double*)Malloc(size * sizeof(double));
  pt->end = (double*)Malloc(size * sizeof(double));
  pt->step = (double*)Malloc(size * sizeof(double));
  pt->invstep = (double*)Malloc(size * sizeof(double));
  pt->first = (int*)Malloc(size * sizeof(int));
  pt->last = (int*)Malloc(size * sizeof(int));

#if defined(APOT)
  apot_table_t* apt = &g_pot.apot_table;

  // allocate memory for analytic potential table
  apt->number = size;
  apt->total_par = 0;

  apt->n_par = (int*)Malloc(size * sizeof(int));
  apt->begin = (double*)Malloc(size * sizeof(double));
  apt->end = (double*)Malloc(size * sizeof(double));
  apt->param_name = (char***)Malloc(size * sizeof(char**));
  apt->fvalue = (fvalue_pointer*)Malloc(size * sizeof(fvalue_pointer));

#if !defined(COULOMB)

#if defined(PAIR)
  if (g_param.enable_cp) {
    apt->values = (double**)Malloc((size + 1) * sizeof(double*));
    apt->values[size] = (double*)Malloc(g_param.ntypes * sizeof(double));
    apt->invar_par = (int**)Malloc(size * sizeof(int*));
    apt->chempot = apt->values[size];
    apt->pmin = (double**)Malloc((size + 1) * sizeof(double*));
    apt->pmin[size] = (double*)Malloc(g_param.ntypes * sizeof(double));
    apt->pmax = (double**)Malloc((size + 1) * sizeof(double*));
    apt->pmax[size] = (double*)Malloc(g_param.ntypes * sizeof(double));
  } else
#endif  // PAIR
  {
    apt->values = (double**)Malloc(size * sizeof(double*));
    apt->invar_par = (int**)Malloc(size * sizeof(int*));
    apt->pmin = (double**)Malloc(size * sizeof(double*));
    apt->pmax = (double**)Malloc(size * sizeof(double*));
  }

#else  // !COULOMB
  apt->ratio = (double*)Malloc(g_param.ntypes * sizeof(double));
  apt->values = (double**)Malloc((size + 5) * sizeof(double*));
  apt->param_name = (char***)Malloc((size + 5) * sizeof(char**));
  apt->pmin = (double**)Malloc((size + 5) * sizeof(double*));
  apt->pmax = (double**)Malloc((size + 5) * sizeof(double*));
  apt->invar_par = (int**)Malloc((size + 5) * sizeof(int*));

  apt->values[size] = (double*)Malloc((g_param.ntypes - 1) * sizeof(double));
  apt->param_name[size] = (char**)Malloc((g_param.ntypes - 1) * sizeof(char*));
  apt->pmin[size] = (double*)Malloc((g_param.ntypes - 1) * sizeof(double));
  apt->pmax[size] = (double*)Malloc((g_param.ntypes - 1) * sizeof(double));
  apt->invar_par[size] = (int*)Malloc((g_param.ntypes - 1) * sizeof(int));

  apt->values[size + 1] = (double*)Malloc(sizeof(double));
  apt->param_name[size + 1] = (char**)Malloc(sizeof(char*));
  apt->pmin[size + 1] = (double*)Malloc(sizeof(double));
  apt->pmax[size + 1] = (double*)Malloc(sizeof(double));
  apt->invar_par[size + 1] = (int*)Malloc(sizeof(int));

  apt->values[size + 2] = (double*)Malloc(g_param.ntypes * sizeof(double));
  apt->param_name[size + 2] = (char**)Malloc(g_param.ntypes * sizeof(char*));
  apt->pmin[size + 2] = (double*)Malloc(g_param.ntypes * sizeof(double));
  apt->pmax[size + 2] = (double*)Malloc(g_param.ntypes * sizeof(double));
  apt->invar_par[size + 2] = (int*)Malloc(g_param.ntypes * sizeof(int));

  for (int i = 3; i < 5; i++) {
    apt->values[size + i] = (double*)Malloc(g_calc.paircol * sizeof(double));
    apt->param_name[size + i] = (char**)Malloc(g_calc.paircol * sizeof(char*));
    apt->pmin[size + i] = (double*)Malloc(g_calc.paircol * sizeof(double));
    apt->pmax[size + i] = (double*)Malloc(g_calc.paircol * sizeof(double));
    apt->invar_par[size + i] = (int*)Malloc(g_calc.paircol * sizeof(int));
  }

  apt->charge = apt->values[size];
  apt->dp_kappa = apt->values[size + 1];
#if defined(DIPOLE)
  apt->dp_alpha = apt->values[size + 2];
  apt->dp_b = apt->values[size + 3];
  apt->dp_c = apt->values[size + 4];
#endif  // DIPOLE
#endif  // !COULOMB

  apt->names = (char**)Malloc(size * sizeof(char*));
  for (int i = 0; i < size; i++)
    apt->names[i] = (char*)Malloc(20 * sizeof(char));
#endif  // APOT
}

/****************************************************************
  calculate_cutoffs
****************************************************************/

void calculate_cutoffs()
{
  const int n = g_param.ntypes;

  pot_table_t* pt = &g_pot.opt_pot;

  g_config.rmin = (double*)Malloc(n * n * sizeof(double));
  g_config.rcut = (double*)Malloc(n * n * sizeof(double));

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      int k = (i <= j) ? i * n + j - ((i * (i + 1)) / 2)
                       : j * n + i - ((j * (j + 1)) / 2);
      g_config.rmin[i * n + j] = pt->begin[k];
      g_config.rcut[i * n + j] = pt->end[k];
    }
  }
#if defined(EAM) || defined(ADP) || defined(MEAM)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      g_config.rcut[i * n + j] =
          MAX(g_config.rcut[i * n + j], pt->end[(n * (n + 1)) / 2 + i]);
      g_config.rcut[i * n + j] =
          MAX(g_config.rcut[i * n + j], pt->end[(n * (n + 1)) / 2 + j]);
      g_config.rmin[i * n + j] =
          MAX(g_config.rmin[i * n + j], pt->begin[(n * (n + 1)) / 2 + i]);
      g_config.rmin[i * n + j] =
          MAX(g_config.rmin[i * n + j], pt->begin[(n * (n + 1)) / 2 + j]);
    }
  }
#endif  // EAM || ADP || MEAM

#if defined(COULOMB)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      g_config.rcut[i * n + j] =
          MAX(g_config.rcut[i * n + j], g_config.dp_cut);
    }
  }
#endif  // COULOMB

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      g_config.rcutmin = MIN(g_config.rcutmin, g_config.rcut[i + n * j]);
      g_config.rcutmax = MAX(g_config.rcutmax, g_config.rcut[i + n * j]);
    }
  }
}

/****************************************************************
  read_maxch_file
****************************************************************/

void read_maxch_file()
{
#if !defined(APOT)
  if (g_param.usemaxch) {
    FILE* pfile = fopen(g_files.maxchfile, "r");

    if (pfile == NULL)
      error(1, "Could not open file %s\n", g_files.maxchfile);

    g_calc.maxchange = (double*)Malloc(g_pot.opt_pot.len * sizeof(double));

    // read maximal changes file
    double* val = g_calc.maxchange;

    for (int i = 0; i < g_pot.opt_pot.len; ++i) {
      if (1 > fscanf(pfile, " %lf\n", val))
        error(1, "Premature end of maxch file %s\n", g_files.maxchfile);
      else
        val++;
    }

    fclose(pfile);
  }
#endif  // !APOT
}

#if defined(APOT)

/****************************************************************
  update_apot_table
    update g_pot.apot_table from g_pot.opt_pot_table, including globals
****************************************************************/

void update_apot_table(double* xi)
{
  for (int i = 0; i < g_calc.ndim; i++)
    g_pot.apot_table
        .values[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]] =
        xi[g_pot.opt_pot.idx[i]];

  if (g_pot.have_globals) {
    for (int i = 0; i < g_pot.apot_table.globals; i++) {
      for (int j = 0; j < g_pot.apot_table.n_glob[i]; j++) {
        int m = g_pot.apot_table.global_idx[i][j][0];
        int n = g_pot.apot_table.global_idx[i][j][1];
        g_pot.apot_table.values[m][n] = *(xi + g_pot.global_idx + i);
      }
    }
  }
}

/****************************************************************
  update_calc_table
    update g_pot.calc_pot.table from g_pot.opt_pot.table, including globals
****************************************************************/

void update_calc_table(double* xi_opt, double* xi_calc, int do_all)
{
  double h = 0.0;
  double* val = xi_opt;
  double* list = g_pot.calc_list + 2;

  // copy global parameters to the right positions
  if (g_pot.have_globals) {
    for (int i = 0; i < g_pot.apot_table.globals; i++) {
      for (int j = 0; j < g_pot.apot_table.n_glob[i]; j++) {
        int m = g_pot.apot_table.global_idx[i][j][0];
        int n = g_pot.apot_table.global_idx[i][j][1];
        *(val + g_pot.opt_pot.first[m] + n) = *(val + g_pot.global_idx + i);
      }
    }
  }

  for (int i = 0; i < g_pot.calc_pot.ncols; i++) {
    if (g_pot.smooth_pot[i] && (do_all || !g_pot.invar_pot[i])) {
      h = *(val + 1 + g_pot.apot_table.n_par[i]);
      if (h == 0.0)
        error(1, "The cutoff parameter for potential %d is 0!\n", i + 1);
    }

    // gradient at rmin
    *val = apot_gradient(g_pot.calc_pot.begin[i], val + 2,
                         g_pot.apot_table.fvalue[i]);

    val += 2;

    // check if something has changed
    int change = 0;

    for (int j = 0; j < g_pot.apot_table.n_par[i]; j++) {
      if (list[j] != val[j]) {
        change = 1;
        list[j] = val[j];
      }
    }

    if (do_all || (change && !g_pot.invar_pot[i])) {
      for (int j = 0; j < APOT_STEPS; j++) {
        double f = 0;
        int k = i * APOT_STEPS + (i + 1) * 2 + j;

        g_pot.apot_table.fvalue[i](g_pot.calc_pot.xcoord[k], val, &f);

        *(xi_calc + k) = g_pot.smooth_pot[i]
                             ? f * apot_cutoff(g_pot.calc_pot.xcoord[k],
                                               g_pot.apot_table.end[i], h)
                             : f;

        if (isnan(f) || isnan(*(xi_calc + k))) {
#if defined(DEBUG)
          error(0, "Potential value was nan or inf. Aborting.\n");
          error(0, "This occured in potential %d (%s)\n", i,
                g_pot.apot_table.names[i]);
          error(0, "at distance r=%f with the parameters:\n",
                g_pot.calc_pot.xcoord[k]);
          for (int m = 0; m < g_pot.apot_table.n_par[i]; m++)
            error(0, "%s %f\n", g_pot.apot_table.param_name[i][m], *(val + m));
          if (g_pot.smooth_pot[i])
            error(0, "h %f\n", h);
#endif  // DEBUG
          error(1, "Potential value was nan or inf. Aborting.\n");
        }
      }
    }
    val += g_pot.apot_table.n_par[i];
    list += g_pot.apot_table.n_par[i] + 2;
  }
}

#endif  // APOT
