/****************************************************************
 *
 * potential_input.c: Routines for reading a potential table
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
#include "potential_input.h"

void read_pot_line_F(char const* pbuf, potential_state* pstate);
void read_pot_line_T(char const* pbuf, potential_state* pstate);
void read_pot_line_I(char* pbuf, potential_state* pstate);
void read_pot_line_G(char* pbuf, potential_state* pstate);
void read_pot_line_C(char const* pbuf, potential_state* pstate);

void allocate_memory_for_potentials(potential_state* pstate);

void calculate_cutoffs();
void read_maxch_file();

/****************************************************************
 *
 * read potential tables
 *
 ****************************************************************/

void read_pot_table(char const* potential_filename)
{
  char buffer[1024];
  char* res = NULL;

  FILE* pfile;

  potential_state pstate;

  memset(&pstate, 0, sizeof(pstate));

  pstate.filename = potential_filename;

  /* set paircol to the number of pair potentials */
  g_calc.paircol = (g_param.ntypes * (g_param.ntypes + 1)) / 2;

  /* open file */
  pfile = fopen(potential_filename, "r");

  if (NULL == pfile) error(1, "Could not open file %s\n", potential_filename);

  printf("Starting to read the potential file:\n");

  /* read the header */
  do {
    /* read one line */
    res = fgets(buffer, 1024, pfile);

    if (NULL == res) error(1, "Unexpected end of file in %s", potential_filename);

    /* check if it is a header line */
    if (buffer[0] != '#') error(1, "Header corrupt in file %s", potential_filename);

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
        warning("Ignoring potential header line starting with #%c.", buffer[1]);
    }

  } while (!pstate.end_header);

  /* do we have a format in the header? */
  if (!pstate.have_format)
    error(1, "Format not specified in header of potential file %s", potential_filename);

  allocate_memory_for_potentials(&pstate);

  switch (g_pot.format) {
    case 0:
      read_pot_table0(potential_filename, pfile);
      break;
    case 3:
      read_pot_table3(potential_filename, pfile, &pstate);
      break;
    case 4:
      read_pot_table4(potential_filename, pfile, &pstate);
      break;
  }

  fclose(pfile);

  printf("Reading potential file >> %s << ... done\n", potential_filename);

  calculate_cutoffs();

  read_maxch_file();

  /* clean up locals and mark globals for later */
  //   reg_for_free(gradient, "gradient");
  //   reg_for_free(invar_pot, "invar_pot");
  // #ifdef APOT
  //   reg_for_free(smooth_pot, "smooth_pot");
  //   reg_for_free(apt->n_par, "apt->n_par");
  //   reg_for_free(apt->begin, "apt->begin");
  //   reg_for_free(apt->end, "apt->end");
  //   reg_for_free(apt->param_name, "apt->param_name");
  //   reg_for_free(apt->fvalue, "apt->fvalue");
  //   reg_for_free(apt->values, "apt->values");
  //   reg_for_free(apt->invar_par, "apt->invar_par");
  //   reg_for_free(apt->pmin, "apt->pmin");
  //   reg_for_free(apt->pmax, "apt->pmax");
  //   reg_for_free(apt->names, "apt->names");
  //   for (i = 0; i < size; i++) {
  //     reg_for_free(apt->names[i], "apt->names[%d]", i);
  //   }
  // #else /* APOT */
  //   reg_for_free(maxchange, "maxchange");
  // #endif /* APOT */
  //   reg_for_free(pt->begin, "pt->begin");
  //   reg_for_free(pt->end, "pt->end");
  //   reg_for_free(pt->step, "pt->step");
  //   reg_for_free(pt->invstep, "pt->invstep");
  //   reg_for_free(pt->first, "pt->first");
  //   reg_for_free(pt->last, "pt->last");
  // #if defined PAIR && defined APOT
  //   if (enable_cp) {
  //     reg_for_free(apt->chempot, "apt->chempot");
  //     reg_for_free(apt->pmin[size], "apt->pmin[%d]", size);
  //     reg_for_free(apt->pmax[size], "apt->pmax[%d]", size);
  //   }
  // #endif /* PAIR && APOT */
  // #ifdef COULOMB
  //   reg_for_free(apt->ratio, "apt->ratio");
  //   reg_for_free(apt->charge, "apt->charge");
  //   reg_for_free(apt->dp_kappa, "apt->dp_kappa");
  //   for (i = 0; i < 5; i++) {
  //     reg_for_free(apt->pmin[size + i], "apt->pmin[%d]", size + i);
  //     reg_for_free(apt->pmax[size + i], "apt->pmax[%d]", size + i);
  //     reg_for_free(apt->invar_par[size + i], "apt->invar_par[%d]", size + i);
  //   }
  // #endif /* COULOMB */
  // #ifdef DIPOLE
  //   reg_for_free(apt->dp_alpha, "apt->dp_alpha");
  //   reg_for_free(apt->dp_b, "apt->dp_b");
  //   reg_for_free(apt->dp_c, "apt->dp_c");
  // #endif /* DIPOLE */
  //   reg_for_free(rcut, "rcut");
  //   reg_for_free(rmin, "rmin");
}

/****************************************************************
 *
 *  read_pot_line_F
 *      bla bla
 *
 ****************************************************************/

void read_pot_line_F(char const* pbuf, potential_state* pstate)
{
  /* format complete? */
  if (2 != sscanf((char const*)(pbuf + 2), "%d %d", &g_pot.format, &pstate->num_pots))
    error(1, "Corrupt format header line in file %s", pstate->filename);

#if !defined(APOT)
  if (g_pot.format == 0)
    error(1, "potfit binary compiled without analytic potential support.\n");
#else
  if (g_pot.format > 0)
    error(1, "potfit binary compiled without tabulated potential support.\n");
#endif /* !APOT */

  if (g_pot.format != 0)
    printf(" - Potential file format %d detected\n", g_pot.format);
  else
    printf(" - Potential file format %d (analytic potentials) detected\n", g_pot.format);

  /* only pair potentials for
    * - pair interactions
    * - coulomb interactions
    * - dipole interactions
    */
  int npots = g_calc.paircol;

/* more potential functions for other interactions */
#if defined(EAM)
  npots = g_calc.paircol + 2 * g_param.ntypes;
#if defined(TBEAM)
  npots += 2 * g_param.ntypes;
#endif /* TBEAM */
#endif /* EAM */

#if defined(ADP)
  npots = 3 * g_calc.paircol + 2 * g_param.ntypes;
#endif /* ADP */

#if defined(MEAM)
  npots = 2 * g_calc.paircol + 3 * g_param.ntypes;
#endif /* MEAM */

#if defined(STIWEB)
  npots = 2 * g_calc.paircol + 1;
#endif /* STIWEB */

#if defined(TERSOFF) && !defined(TERSOFFMOD)
  npots = g_param.ntypes * g_param.ntypes;
#endif /* TERSOFF && !TERSOFFMOD */

  if (pstate->num_pots == npots) {
    printf(" - Using %d %s potentials to calculate forces\n", npots,
           g_todo.interaction_name);
    fflush(stdout);
  } else {
    error(0, "Wrong number of data columns in %s potential file \"%s\".\n",
          g_todo.interaction_name, pstate->filename);
    error(1, "For g_param.ntypes=%d there should be %d, but there are %d.",
          g_param.ntypes, npots, pstate->num_pots);
  }
  /* recognized format? */
  if ((g_pot.format != 0) && (g_pot.format != 3) && (g_pot.format != 4))
    error(1, "Unrecognized potential format specified for file %s", pstate->filename);

  g_pot.gradient = (int*)malloc(npots * sizeof(int));
  g_pot.invar_pot = (int*)malloc(npots * sizeof(int));
#if defined(APOT)
  g_pot.smooth_pot = (int*)malloc(npots * sizeof(int));
#endif /* APOT */

  memset(g_pot.gradient, 0, npots * sizeof(int));
  memset(g_pot.invar_pot, 0, npots * sizeof(int));
#if defined(APOT)
  memset(g_pot.smooth_pot, 0, npots * sizeof(int));
#endif /* APOT */

  pstate->have_format = 1;
}

/****************************************************************
 *
 *  read_pot_line_T
 *      bla bla
 *
 ****************************************************************/

void read_pot_line_T(char const* pbuf, potential_state* pstate)
{
  char* pchar = strchr(pbuf + 3, '\n');

  if (pchar != NULL) *pchar = '\0';

  if (strcmp(pbuf + 3, g_todo.interaction_name) != 0) {
    error(0, "Wrong potential type found in potential file!\n");
    error(0, "This binary only supports %s potentials.\n", g_todo.interaction_name);
    error(1, "Your potential file contains a %s potential.\n", pbuf + 3);
  }
}

/****************************************************************
 *
 *  read_pot_line_I
 *      bla bla
 *
 ****************************************************************/

void read_pot_line_I(char* pbuf, potential_state* pstate)
{
  int have_invar = 0;
  int i = 0;
  char* str = NULL;

  if (pstate->have_format) {
#if defined(APOT)
    g_pot.apot_table.invar_pots = 0;
#endif /* APOT */

    for (i = 0; i < pstate->num_pots; i++) {
      str = strtok(((i == 0) ? pbuf + 2 : NULL), " \t\r\n");

      if (str == NULL) {
        error(1, "Not enough items in #I header line.");
      } else {
        g_pot.invar_pot[i] = atoi(str);
        if (g_pot.invar_pot[i] == 1) {
          have_invar = 1;
#if defined(APOT)
          g_pot.apot_table.invar_pots++;
#endif /* APOT */
        }
      }
    }

    g_pot.have_invar = have_invar;
  } else
    error(1, "#I needs to be specified after #F in file %s", pstate->filename);
}

/****************************************************************
 *
 *  read_pot_line_G
 *      bla bla
 *
 ****************************************************************/

void read_pot_line_G(char* pbuf, potential_state* pstate)
{
#if !defined(APOT)
  int i = 0;

  if (pstate->have_format) {
    /* gradient complete */
    for (i = 0; i < pstate->num_pots; i++) {
      char* str = strtok(((i == 0) ? pbuf + 2 : NULL), " \t\r\n");

      if (str == NULL)
        error(1, "Not enough items in #G header line.");
      else
        g_pot.gradient[i] = atoi(str);
    }
    pstate->have_gradient = 1;
  } else
    error(1, "#G needs to be specified after #F in file %s", pstate->filename);
#endif  // !APOT
}

/****************************************************************
 *
 *  read_pot_line_C
 *      bla bla
 *
 ****************************************************************/

void read_pot_line_C(char const* pbuf, potential_state* pstate)
{
  // TODO
}

/****************************************************************
 *
 *  read_pot_line_C
 *      bla bla
 *
 ****************************************************************/

void allocate_memory_for_potentials(potential_state* pstate)
{
  pot_table_t* pt = &g_pot.opt_pot;
  int size = pstate->num_pots;

  /* allocate info block of function table */
  pt->len = 0;
  pt->ncols = size;
  pt->begin = (double*)malloc(size * sizeof(double));
  pt->end = (double*)malloc(size * sizeof(double));
  pt->step = (double*)malloc(size * sizeof(double));
  pt->invstep = (double*)malloc(size * sizeof(double));
  pt->first = (int*)malloc(size * sizeof(int));
  pt->last = (int*)malloc(size * sizeof(int));
  if ((pt->begin == NULL) || (pt->end == NULL) || (pt->step == NULL) ||
      (pt->invstep == NULL) || (pt->first == NULL) || (pt->last == NULL))
    error(1, "Cannot allocate info block for potential table %s", pstate->filename);

#if defined(APOT)
  apot_table_t* apt = &g_pot.apot_table;

  /* allocate memory for analytic potential table */
  apt->number = size;
  apt->total_par = 0;

  apt->n_par = (int*)malloc(size * sizeof(int));
  apt->begin = (double*)malloc(size * sizeof(double));
  apt->end = (double*)malloc(size * sizeof(double));
  apt->param_name = (char***)malloc(size * sizeof(char**));
  apt->fvalue = (fvalue_pointer*)malloc(size * sizeof(fvalue_pointer));

#if defined(PAIR)
  if (g_param.enable_cp) {
    apt->values = (double**)malloc((size + 1) * sizeof(double*));
    apt->values[size] = (double*)malloc(g_param.ntypes * sizeof(double));
    apt->invar_par = (int**)malloc(size * sizeof(int*));
    apt->chempot = apt->values[size];
    apt->pmin = (double**)malloc((size + 1) * sizeof(double*));
    apt->pmin[size] = (double*)malloc(g_param.ntypes * sizeof(double));
    apt->pmax = (double**)malloc((size + 1) * sizeof(double*));
    apt->pmax[size] = (double*)malloc(g_param.ntypes * sizeof(double));
  } else {
    apt->values = (double**)malloc(size * sizeof(double*));
    apt->invar_par = (int**)malloc(size * sizeof(int*));
    apt->pmin = (double**)malloc(size * sizeof(double*));
    apt->pmax = (double**)malloc(size * sizeof(double*));
  }
#endif

#if defined(COULOMB)
  apt->ratio = (double*)malloc(g_param.ntypes * sizeof(double));
  apt->values = (double**)malloc((size + 5) * sizeof(double*));
  apt->param_name = (char***)malloc((size + 5) * sizeof(char**));
  apt->pmin = (double**)malloc((size + 5) * sizeof(double*));
  apt->pmax = (double**)malloc((size + 5) * sizeof(double*));
  apt->invar_par = (int**)malloc((size + 5) * sizeof(int*));

  apt->values[size] = (double*)malloc((g_param.ntypes - 1) * sizeof(double));
  apt->param_name[size] = (char**)malloc((g_param.ntypes - 1) * sizeof(char*));
  apt->pmin[size] = (double*)malloc((g_param.ntypes - 1) * sizeof(double));
  apt->pmax[size] = (double*)malloc((g_param.ntypes - 1) * sizeof(double));
  apt->invar_par[size] = (int*)malloc((g_param.ntypes - 1) * sizeof(int));

  apt->values[size + 1] = (double*)malloc(sizeof(double));
  apt->param_name[size + 1] = (char**)malloc(sizeof(char*));
  apt->pmin[size + 1] = (double*)malloc(sizeof(double));
  apt->pmax[size + 1] = (double*)malloc(sizeof(double));
  apt->invar_par[size + 1] = (int*)malloc(sizeof(int));

  apt->values[size + 2] = (double*)malloc(g_param.ntypes * sizeof(double));
  apt->param_name[size + 2] = (char**)malloc(g_param.ntypes * sizeof(char*));
  apt->pmin[size + 2] = (double*)malloc(g_param.ntypes * sizeof(double));
  apt->pmax[size + 2] = (double*)malloc(g_param.ntypes * sizeof(double));
  apt->invar_par[size + 2] = (int*)malloc(g_param.ntypes * sizeof(int));

  for (int i = 3; i < 5; i++) {
    apt->values[size + i] = (double*)malloc(g_calc.paircol * sizeof(double));
    apt->param_name[size + i] = (char**)malloc(g_calc.paircol * sizeof(char*));
    apt->pmin[size + i] = (double*)malloc(g_calc.paircol * sizeof(double));
    apt->pmax[size + i] = (double*)malloc(g_calc.paircol * sizeof(double));
    apt->invar_par[size + i] = (int*)malloc(g_calc.paircol * sizeof(int));
  }

  apt->charge = apt->values[size];
  apt->dp_kappa = apt->values[size + 1];
#ifdef DIPOLE
  apt->dp_alpha = apt->values[size + 2];
  apt->dp_b = apt->values[size + 3];
  apt->dp_c = apt->values[size + 4];
#endif /* DIPOLE */
#endif /* COULOMB */

  apt->names = (char**)malloc(size * sizeof(char*));

  for (int i = 0; i < size; i++) apt->names[i] = (char*)malloc(20 * sizeof(char));

  if ((apt->n_par == NULL) || (apt->begin == NULL) || (apt->end == NULL) ||
      (apt->fvalue == NULL) || (apt->names == NULL) || (apt->pmin == NULL) ||
      (apt->pmax == NULL) || (apt->param_name == NULL) || (apt->values == NULL))
    error(1, "Cannot allocate info block for analytic potential table %s",
          pstate->filename);
#endif /* APOT */
}

/****************************************************************
 *
 *  calculate_cutoffs
 *      bla bla
 *
 ****************************************************************/

void calculate_cutoffs()
{
  const int n = g_param.ntypes;

  pot_table_t* pt = &g_pot.opt_pot;

  g_config.rmin = (double*)malloc(n * n * sizeof(double));

  if (NULL == g_config.rmin) error(1, "Cannot allocate rmin");

  g_config.rcut = (double*)malloc(n * n * sizeof(double));

  if (NULL == g_config.rcut) error(1, "Cannot allocate rcut");

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      int k =
          (i <= j) ? i * n + j - ((i * (i + 1)) / 2) : j * n + i - ((j * (j + 1)) / 2);
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
#endif /* EAM || ADP || MEAM */

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      g_config.rcutmin = MIN(g_config.rcutmin, g_config.rcut[i + n * j]);
      g_config.rcutmax = MAX(g_config.rcutmax, g_config.rcut[i + n * j]);
    }
  }
}

/****************************************************************
 *
 *  read_maxch_file
 *      bla bla
 *
 ****************************************************************/

void read_maxch_file()
{
#if !defined(APOT)
  if (g_param.usemaxch) {
    int i = 0;

    FILE* pfile = NULL;

    /* open file */
    pfile = fopen(g_files.maxchfile, "r");

    if (pfile == NULL) error(1, "Could not open file %s\n", g_files.maxchfile);

    /* read maximal changes file */
    g_todo.maxchange = (double*)malloc(g_pot.opt_pot.len * sizeof(double));

    double* val = g_todo.maxchange;

    for (i = 0; i < g_pot.opt_pot.len; i++) {
      if (1 > fscanf(pfile, " %lf\n", val))
        error(1, "Premature end of maxch file %s", g_files.maxchfile);
      else
        val++;
    }

    fclose(pfile);
  }
#endif /* !APOT */
}

#ifdef APOT

/****************************************************************
 *
 * update g_pot.apot_table from g_pot.opt_pot_table, including globals
 *
 ****************************************************************/

void update_apot_table(double* xi)
{
  for (int i = 0; i < g_calc.ndim; i++)
    g_pot.apot_table.values[g_pot.apot_table.idxpot[i]][g_pot.apot_table.idxparam[i]] =
        xi[g_todo.idx[i]];
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
 *
 * update g_pot.calc_pot.table from g_pot.opt_pot.table, including globals
 *
 ****************************************************************/

void update_calc_table(double* xi_opt, double* xi_calc, int do_all)
{
  //   int   i, j, k, m, n, change;
  //   double f, h = 0;

  double* val = xi_opt;
  double* list = g_pot.calc_list + 2;

  int h = 0;

  /* copy global parameters to the right positions */
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
      if (h == 0) error(1, "The cutoff parameter for potential %d is 0!", i);
    }

    (*val) = apot_grad(g_pot.calc_pot.begin[i], val + 2, g_pot.apot_table.fvalue[i]);

    val += 2;

    /* check if something has changed */
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

        *(xi_calc + k) =
            g_pot.smooth_pot[i]
                ? f * cutoff(g_pot.calc_pot.xcoord[k], g_pot.apot_table.end[i], h)
                : f;

        if (isnan(f) || isnan(*(xi_calc + k))) {
#if defined(DEBUG)
          error(0, "Potential value was nan or inf. Aborting.\n");
          error(0, "This occured in potential %d (%s)\n", i, g_pot.apot_table.names[i]);
          error(0, "at distance r=%f with the parameters:\n", g_pot.calc_pot.xcoord[k]);
          for (int m = 0; m < g_pot.apot_table.n_par[i]; m++)
            error(0, "%s %f\n", g_pot.apot_table.param_name[i][m], *(val + m));
          if (g_pot.smooth_pot[i]) error(0, "h %f\n", h);
#endif /* DEBUG */
          error(1, "Potential value was nan or inf. Aborting.\n");
        }
      }
    }
    val += g_pot.apot_table.n_par[i];
    list += g_pot.apot_table.n_par[i] + 2;
  }
}

#endif /* APOT */
