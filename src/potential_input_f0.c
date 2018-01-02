/****************************************************************
 *
 * potential_input_f0.c: Routines for reading a potential table
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

#include "chempot.h"
#include "functions.h"
#include "memory.h"
#include "potential_input.h"
#include "utils.h"

#if !defined(APOT)

void read_pot_table0(char const* potential_filename, FILE* pfile)
{
  error(1, "Unsupported potential format in %s", potential_filename);
}

#else

typedef struct {
  char const* filename;
  FILE* pfile;
  fpos_t startpos;
} apot_state;

void read_chemical_potentials(apot_state* pstate);
void read_elstat_table(apot_state* pstate);
void read_global_parameters(apot_state* pstate);
void read_analytic_potentials(apot_state* pstate);

void init_calc_table0();

/****************************************************************
 *
 *  read potential in analytic format:
 *    for more information an how to specify an analytic potential
 *    please check the documentation
 *
 *  parameters:
 *    pot_table_t * ... pointer to the potential table
 *    apot_table_t * ... pointer to the analytic potential table
 *    char * ... name of the potential file (for error messages)
 *    FILE * ... open file handle of the potential file
 *
 ****************************************************************/

void read_pot_table0(char const* potential_filename, FILE* pfile)
{
  apot_state state;
  apot_table_t* apt = &g_pot.apot_table;
  pot_table_t* pt = &g_pot.opt_pot;

  state.filename = potential_filename;
  state.pfile = pfile;

  // save starting position
  fgetpos(pfile, &state.startpos);

  // initialize the function table for analytic potentials
  initialize_analytic_potentials();

  read_chemical_potentials(&state);

  read_elstat_table(&state);

  read_global_parameters(&state);

  read_analytic_potentials(&state);

#if defined(COULOMB)
  apt->total_ne_par = apt->total_par;
#endif  // COULOMB

  /* if we have global parameters, are they actually used ? */
  if (g_pot.have_globals) {
    int use_count = 0;
    for (int i = 0; i < apt->globals; i++)
      use_count += apt->n_glob[i];
    if (use_count == 0) {
      g_pot.have_globals = 0;
      printf("You defined global parameters but did not use them.\n");
      printf("Disabling global parameters.\n\n");
    }
  }

  /* assign the potential functions to the function pointers */
  if (apot_assign_function_pointers(apt) == -1)
    error(1, "Could not assign the function pointers.\n");
#if defined(PAIR)
  if (g_param.enable_cp) {
    g_pot.cp_start =
        apt->total_par - apt->globals + g_param.ntypes * (g_param.ntypes + 1);
    apt->total_par += (g_param.ntypes + g_param.compnodes -
                       apt->invar_par[apt->number][g_param.ntypes]);
  }
#endif  // PAIR

#if defined(COULOMB)
  apt->total_par += g_param.ntypes;
#endif  // COULOMB

#if defined(DIPOLE)
  apt->total_par += g_param.ntypes * (g_param.ntypes + 2);
#endif  // DIPOLE

  /* initialize function table and write indirect index */
  for (int i = 0; i < apt->number; i++) {
    pt->begin[i] = apt->begin[i];
    pt->end[i] = apt->end[i];
    pt->step[i] = 0;
    pt->invstep[i] = 0;
    if (i == 0)
      pt->first[i] = 2;
    else
      pt->first[i] = pt->last[i - 1] + 3;
    pt->last[i] = pt->first[i] + apt->n_par[i] - 1;
  }
  pt->len = pt->first[apt->number - 1] + apt->n_par[apt->number - 1];
  if (g_pot.have_globals)
    pt->len += apt->globals;

#if defined(PAIR)
  if (g_param.enable_cp) {
    pt->len += (g_param.ntypes + g_param.compnodes);
  }
#endif  // PAIR

#if defined(COULOMB)
  pt->len += 2 * g_param.ntypes - 1;
#endif  // COULOMB
#if defined(DIPOLE)
  pt->len += g_param.ntypes * (g_param.ntypes + 2);
#endif  // DIPOLE

  pt->table = (double*)Malloc(pt->len * sizeof(double));

  g_pot.calc_list = (double*)Malloc(pt->len * sizeof(double));

  pt->idx = (int*)Malloc(pt->len * sizeof(int));

  apt->idxpot = (int*)Malloc(apt->total_par * sizeof(int));

  apt->idxparam = (int*)Malloc(apt->total_par * sizeof(int));

  /* this is the indirect index */
  int k = 0;
  int l = 0;
  double* val = pt->table;
  double* list = g_pot.calc_list;

  for (int i = 0; i < apt->number; i++) { /* loop over potentials */
    val += 2;
    list += 2;
    l += 2;
    for (int j = 0; j < apt->n_par[i]; j++) { /* loop over parameters */
      *val = apt->values[i][j];
      *list = apt->values[i][j];
      val++;
      list++;
      if (!g_pot.invar_pot[i] && !apt->invar_par[i][j]) {
        pt->idx[k] = l++;
        apt->idxpot[k] = i;
        apt->idxparam[k++] = j;
      } else
        l++;
    }
    if (!g_pot.invar_pot[i])
      pt->idxlen += apt->n_par[i] - apt->invar_par[i][apt->n_par[i]];
    apt->total_par -= apt->invar_par[i][apt->n_par[i]];
  }

#if defined(PAIR) && !defined(KIM)
  if (g_param.enable_cp) {
    init_chemical_potential(g_param.ntypes);
    int i = apt->number;
    for (int j = 0; j < (g_param.ntypes + g_param.compnodes); j++) {
      *val = apt->values[i][j];
      val++;
      if (!apt->invar_par[i][j]) {
        pt->idx[k] = l++;
        apt->idxpot[k] = i;
        apt->idxparam[k++] = j;
      }
    }
    pt->idxlen += (g_param.ntypes + g_param.compnodes -
                   apt->invar_par[apt->number][g_param.ntypes]);
    g_pot.global_idx += (g_param.ntypes + g_param.compnodes -
                         apt->invar_par[apt->number][g_param.ntypes]);
  }
#endif  // PAIR

#if defined(COULOMB)
  int i = apt->number;
  for (int j = 0; j < (g_param.ntypes - 1); j++) {
    *val = apt->values[i][j];
    val++;
    if (!apt->invar_par[i][j]) {
      pt->idx[k] = l++;
      apt->idxpot[k] = i;
      apt->idxparam[k++] = j;
    } else {
      l++;
      apt->total_par -= apt->invar_par[i][j];
      pt->idxlen -= apt->invar_par[i][j];
    }
  }
  i = apt->number + 1;
  *val = apt->values[i][0];
  val++;
  if (!apt->invar_par[i][0]) {
    pt->idx[k] = l++;
    apt->idxpot[k] = i;
    apt->idxparam[k++] = 0;
  } else {
    l++;
    apt->total_par -= apt->invar_par[i][0];
    pt->idxlen -= apt->invar_par[i][0];
  }
  pt->idxlen += g_param.ntypes;
#endif  // COULOMB

#if defined(DIPOLE)
  i = apt->number + 2;
  for (int j = 0; j < (g_param.ntypes); j++) {
    *val = apt->values[i][j];
    val++;
    if (!apt->invar_par[i][j]) {
      pt->idx[k] = l++;
      apt->idxpot[k] = i;
      apt->idxparam[k++] = j;
    } else {
      l++;
      apt->total_par -= apt->invar_par[i][j];
      pt->idxlen -= apt->invar_par[i][j];
    }
  }
  for (i = apt->number + 3; i < apt->number + 5; i++) {
    for (int j = 0; j < (g_param.ntypes * (g_param.ntypes + 1) / 2); j++) {
      *val = apt->values[i][j];
      val++;
      if (!apt->invar_par[i][j]) {
        pt->idx[k] = l++;
        apt->idxpot[k] = i;
        apt->idxparam[k++] = j;
      } else {
        l++;
        apt->total_par -= apt->invar_par[i][j];
        pt->idxlen -= apt->invar_par[i][j];
      }
    }
  }
  pt->idxlen += g_param.ntypes * (g_param.ntypes + 2);
#endif  // DIPOLE

  if (g_pot.have_globals) {
    int i = g_pot.global_pot;
    for (int j = 0; j < apt->globals; j++) {
      *val = apt->values[i][j];
      *list = apt->values[i][j];
      val++;
      list++;
      if (!apt->invar_par[i][j]) {
        pt->idx[k] = l++;
        apt->idxpot[k] = i;
        apt->idxparam[k++] = j;
      } else
        l++;
    }
    pt->idxlen += apt->globals - apt->invar_par[i][apt->globals];
    apt->total_par -= apt->invar_par[i][apt->globals];
  }
  g_pot.global_idx += pt->last[apt->number - 1] + 1;

#if defined(NOPUNISH)
  if (g_param.opt)
    warning("Gauge degrees of freedom are NOT fixed!\n");
#endif  // NOPUNISH

  check_correct_apot_functions();

  init_calc_table0();

  return;
}

/****************************************************************
 *
 *  read_chemical_potentials:
 *
 ****************************************************************/

void read_chemical_potentials(apot_state* pstate)
{
#if defined(PAIR)
  apot_table_t* apt = &g_pot.apot_table;

  char buffer[255];
  char name[255] = {0};

  fpos_t filepos;

  if (g_param.enable_cp) {
    /* search for cp */
    do {
      fgetpos(pstate->pfile, &filepos);
      if (1 != fscanf(pstate->pfile, "%s", buffer))
        error(1, "Error while searching for chemical potentials\n");
    } while (strncmp(buffer, "cp", 2) != 0 && !feof(pstate->pfile));

    /* and save the position */
    fsetpos(pstate->pfile, &filepos);

    /* shortcut for apt->number */
    int i = apt->number;

    /* allocate memory for global parameters */
    apt->names = (char**)Realloc(apt->names, (i + 1) * sizeof(char*));
    apt->names[i] = (char*)Malloc(20 * sizeof(char));
    strcpy(apt->names[i], "chemical potentials");

    apt->invar_par = (int**)Realloc(apt->invar_par, (i + 1) * sizeof(int*));
    apt->invar_par[i] = (int*)Malloc((g_param.ntypes + 1) * sizeof(int));

    apt->param_name =
        (char***)Realloc(apt->param_name, (i + 1) * sizeof(char**));
    apt->param_name[i] = (char**)Malloc(g_param.ntypes * sizeof(char*));

    /* loop over all atom types */
    for (int j = 0; j < g_param.ntypes; j++) {
      /* allocate memory for parameter name */
      apt->param_name[i][j] = (char*)Malloc(30 * sizeof(char));

      /* read one line */
      if (4 > fscanf(pstate->pfile, "%s %lf %lf %lf", buffer, &apt->chempot[j],
                     &apt->pmin[i][j], &apt->pmax[i][j]))
        error(1, "Could not read chemical potential for %d. atomtype.\n", j);

      /* split cp and _# */
      char* token = strchr(buffer, '_');

      if (token != NULL) {
        strncpy(name, buffer, strlen(buffer) - strlen(token));
        name[strlen(buffer) - strlen(token)] = '\0';
      }

      if (strcmp("cp", name) != 0) {
        if (strlen(name))
          fprintf(stderr, "Found \"%s\" instead of \"cp\"\n", name);
        error(1, "No chemical potentials found in %s.\n", pstate->filename);
      }

      /* check for invariance and proper value (respect boundaries) */
      apt->invar_par[i][j] = 0.0;
      /* parameter will not be optimized if min==max */
      if (apt->pmin[i][j] == apt->pmax[i][j]) {
        apt->invar_par[i][j] = 1;
        apt->invar_par[i][g_param.ntypes]++;
        /* swap min and max if max<min */
      } else if (apt->pmin[i][j] > apt->pmax[i][j]) {
        double temp = apt->pmin[i][j];
        apt->pmin[i][j] = apt->pmax[i][j];
        apt->pmax[i][j] = temp;
        /* reset value if >max or <min */
      } else if ((apt->values[i][j] < apt->pmin[i][j]) ||
                 (apt->values[i][j] > apt->pmax[i][j])) {
        /* Only print warning if we are optimizing */
        if (g_param.opt) {
          if (apt->values[i][j] < apt->pmin[i][j])
            apt->values[i][j] = apt->pmin[i][j];
          if (apt->values[i][j] > apt->pmax[i][j])
            apt->values[i][j] = apt->pmax[i][j];
          warning("Starting value for chemical potential #%d is ", j + 1);
          warning("outside of specified adjustment range.\n");
          warning("Resetting it to %f.\n", j + 1, apt->values[i][j]);
          if (apt->values[i][j] == 0)
            warning("New value is 0 ! Please be careful about this.\n");
        }
      }
      strcpy(apt->param_name[i][j], buffer);
    }
    printf(" - Enabled %d chemical potential(s)\n", g_param.ntypes);

/* disable composition nodes for now */
#if defined(CN)
    /* read composition nodes */
    if (2 > fscanf(infile, "%s %d", buffer, &compnodes)) {
      if (strcmp("type", buffer) == 0)
        compnodes = -1;
      else
        error(1,
              "Could not read number of composition nodes from potential "
              "file.\n");
    }
    if (strcmp(buffer, "cn") != 0 && g_param.ntypes > 1 &&
        g_param.compnodes != -1)
      error(1, "No composition nodes found in %s.\nUse \"cn 0\" for none.\n",
            pstate->filename);
    if (g_param.ntypes == 1) {
      compnodes = 0;
    }
    if (compnodes != -1) {
      apt->values[apt->number] = (double*)Realloc(
          apt->values[apt->number],
          (g_param.ntypes + g_param.compnodes) * sizeof(double));
      apt->pmin[apt->number] = (double*)Realloc(
          apt->pmin[apt->number],
          (g_param.ntypes + g_param.compnodes) * sizeof(double));
      apt->pmax[apt->number] = (double*)Realloc(
          apt->pmax[apt->number],
          (g_param.ntypes + g_param.compnodes) * sizeof(double));
      apt->chempot = apt->values[apt->number];
      compnodelist = (double*)Malloc((g_param.ntypes + g_param.compnodes) *
                                     sizeof(double));

      for (j = 0; j < compnodes; j++) {
        if (4 > fscanf(infile, "%lf %lf %lf %lf", &compnodelist[j],
                       &apt->chempot[ntypes + j],
                       &apt->pmin[apt->number][ntypes + j],
                       &apt->pmax[apt->number][ntypes + j]))
          error(1, "Could not read composition node %d\n", j + 1);
        if (apt->pmin[apt->number][ntypes + j] > apt->chempot[ntypes + j] ||
            apt->pmax[apt->number][ntypes + j] < apt->chempot[ntypes + j])
          error(1, "composition node %d is out of bounds.\n", j + 1);
      }

      /* check compnodes for valid values */
      if (g_param.ntypes == 2) {
        for (j = 0; j < compnodes; j++)
          if (compnodelist[j] > 1 || compnodelist[j] < 0)
            error(1, "Composition node %d is %f but should be inside [0,1].\n",
                  j + 1, compnodelist[j]);
      }
    }
    if (compnodes != -1)
      printf("Enabled chemical potentials with %d extra composition node(s).\n",
             compnodes);
    if (compnodes == -1)
      compnodes = 0;
#endif  // CN
  }
#endif  // PAIR
}

/****************************************************************
 *
 *  read_elstat_table:
 *      bla bla
 *
 ****************************************************************/

void read_elstat_table(apot_state* pstate)
{
#if defined(COULOMB)
  char buffer[255];
  fpos_t filepos;

  fsetpos(pstate->pfile, &pstate->startpos);
  /* skip to electrostatic section */
  do {
    fgetpos(pstate->pfile, &filepos);
    fscanf(pstate->pfile, "%s", buffer);
  } while (strcmp(buffer, "elstat") != 0 && !feof(pstate->pfile));

  /* check for elstat keyword */
  if (strcmp("elstat", buffer) != 0) {
    error(1, "No elstat option found in %s.\n", pstate->filename);
  }

  /* read electrostatic parameters */
  fscanf(pstate->pfile, " %s", buffer);
  if (strcmp("ratio", buffer) != 0) {
    error(1, "Could not read ratio\n");
  }

  apot_table_t* apt = &g_pot.apot_table;

  for (int i = 0; i < g_param.ntypes; i++) {
    if (1 > fscanf(pstate->pfile, "%lf", &apt->ratio[i])) {
      error(1, "Could not read ratio for atomtype #%d\n", i);
    }
  }
  for (int i = 0; i < g_param.ntypes - 1; i++) {
    apt->param_name[apt->number][i] = (char*)Malloc(30 * sizeof(char));
    if (4 > fscanf(pstate->pfile, "%s %lf %lf %lf",
                   apt->param_name[apt->number][i], &apt->charge[i],
                   &apt->pmin[apt->number][i], &apt->pmax[apt->number][i])) {
      error(1, "Could not read charge for atomtype #%d\n", i);
    }
    apt->invar_par[apt->number][i] = 0;
    if (apt->pmin[apt->number][i] == apt->pmax[apt->number][i]) {
      apt->invar_par[apt->number][i]++;
    }
  }
  apt->param_name[apt->number + 1][0] = (char*)Malloc(30 * sizeof(char));
  if (4 > fscanf(pstate->pfile, "%s %lf %lf %lf",
                 apt->param_name[apt->number + 1][0], &apt->dp_kappa[0],
                 &apt->pmin[apt->number + 1][0],
                 &apt->pmax[apt->number + 1][0])) {
    error(1, "Could not read kappa\n");
  }
  apt->invar_par[apt->number + 1][0] = 0;
  if (apt->pmin[apt->number + 1][0] == apt->pmax[apt->number + 1][0]) {
    apt->invar_par[apt->number + 1][0]++;
  }
  apt->sw_kappa = apt->invar_par[apt->number + 1][0];
#if !defined(DIPOLE)
  printf(" - Read elstat table\n");
#endif  // !DIPOLE
#endif  // COULOMB

#if defined(DIPOLE)
  int ncols = g_param.ntypes * (g_param.ntypes + 1) / 2;

  for (int i = 0; i < g_param.ntypes; i++) {
    apt->param_name[apt->number + 2][i] = (char*)Malloc(30 * sizeof(char));
    if (4 > fscanf(pstate->pfile, "%s %lf %lf %lf",
                   apt->param_name[apt->number + 2][i], &apt->dp_alpha[i],
                   &apt->pmin[apt->number + 2][i],
                   &apt->pmax[apt->number + 2][i])) {
      error(1, "Could not read polarisability for atomtype #%d\n", i);
    }
    apt->invar_par[apt->number + 2][i] = 0;
    if (apt->pmin[apt->number + 2][i] == apt->pmax[apt->number + 2][i]) {
      apt->invar_par[apt->number + 2][i]++;
    }
  }
  for (int i = 0; i < ncols; i++) {
    apt->param_name[apt->number + 3][i] = (char*)Malloc(30 * sizeof(char));
    if (4 > fscanf(pstate->pfile, "%s %lf %lf %lf",
                   apt->param_name[apt->number + 3][i], &apt->dp_b[i],
                   &apt->pmin[apt->number + 3][i],
                   &apt->pmax[apt->number + 3][i])) {
      error(1, "Could not read parameter dp_b for potential #%d\n", i);
    }
    apt->invar_par[apt->number + 3][i] = 0;
    if (apt->pmin[apt->number + 3][i] == apt->pmax[apt->number + 3][i]) {
      apt->invar_par[apt->number + 3][i]++;
    }
  }
  for (int i = 0; i < ncols; i++) {
    apt->param_name[apt->number + 4][i] = (char*)Malloc(30 * sizeof(char));
    if (4 > fscanf(pstate->pfile, "%s %lf %lf %lf",
                   apt->param_name[apt->number + 4][i], &apt->dp_c[i],
                   &apt->pmin[apt->number + 4][i],
                   &apt->pmax[apt->number + 4][i])) {
      error(1, "Could not read parameter dp_c for potential #%d\n", i);
    }
    apt->invar_par[apt->number + 4][i] = 0;
    if (apt->pmin[apt->number + 4][i] == apt->pmax[apt->number + 4][i]) {
      apt->invar_par[apt->number + 4][i]++;
    }
  }

  printf(" - Read elstat table\n");
#endif  // DIPOLE
}

/****************************************************************
 *
 *  read_global_parameters:
 *      bla bla
 *
 ****************************************************************/

void read_global_parameters(apot_state* pstate)
{
  apot_table_t* apt = &g_pot.apot_table;

  char buffer[255];

  fpos_t filepos;

  pot_table_t* pt = &g_pot.opt_pot;

  /* skip to global section */
  fsetpos(pstate->pfile, &pstate->startpos);
  do {
    fgetpos(pstate->pfile, &filepos);
    int ret = fscanf(pstate->pfile, "%s", buffer);
    if (feof(pstate->pfile))
      return;
    else if (ret != 1)
      error(1, "Error while searching for global parameters\n");
  } while (strcmp(buffer, "global") != 0 && !feof(pstate->pfile));
  fsetpos(pstate->pfile, &filepos);

  /* check for global keyword */
  if (strcmp(buffer, "global") == 0) {
    if (2 > fscanf(pstate->pfile, "%s %d", buffer, &apt->globals))
      error(1, "Premature end of potential file %s\n", pstate->filename);
    g_pot.have_globals = 1;
    apt->total_par += apt->globals;

    int i = apt->number + g_param.enable_cp;
#if defined(COULOMB)
    i += 5;
#endif
    int j = apt->globals;
    g_pot.global_pot = i;

    /* allocate memory for global parameters */
    apt->names =
        (char**)Realloc(apt->names, (g_pot.global_pot + 1) * sizeof(char*));
    apt->names[g_pot.global_pot] = (char*)Malloc(20 * sizeof(char));
    strcpy(apt->names[g_pot.global_pot], "global parameters");

    apt->n_glob = (int*)Malloc(apt->globals * sizeof(int));

    apt->global_idx = (int***)Malloc(apt->globals * sizeof(int**));

    apt->values = (double**)Realloc(apt->values,
                                    (g_pot.global_pot + 1) * sizeof(double*));
    apt->values[g_pot.global_pot] = (double*)Malloc(j * sizeof(double));

    apt->invar_par =
        (int**)Realloc(apt->invar_par, (g_pot.global_pot + 1) * sizeof(int*));
    apt->invar_par[g_pot.global_pot] = (int*)Malloc((j + 1) * sizeof(int));

    apt->pmin =
        (double**)Realloc(apt->pmin, (g_pot.global_pot + 1) * sizeof(double*));
    apt->pmin[g_pot.global_pot] = (double*)Malloc(j * sizeof(double));

    apt->pmax =
        (double**)Realloc(apt->pmax, (g_pot.global_pot + 1) * sizeof(double*));
    apt->pmax[g_pot.global_pot] = (double*)Malloc(j * sizeof(double));

    apt->param_name = (char***)Realloc(apt->param_name,
                                       (g_pot.global_pot + 1) * sizeof(char**));
    apt->param_name[g_pot.global_pot] = (char**)Malloc(j * sizeof(char*));

    pt->first = (int*)Realloc(pt->first, (g_pot.global_pot + 1) * sizeof(int));

    /* read the global parameters */
    for (j = 0; j < apt->globals; j++) {
      apt->param_name[g_pot.global_pot][j] = (char*)Malloc(30 * sizeof(char));

      int ret_val = fscanf(
          pstate->pfile, "%s %lf %lf %lf", apt->param_name[g_pot.global_pot][j],
          &apt->values[g_pot.global_pot][j], &apt->pmin[g_pot.global_pot][j],
          &apt->pmax[g_pot.global_pot][j]);
      if (4 > ret_val)
        if (strcmp(apt->param_name[g_pot.global_pot][j], "type") == 0) {
          error(0, "Not enough global parameters!\n");
          error(1, "You specified %d parameter(s), but needed are %d.\n", j, apt->globals);
        }

      /* check for duplicate names */
      for (int k = j - 1; k >= 0; k--) {
        if (strcmp(apt->param_name[g_pot.global_pot][j],
                   apt->param_name[g_pot.global_pot][k]) == 0) {
          error(0, "\nFound duplicate global parameter name!\n");
          error(1, "Parameter #%d (%s) is the same as #%d (%s)\n", j + 1,
                apt->param_name[g_pot.global_pot][j], k + 1,
                apt->param_name[g_pot.global_pot][k]);
        }
      }

      apt->n_glob[j] = 0;

      /* check for invariance and proper value (respect boundaries) */
      /* parameter will not be optimized if min==max */
      apt->invar_par[i][j] = 0;

      if (apt->pmin[i][j] == apt->pmax[i][j]) {
        apt->invar_par[i][j] = 1;
        apt->invar_par[i][apt->globals]++;
      } else if (apt->pmin[i][j] > apt->pmax[i][j]) {
        double temp = apt->pmin[i][j];
        apt->pmin[i][j] = apt->pmax[i][j];
        apt->pmax[i][j] = temp;
      } else if ((apt->values[i][j] < apt->pmin[i][j]) ||
                 (apt->values[i][j] > apt->pmax[i][j])) {
        /* Only print warning if we are optimizing */
        if (g_param.opt) {
          if (apt->values[i][j] < apt->pmin[i][j])
            apt->values[i][j] = apt->pmin[i][j];
          if (apt->values[i][j] > apt->pmax[i][j])
            apt->values[i][j] = apt->pmax[i][j];
          warning("Starting value for global parameter #%d is ", j + 1);
          warning("outside of specified adjustment range.\n");
          warning("Resetting it to %f.\n", j + 1, apt->values[i][j]);
          if (apt->values[i][j] == 0)
            warning("New value is 0 ! Please be careful about this.\n");
        }
      }
    }

    printf(" - Read %d global parameter(s)\n", apt->globals);
  }
}

/****************************************************************
  read_analytic_potentials:
****************************************************************/

void read_analytic_potentials(apot_state* pstate)
{
  apot_table_t* apt = &g_pot.apot_table;

  char buffer[255];
  char name[255];

  fpos_t filepos;

  // skip to actual potentials
  fsetpos(pstate->pfile, &pstate->startpos);
  do {
    fgetpos(pstate->pfile, &filepos);
    if (1 != fscanf(pstate->pfile, "%s", buffer))
      error(1, "Error while searching for analytic potentials\n");
  } while (strcmp(buffer, "type") != 0 && !feof(pstate->pfile));
  fsetpos(pstate->pfile, &filepos);

  for (int i = 0; i < apt->number; i++) {
    // scan for next "type" keyword
    do {
      fgetpos(pstate->pfile, &filepos);
      if (1 != fscanf(pstate->pfile, "%s", buffer))
        error(1, "Error while searching for analytic potentials\n");
    } while (strcmp(buffer, "type") != 0 && !feof(pstate->pfile));
    fsetpos(pstate->pfile, &filepos);

    // read type
    if (2 > fscanf(pstate->pfile, "%s %s", buffer, name))
      error(1, "Premature end of potential file %s\n", pstate->filename);
    if (strcmp(buffer, "type") != 0)
      error(1,
            "Unknown keyword in file %s, expected \"type\" but found \"%s\".\n",
            pstate->filename, buffer);

    // split name and _sc
    char* token = strrchr(name, '_');
    if (token != NULL && strcmp(token + 1, "sc") == 0) {
      strncpy(buffer, name, strlen(name) - 3);
      buffer[strlen(name) - 3] = '\0';
      strcpy(name, buffer);
      g_pot.smooth_pot[i] = 1;
    }

    // check if potential is "pohlong" and change it to bjs
    if (strcmp(name, "pohlong") == 0)
      strcpy(name, "bjs");

    if (apot_get_num_parameters(name) == -1)
      error(1,
            "Unknown function type in file %s, please define \"%s\" in "
            "functions.c.\n",
            pstate->filename, name);

    apt->names[i] = (char*)Malloc((strlen(name) + 1) * sizeof(char));
    strncpy(apt->names[i], name, strlen(name) + 1);
    apt->n_par[i] = apot_get_num_parameters(name);

    // add one parameter for cutoff function if _sc is found
    if (g_pot.smooth_pot[i] == 1)
      apt->n_par[i]++;
    apt->total_par += apt->n_par[i];

    // read cutoff
    if (2 > fscanf(pstate->pfile, "%s %lf", buffer, &apt->end[i]))
      error(1, "Could not read cutoff for potential #%d in file %s\n",
            i, pstate->filename);
    if (strcmp(buffer, "cutoff") != 0)
      error(1, "No cutoff found for the %d. potential (%s) after \"type\" in file " "%s.\n",
            i + 1, apt->names[i], pstate->filename);

    // set very small begin, needed for EAM embedding function
    apt->begin[i] = 0.0001;

    // allocate memory for this parameter
    apt->values[i] = (double*)Malloc(apt->n_par[i] * sizeof(double));
    apt->invar_par[i] = (int*)Malloc((apt->n_par[i] + 1) * sizeof(int));
    apt->pmin[i] = (double*)Malloc(apt->n_par[i] * sizeof(double));
    apt->pmax[i] = (double*)Malloc(apt->n_par[i] * sizeof(double));
    apt->param_name[i] = (char**)Malloc(apt->n_par[i] * sizeof(char*));

    char c = 0;

    // check for comments
    do {
      c = fgetc(pstate->pfile);
    } while (c != 10);

    fgetpos(pstate->pfile, &filepos);

    if (NULL == fgets(buffer, 255, pstate->pfile))
      error(1, "Error reading analytic potentials\n");
    while (buffer[0] == '#') {
      fgetpos(pstate->pfile, &filepos);
      if (NULL == fgets(buffer, 255, pstate->pfile))
        error(1, "Error reading analytic potentials\n");
    }
    fsetpos(pstate->pfile, &filepos);

    /* read parameters */
    apt->invar_par[i][apt->n_par[i]] = 0;

    for (int j = 0; j < apt->n_par[i]; j++) {
      apt->param_name[i][j] = (char*)Malloc(30 * sizeof(char));

      strcpy(apt->param_name[i][j], "empty");

      fgetpos(pstate->pfile, &filepos);

      if (NULL == fgets(name, 255, pstate->pfile))
        error(1, "Error reading analytic potentials: not enough parameters for potential %s\n", apt->names[i]);
      while (name[0] == '#' && !feof(pstate->pfile) &&
             (j != apt->n_par[i] - 1)) {
        if (NULL == fgets(name, 255, pstate->pfile))
          error(1, "Error reading analytic potentials\n");
      }

      if ((j != (apt->n_par[i] - 1)) &&
          (feof(pstate->pfile) || name[0] == '\0')) {
        error(0, "Premature end of potential definition or file.\n");
        error(
            1,
            "Probably your potential definition is missing some parameters.\n");
      }

      buffer[0] = '\0';

      int ret_val = sscanf(name, "%s %lf %lf %lf", buffer, &apt->values[i][j],
                           &apt->pmin[i][j], &apt->pmax[i][j]);

      if (strlen(buffer))
        strncpy(apt->param_name[i][j], buffer, 30);

      /* if last char of name is "!" we have a global parameter */
      if (strrchr(apt->param_name[i][j], '!') != NULL) {
        apt->param_name[i][j][strlen(apt->param_name[i][j]) - 1] = '\0';
        int l = -1;
        for (int k = 0; k < apt->globals; k++) {
          if (strcmp(apt->param_name[i][j],
                     apt->param_name[g_pot.global_pot][k]) == 0)
            l = k;
        }

        if (-1 == l)
          error(1, "Could not find global parameter %s!\n",
                apt->param_name[i][j]);

        sprintf(apt->param_name[i][j], "%s!", apt->param_name[i][j]);

        /* write index array for global parameters */
        apt->n_glob[l]++;

        apt->global_idx[l] =
            (int**)Realloc(apt->global_idx[l], apt->n_glob[l] * sizeof(int*));

        apt->global_idx[l][apt->n_glob[l] - 1] = (int*)Malloc(2 * sizeof(int));
        apt->global_idx[l][apt->n_glob[l] - 1][0] = i;
        apt->global_idx[l][apt->n_glob[l] - 1][1] = j;

        apt->values[i][j] = apt->values[g_pot.global_pot][l];
        apt->pmin[i][j] = apt->pmin[g_pot.global_pot][l];
        apt->pmax[i][j] = apt->pmax[g_pot.global_pot][l];
        apt->invar_par[i][j] = 1;
        apt->invar_par[i][apt->n_par[i]]++;
      } else {
        /* this is no global parameter */
        if (4 > ret_val) {
          if (g_pot.smooth_pot[i] &&
              j == apot_get_num_parameters(apt->names[i])) {
            if (0 == strcmp(apt->param_name[i][j], "type") ||
                0 == strcmp(apt->param_name[i][j], "empty") ||
                feof(pstate->pfile)) {
              warning(
                  "No cutoff parameter given for potential #%d: adding one "
                  "parameter.\n",
                  i);
              strcpy(apt->param_name[i][j], "h");
              apt->values[i][j] = 1;
              apt->pmin[i][j] = 0.5;
              apt->pmax[i][j] = 2;
              fsetpos(pstate->pfile, &filepos);
            }
          } else {
            if (strcmp(apt->param_name[i][j], "type") == 0 || ret_val == EOF) {
              error(
                  0,
                  "Not enough parameters for potential #%d (%s) in file %s!\n",
                  i + 1, apt->names[i], pstate->filename);
              error(1, "You specified %d parameter(s), but required are %d.\n",
                    j, apt->n_par[i]);
            }
            error(1, "Could not read parameter #%d of potential #%d in file %s\n",
                  j + 1, i + 1, pstate->filename);
          }
        }

        /* check for invariance and proper value (respect boundaries) */
        /* parameter will not be optimized if min==max */
        apt->invar_par[i][j] = 0;

        if (apt->pmin[i][j] == apt->pmax[i][j]) {
          apt->invar_par[i][j] = 1;
          apt->invar_par[i][apt->n_par[i]]++;
        } else if (apt->pmin[i][j] > apt->pmax[i][j]) {
          double temp = apt->pmin[i][j];
          apt->pmin[i][j] = apt->pmax[i][j];
          apt->pmax[i][j] = temp;
        } else if ((apt->values[i][j] < apt->pmin[i][j]) ||
                   (apt->values[i][j] > apt->pmax[i][j])) {
          /* Only print warning if we are optimizing */
          if (g_param.opt) {
            if (apt->values[i][j] < apt->pmin[i][j])
              apt->values[i][j] = apt->pmin[i][j];
            if (apt->values[i][j] > apt->pmax[i][j])
              apt->values[i][j] = apt->pmax[i][j];
            warning("Starting value for parameter #%d in potential #%d is ",
                    j + 1, i + 1);
            warning("outside of specified adjustment range.\n");
            warning("Resetting it to %f.\n", apt->values[i][j]);
            if (apt->values[i][j] == 0)
              warning("New value is 0 ! Please be careful about this.\n");
          }
        }
      }
    }
  }
  printf(" - Successfully read %d potential table(s)\n", apt->number);
}

/****************************************************************
 *
 *  init_calc_table0: Initialize table used for calculation.
 *
 ****************************************************************/

void init_calc_table0()
{
  const int size = g_pot.apot_table.number;

  pot_table_t* calc = &g_pot.calc_pot;
  pot_table_t* opt = &g_pot.opt_pot;

  calc->len =
      size * APOT_STEPS + 2 * opt->ncols + g_param.ntypes + g_param.compnodes;
  calc->idxlen = APOT_STEPS;
  calc->ncols = opt->ncols;
  calc->begin = opt->begin;
  calc->end = opt->end;
  calc->first = (int*)Malloc(size * sizeof(int));
  calc->last = (int*)Malloc(size * sizeof(int));
  calc->step = (double*)Malloc(size * sizeof(double));
  calc->invstep = (double*)Malloc(size * sizeof(double));
  calc->xcoord = (double*)Malloc(calc->len * sizeof(double));
  calc->table = (double*)Malloc(calc->len * sizeof(double));
  calc->d2tab = (double*)Malloc(calc->len * sizeof(double));
  calc->idx = (int*)Malloc(calc->len * sizeof(int));

  double f = 0;
  int x = 0;

  /* initialize the g_pot.calc_pot table */
  for (int i = 0; i < size; i++) {
    double* val = g_pot.apot_table.values[i];
    double h = g_pot.apot_table.values[i][g_pot.apot_table.n_par[i] - 1];
    calc->table[i * APOT_STEPS + i * 2] = 10e30;
    calc->table[i * APOT_STEPS + i * 2 + 1] = 0;
    calc->first[i] = (x += 2);
    calc->last[i] = (x += APOT_STEPS - 1);
    x++;
    calc->step[i] = (calc->end[i] - calc->begin[i]) / (APOT_STEPS - 1);
    calc->invstep[i] = 1.0 / calc->step[i];

    for (int j = 0; j < APOT_STEPS; j++) {
      int index = i * APOT_STEPS + (i + 1) * 2 + j;
      calc->xcoord[index] = calc->begin[i] + j * calc->step[i];

      g_pot.apot_table.fvalue[i](calc->xcoord[index], val, &f);
      calc->table[index] =
          g_pot.smooth_pot[i]
              ? f * apot_cutoff(calc->xcoord[index], calc->begin[i], h)
              : f;
      calc->idx[i * APOT_STEPS + j] = index;
    }
  }
}

#endif  // APOT
