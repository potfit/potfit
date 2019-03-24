/****************************************************************
 *
 * potential_input_f5.c: Routines for reading a KIM potential table
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

#include <assert.h>

#include "potfit.h"

#include "kim.h"
#include "memory.h"
#include "potential_input.h"

#if !defined(KIM)

void read_pot_table5(char const* potential_filename, FILE* pfile)
{
  error(1, "Unsupported potential format in %s", potential_filename);
}

#else

void print_string_space(const char* string, int len)
{
  printf("%s%*s", string, len - (int)strlen(string), " ");
}

void print_int_space(const int value, int len)
{
  printf("%*d", len, value);
}

void print_description(const char* string, int offset, int desc_len)
{
  const char* substr = string;

  while (*substr) {
    if (strlen(substr) < desc_len) {
      printf("%s\n", substr);
      return;
    }

    // find space before desc_len
    size_t pos = 0;
    for (size_t i = 0; i < strlen(substr) && i < desc_len; ++i) {
      if (substr[i] == ' ')
        pos = i;
    }

    fwrite(substr, sizeof(char), pos, stdout);

    printf("\n%*s", offset, " ");

    substr += pos + 1;
  }
}

/****************************************************************
  print_parameter_info
****************************************************************/

void print_parameter_info()
{
  int name_len = 0;

  for (int i = 0; i < g_kim.nparams; ++i)
    if (strlen(g_kim.params[i].name) > name_len)
      name_len = strlen(g_kim.params[i].name);

  name_len = ((name_len + 7) & ~7);

  const int type_len = 8;
  const int extent_len = 8;
  const int desc_len = 60;

  print_string_space("Name", name_len);
  print_string_space("Type", type_len);
  print_string_space("Extent", extent_len);
  print_string_space("Description", desc_len);
  printf("\n");
  for (int i = 0; i < name_len + type_len + extent_len + desc_len; ++i)
    printf("-");
  printf("\n");

  for (int i = 0; i < g_kim.nparams; ++i) {
    print_string_space(g_kim.params[i].name, name_len);
    print_string_space(KIM_DataType_ToString(g_kim.params[i].type), type_len);
    print_int_space(g_kim.params[i].extent, extent_len - 2);
    printf("  ");
    print_description(g_kim.params[i].desc, name_len + type_len + extent_len, desc_len);
  }
  printf("\n");
}

/****************************************************************
  print_model_parameters_impl
****************************************************************/

void print_model_parameters_impl(FILE* file)
{
  warning("\"kim_model_params\" has been set in the parameter file!\n");
  warning("All possible parameters will be listed and potfit will be terminated!\n\n");
  printf("Available species (%d):\n    ", g_kim.nspecies);

  for (int i = 0; i < g_kim.nspecies; ++i) {
    if ((i + 1) % 13 == 0)
      printf("\n    ");
    printf("%s", KIM_SpeciesName_ToString(g_kim.species[i]));
    if (i != (g_kim.nspecies - 1))
      printf(", ");
  }
  printf("\n\n");

  printf("The following potential parameters are available from the KIM model\n\t%s:\n\n", g_kim.model_name);

  print_parameter_info();

  fprintf(file, "#F 5 1\n");
  fprintf(file, "#T %s\n", g_kim.model_name);
  fprintf(file, "#C");

  for (int i = 0; i < g_kim.nspecies; ++i)
    fprintf(file, " %s", KIM_SpeciesName_ToString(g_kim.species[i]));

  fprintf(file, "\n#E\n\n");

  for (int i = 0; i < g_kim.nparams; ++i) {
    fprintf(file, "KIM_%s %s\n", g_kim.params[i].is_cutoff ? "CUTOFF" : "PARAM", g_kim.params[i].name);
    for (int j = 0; j < g_kim.params[i].extent; ++j)
      if (KIM_DataType_Equal(g_kim.params[i].type, KIM_DATA_TYPE_Integer))
        fprintf(file, "%d %d %d\n", g_kim.params[i].values[j].i, g_kim.params[i].values[j].i, g_kim.params[i].values[j].i);
      else if (KIM_DataType_Equal(g_kim.params[i].type, KIM_DATA_TYPE_Double)) {
        if (g_kim.params[i].is_cutoff)
          fprintf(file, "%f\n", g_kim.params[i].values[j].d);
        else
          fprintf(file, "%f %f %f\n", g_kim.params[i].values[j].d, g_kim.params[i].values[j].d, g_kim.params[i].values[j].d);
      }
    fprintf(file, "\n");
  }

  if (file != stdout)
    printf("A potfit potential file has been written to %s.default\n\n", g_kim.model_name);

  warning("All parameters are fixed in this template by the min = max constraint.\n");
  warning("Please enable those you would like to optimize by changing their min/max values.\n");
}

/****************************************************************
  print_model_parameters
****************************************************************/

void print_model_parameters(KIM_MODEL_PARAMS_ENUM dump)
{
  switch (dump) {
    case KIM_MODEL_PARAMS_NONE:
    case KIM_MODEL_PARAMS_USE_DEFAULT:
      return;
    case KIM_MODEL_PARAMS_DUMP_FILE: {
      size_t len = strlen(g_kim.model_name) + 9;
      char* name = (char*)malloc(len * sizeof(char));
      sprintf(name, "%s.default", g_kim.model_name);
      FILE* pout = fopen(name, "w");
      if (!pout)
        error(1, "Error creating default KIM model potential file");
      print_model_parameters_impl(pout);
      fclose(pout);
      break;
    }
    case KIM_MODEL_PARAMS_DUMP:
      print_model_parameters_impl(stdout);
      break;
  }

  exit_potfit(EXIT_SUCCESS);
}

/****************************************************************
  read_global_parameters:
****************************************************************/

void read_global_parameters(FILE* pfile, const char* filename)
{
  apot_table_t* apt = &g_pot.apot_table;
  char buffer[255];
  fpos_t filepos;
  pot_table_t* pt = &g_pot.opt_pot;

  // skip to global section

  do {
    fgetpos(pfile, &filepos);
    int ret = fscanf(pfile, "%s", buffer);
    if (feof(pfile))
      return;
    else if (ret != 1)
      error(1, "Error while searching for global parameters\n");
  } while (strncmp(buffer, "global", 6) != 0 && !feof(pfile));
  fsetpos(pfile, &filepos);

  // check for global keyword
  if (strncmp(buffer, "global", 6) != 0)
    return;

  if (2 > fscanf(pfile, "%s %d", buffer, &apt->globals))
    error(1, "Premature end of potential file %s\n", filename);

  g_pot.have_globals = 1;
  apt->total_par += apt->globals;
  g_pot.global_pot = apt->number;

  /* allocate memory for global parameters */
  apt->names = (char**)Realloc(apt->names, 2 * sizeof(char*));
  apt->names[1] = (char*)Malloc(18 * sizeof(char));
  snprintf(apt->names[1], 18, "global parameters");

  apt->n_glob = (int*)Malloc(apt->globals * sizeof(int));

  apt->global_idx = (int***)Malloc(apt->globals * sizeof(int**));

  apt->values = (double**)Realloc(apt->values, 2 * sizeof(double*));
  apt->values[1] = (double*)Malloc(apt->globals * sizeof(double));

  apt->invar_par = (int**)Realloc(apt->invar_par, 2 * sizeof(int*));
  apt->invar_par[1] = (int*)Malloc((apt->globals + 1) * sizeof(int));

  apt->pmin = (double**)Realloc(apt->pmin, 2 * sizeof(double*));
  apt->pmin[1] = (double*)Malloc(apt->globals * sizeof(double));

  apt->pmax = (double**)Realloc(apt->pmax, 2 * sizeof(double*));
  apt->pmax[1] = (double*)Malloc(apt->globals * sizeof(double));

  apt->param_name = (const char***)Realloc(apt->param_name, 2 * sizeof(char**));
  apt->param_name[1] = (const char**)Malloc(apt->globals * sizeof(char*));

  pt->first = (int*)Realloc(pt->first, 2 * sizeof(int));

  // read the global parameters
  for (int j = 0; j < apt->globals; j++) {
    int ret_val = fscanf(pfile, "%s %lf %lf %lf", buffer, &apt->values[1][j], &apt->pmin[1][j], &apt->pmax[1][j]);
    if (4 > ret_val) {
      if (strcmp(apt->param_name[1][j], "type") == 0) {
        error(0, "Not enough global parameters!\n");
        error(1, "You specified %d parameter(s), but needed are %d.\n", j, apt->globals);
      }
      error(1, "Error reading global parameters\n");
    }

    if ((buffer[0] > 65) || (buffer[0] > 90 && buffer[0] < 97) || (buffer[0] > 122))
      error(1, "Global parameters for KIM need to start with letters: %s is invalid\n", buffer);
    apt->param_name[1][j] = (const char*)Malloc((strlen(buffer) + 1) * sizeof(char));
    snprintf((char*)apt->param_name[1][j], strlen(buffer) + 1, buffer);

    /* check for duplicate names */
    for (int k = j - 1; k >= 0; k--) {
      if (strcmp(apt->param_name[1][j], apt->param_name[1][k]) == 0) {
        error(0, "\nFound duplicate global parameter name!\n");
        error(1, "Parameter #%d (%s) is the same as #%d (%s)\n", j + 1,
              apt->param_name[1][j], k + 1, apt->param_name[1][k]);
      }
    }

    apt->n_glob[j] = 0;

    // check for invariance and proper value (respect boundaries)
    // parameter will not be optimized if min==max
    apt->invar_par[1][j] = 0;

    if (apt->pmin[1][j] == apt->pmax[1][j]) {
      apt->invar_par[1][j] = 1;
      apt->invar_par[1][apt->globals]++;
    } else if (apt->pmin[1][j] > apt->pmax[1][j]) {
      double temp = apt->pmin[1][j];
      apt->pmin[1][j] = apt->pmax[1][j];
      apt->pmax[1][j] = temp;
    } else if ((apt->values[1][j] < apt->pmin[1][j]) || (apt->values[1][j] > apt->pmax[1][j])) {
      // Only print warning if we are optimizing
      if (g_param.opt) {
        if (apt->values[1][j] < apt->pmin[1][j])
          apt->values[1][j] = apt->pmin[1][j];
        if (apt->values[1][j] > apt->pmax[1][j])
          apt->values[1][j] = apt->pmax[1][j];
        warning("Starting value for global parameter #%d is ", j + 1);
        warning("outside of specified adjustment range.\n");
        warning("Resetting it to %f.\n", j + 1, apt->values[1][j]);
        if (apt->values[1][j] == 0)
          warning("New value is 0! Please be careful about this.\n");
      }
    }
  }

  printf(" - Read %d global parameter(s)\n", apt->globals);
}

/****************************************************************
  read_cutoff_parameter
****************************************************************/

void read_cutoff_parameter(FILE* pfile, const char* filename)
{
  char buffer[255], name[255];

  // skip to cutoff section

  do {
    if (NULL == fgets(buffer, 255, pfile) && !feof(pfile))
      error(1, "Error reading 'KIM_CUTOFF' keyword from potential file %s\n", filename);
    if (strlen(buffer) > 1 && 1 != sscanf(buffer, "%s", name)) {
      error(1, "Error reading 'KIM_CUTOFF ...' keyword from potential file %s\n", filename);
    }
  } while (strncmp(name, "KIM_CUTOFF", 10) != 0 && !feof(pfile));


  // check for global keyword
  if (strncmp(buffer, "KIM_CUTOFF", 10) != 0) {
    warning("No cutoff parameter specified! All parameters will be optimized!\n");
    return;
  }

  int ret = sscanf(buffer, "%*s %s", name);
  if (ret != 1)
    error(1, "Error reading KIM_CUTOFF line!\n");

  // find that parameter in the KIM parameters

  for (int i = 0; i < g_kim.nparams; ++i) {
    if (strcmp(g_kim.params[i].name, name) == 0) {
      g_kim.params[i].is_cutoff = 1;
      return;
    }
  }

  error(1, "Could not associate KIM_CUTOFF name \"%s\" with any parameter!\n", name);
}

void read_keywords_and_parameters(char const* filename, FILE* infile, fpos_t startpos)
{
  print_model_parameters(g_param.kim_model_params);

  fsetpos(infile, &startpos);
  read_global_parameters(infile, filename);

  fsetpos(infile, &startpos);
  if (g_param.kim_model_params == KIM_MODEL_PARAMS_NONE)
    read_cutoff_parameter(infile, filename);
}

double read_max_cutoff()
{
  int num_lists = 0;
  const double* cutoffs = NULL;
  const int* data = NULL; // modelWillNotRequestNeighborsOfNoncontributingParticles
  double max_cutoff = 0.0;

  KIM_Model_GetInfluenceDistance(g_kim.model, &max_cutoff);

  KIM_Model_GetNeighborListPointers(g_kim.model, &num_lists, &cutoffs, &data);

  g_kim.cutoffs = (double*)Malloc((num_lists + 1) * sizeof(double));

  g_kim.cutoffs[num_lists] = max_cutoff;

  for (int i = 0; i < num_lists; ++i) {
    if (!data[i])
      error(1, "KIM model does request neighbors of ghost atoms - NYI in potfit\n");
    if (cutoffs[i] > max_cutoff)
      max_cutoff = cutoffs[i];
    g_kim.cutoffs[i] = cutoffs[i];
  }

  return max_cutoff;
}

void read_kim_parameters_default()
{
  const apot_table_t* apt = &g_pot.apot_table;

  int pos = 0;

  for (int i = 0; i < g_kim.nparams; i++) {
    if (KIM_DataType_Equal(g_kim.params[i].type, KIM_DATA_TYPE_Integer)) {
      for (int j = 0; j < g_kim.params[i].extent; ++j) {
        apt->values[0][pos] = g_kim.params[i].values[j].i;
        apt->pmin[0][pos] = g_kim.params[i].values[j].i;
        apt->pmax[0][pos++] = g_kim.params[i].values[j].i;
      }
    } else {
      for (int j = 0; j < g_kim.params[i].extent; ++j) {
        apt->values[0][pos] = g_kim.params[i].values[j].d;
        apt->pmin[0][pos] = g_kim.params[i].values[j].d;
        apt->pmax[0][pos++] = g_kim.params[i].values[j].d;
      }
    }
  }

  assert(pos == g_kim.total_params);
}

void read_single_kim_parameter(const char* filename, FILE* pfile, int idx, int pos)
{
  const apot_table_t* apt = &g_pot.apot_table;

  char buffer[255];
  char name[255];

  memset(buffer, 0, sizeof(buffer));
  memset(name, 0, sizeof(name));

  do {
    if (NULL == fgets(buffer, 255, pfile) && !feof(pfile))
      error(1, "Error reading 'KIM_' keyword #%d from potential file %s\n", idx + 1, filename);
    if (strlen(buffer) > 1 && 1 != sscanf(buffer, "%s", name))
      error(1, "Error reading 'KIM_' keyword #%d from potential file %s\n", idx + 1, filename);
  } while (strncmp(name, "KIM_", 4) != 0 && !feof(pfile));

  if (strncmp(name, "KIM_", 4) != 0)
    error(1, "Could not find 'KIM_' keyword for parameter \"%s\"\n", g_kim.params[idx].name);

  const int is_cutoff = (strncmp("KIM_CUTOFF", name, 10) == 0);

  int ret = sscanf(buffer, "%*s %s", name);
  if (ret != 1)
    error(1, "Error reading KIM_ line for parameter %d\n", idx + 1);

  // make sure we have the same order as before
  if (strncmp(g_kim.params[idx].name, name, strlen(name)) != 0)
    error(1, "Parameter order mismatch, expected \"%s\" but found \"%s\"\n", g_kim.params[idx].name, name);

  double val = 0.0;
  double min = 0.0;
  double max = 0.0;

  // read all lines for this parameter
  for (int j = 0; j < g_kim.params[idx].extent; ++j) {
    const int offset = pos + j;
    if (NULL == fgets(buffer, 255, pfile))
      error(1, "Error reading '%s' entry #%d from potential file %s\n", name, j + 1, filename);

    if (KIM_DataType_Equal(g_kim.params[idx].type, KIM_DATA_TYPE_Integer)) {
      if (strncmp(buffer, "default", 7) != 0) {
        // read initial value and not optimize it
        if (buffer[0] == '!')
          error(1, "Integer options cannot use global parameters!\n");
        if (!(buffer[0] > 0x2f && buffer[0] < 0x39))
          error(1, "Invalid integer for parameter \"%s\", #%d: \"%s\"!\n", g_kim.params[idx].name, j + 1, buffer);
        ret = sscanf(buffer, "%d", &g_kim.params[idx].values[j].i);
        if (ret != 1)
          error(1, "Error reading default line! TODO\n");
      }
      // store in KIM and forget about it
      ret = KIM_Model_SetParameterInteger(g_kim.model, idx, j, g_kim.params[idx].values[j].i);
      if (ret)
        error(1, "Error setting parameter\n");
      apt->invar_par[0][offset] = 1;
      apt->invar_par[0][apt->n_par[0]] += 1;
    } else if (KIM_DataType_Equal(g_kim.params[idx].type, KIM_DATA_TYPE_Double)) {
      if (is_cutoff) {
        ret = sscanf(buffer, "%lf", &val);
        if (ret != 1)
          error(1, "Error reading parameters line!\n");
        apt->invar_par[0][offset] = 1;
        apt->invar_par[0][apt->n_par[0]] += 1;
        apt->pmin[0][offset] = val;
        apt->pmax[0][offset] = val;
        apt->values[0][offset] = val;
          // store cutoff in KIM
        ret = KIM_Model_SetParameterDouble(g_kim.model, idx, j, val);
        if (ret)
          error(1, "Error setting parameter\n");
        continue;
      }

      // global parameter
      if (buffer[0] == '!') {
          // TODO
          continue;
      }

      if (strncmp(buffer, "default", 7) == 0) {
          // use default value from KIM
          ret = sscanf(buffer, "%*s %lf %lf", &min, &max);
          if (ret != 2)
            error(1, "Error reading default line!\n");
          val = g_kim.params[idx].values[j].d;
      } else {
        ret = sscanf(buffer, "%lf %lf %lf", &val, &min, &max);
        if (ret != 3)
          error(1, "Error reading parameters line!\n");
      }

      if (min > max) {
        double temp = 0.0;
        SWAP(min, max, temp);
      }

      // fixed param or not
      if (min == max) {
        apt->invar_par[0][offset] = 1;
        apt->invar_par[0][apt->n_par[0]] += 1;
      }

      apt->values[0][offset] = val;
      apt->pmin[0][offset] = min;
      apt->pmax[0][offset] = max;

      if(!apt->invar_par[0][offset] && (min > apt->values[0][offset] || max < apt->values[0][offset]))
        error(1, "Value %d of '%s' is not within its limits in file '%s'.\n", j + 1, g_kim.params[idx].name, filename);
    } else {
      error(1, "Unsupported DataType detected!\n");
    }
  }
}

void read_kim_parameters_file(const char* filename, FILE* pfile)
{
  int pos = 0;

  for (int i = 0; i < g_kim.nparams; ++i) {
    read_single_kim_parameter(filename, pfile, i, pos);
    pos += g_kim.params[i].extent;
  }
}

/****************************************************************
 *
 *  read KIM potential with lower and upper limits for each parameter.
 *
 *  parameters:
 *    char * ... name of the potential file (for error messages)
 *    FILE * ... open file handle of the potential file
 *
 ****************************************************************/

void read_pot_table5(char const* filename, FILE* pfile)
{
  apot_table_t* apt = &g_pot.apot_table;
  pot_table_t* pt = &g_pot.opt_pot;

  fpos_t startpos;

  // save starting position
  fgetpos(pfile, &startpos);

  // read the keywords and the names of parameters that will be optimized
  read_keywords_and_parameters(filename, pfile, startpos);

  const double max_cutoff = read_max_cutoff();

  apt->n_par[0] = g_kim.total_params;
  apt->total_par += apt->n_par[0];
  // set begin and end
  apt->begin[0] = 0.0;
  apt->end[0] = max_cutoff;
  // allocate memory
  apt->values[0] = (double *)Malloc(apt->n_par[0] * sizeof(double));
  // last slot stores the total number of invariant parameters
  apt->invar_par[0] = (int *)Malloc((apt->n_par[0] + 1) * sizeof(int));
  apt->pmin[0] = (double *)Malloc(apt->n_par[0] * sizeof(double));
  apt->pmax[0] = (double *)Malloc(apt->n_par[0] * sizeof(double));

  // read parameter values

  fsetpos(pfile, &startpos);

  if (g_param.kim_model_params == KIM_MODEL_PARAMS_USE_DEFAULT)
    read_kim_parameters_default();
  else
    read_kim_parameters_file(filename, pfile);

  printf("Successfully read %d potential parameters\n", apt->total_par);

  // initialize optimization table
  pt->begin[0] = apt->begin[0];
  pt->end[0] = apt->end[0];
  pt->first[0] = 2;
  pt->last[0] = pt->first[0] + apt->n_par[0] - 1;
  pt->len = pt->last[0] - pt->first[0] + 1;

  if (g_pot.have_globals)
    pt->len += apt->globals;

  pt->table = (double *)Malloc(pt->len * sizeof(double));
  pt->idx = (int *)Malloc(pt->len * sizeof(int));
  apt->idxpot = (int *)Malloc(apt->total_par * sizeof(int));
  apt->idxparam = (int *)Malloc(apt->total_par * sizeof(int));

  double* val = pt->table;
  int k = 0;
  int l = 0;

  for (int i = 0; i < apt->n_par[0]; ++i) {
    *val = apt->values[0][i];
    if (!apt->invar_par[0][i]) {
      pt->idx[k] = l;
      apt->idxpot[k] = 0;
      apt->idxparam[k++] = i;
    }
    ++val;
    ++l;
  }

  pt->idxlen += apt->n_par[0] - apt->invar_par[0][apt->n_par[0]];
  apt->total_par -= apt->invar_par[0][apt->n_par[0]];
}

#endif // !KIM
