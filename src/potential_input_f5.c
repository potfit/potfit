/****************************************************************
 *
 * potential_input_f5.c: Routines for reading a potential table
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

#include <ctype.h>

#include "potfit.h"

#include "chempot.h"
#include "functions.h"
#include "kim.h"
#include "memory.h"
#include "potential_input.h"
#include "utils.h"

#if !defined(APOT) && !defined(KIM)

void read_pot_table5(char const* potential_filename, FILE* pfile)
{
  error(1, "Unsupported potential format in %s", potential_filename);
}

#else

#define KIM_INPUT_FILE_TEMPLATE \
  "KIM_API_Version := 1.7.3\n" \
  "Unit_length      := A\n" \
  "Unit_energy      := eV\n" \
  "Unit_charge      := e\n" \
  "Unit_temperature := K\n" \
  "Unit_time        := ps\n" \
  "PARTICLE_SPECIES:\n" \
  "%s                           spec                0\n" \
  "CONVENTIONS:\n" \
  "ZeroBasedLists               flag\n" \
  "Neigh_LocaAccess             flag\n" \
  "NEIGH_RVEC_H                 flag\n" \
  "NEIGH_RVEC_F                 flag\n" \
  "MI_OPBC_H                    flag\n" \
  "MI_OPBC_F                    flag\n" \
  "MODEL_INPUT:\n" \
  "numberOfParticles            integer   none        []\n" \
  "numberOfSpecies              integer   none        []\n" \
  "particleSpecies              integer   none        [numberOfParticles]\n" \
  "coordinates                  double    length      [numberOfParticles,3]\n" \
  "boxSideLengths               double    length      [3]\n" \
  "numberContributingParticles  integer   none        []\n" \
  "get_neigh                    method    none        []\n" \
  "neighObject                  pointer   none        []\n" \
  "MODEL_OUTPUT:\n" \
  "destroy                      method    none        []\n" \
  "compute                      method    none        []\n" \
  "reinit                       method    none        []           optional\n" \
  "cutoff                       double    length      []\n" \
  "energy                       double    energy      []\n" \
  "forces                       double    force       [numberOfParticles,3]\n" \
  "virial                       double    energy      [6]"

void read_kim_model_params()
{
  void* pkim = NULL;

  // create a temporary KIM object
  int status = KIM_API_model_info(&pkim, g_kim.model_name);
  if (KIM_STATUS_OK > status)
    error(1, "KIM_API_model_info failed");

  // get the first species supported by the model
  char const* species;

  status = KIM_API_get_model_species(pkim, 0, &species);
  if (KIM_STATUS_OK > status)
    error(1, "KIM_API_get_model_species failed");

  // create buffer for quering the model parameters

  char desc[4096];
  memset(desc, 0, sizeof(desc));

  snprintf(desc, sizeof(desc), KIM_INPUT_FILE_TEMPLATE, species);

  // free the temporary object
  KIM_API_free(&pkim, &status);
  if (KIM_STATUS_OK > status)
    error(1, "KIM_API_free failed");

  // initialize KIM API object
  status = KIM_API_string_init(&pkim, desc, g_kim.model_name);
  if (KIM_STATUS_OK > status)
    error(1, "KIM_API_string_init failed");

  // Allocate memory for each data argument of initialized KIM object
  KIM_API_allocate(pkim, 1, 1, &status);
  if (KIM_STATUS_OK > status)
    error(1, "KIM_API_allocate failed");

  // call Model's init routine
  status = KIM_API_model_init(pkim);
  if (KIM_STATUS_OK > status)
    error(1, "KIM_API_model_init failed");

  int descriptor_str_len = 0;

  status = KIM_API_get_model_kim_str_len(g_kim.model_name, &descriptor_str_len);
  if (status < KIM_STATUS_OK)
    error(1, "KIM_API_get_model_kim_str_len failed");

  char* descriptor_str = NULL;

  status = KIM_API_get_model_kim_str(g_kim.model_name, &descriptor_str);
  if (status < KIM_STATUS_OK)
    error(1, "KIM_API_get_model_kim_str failed");

  int numfreeparams = 0;
  int maxstringlength = 0;

  // get the maxStringLength of free parameters
  status = KIM_API_get_num_free_params(pkim, &numfreeparams, &maxstringlength);
  if (KIM_STATUS_OK > status)
    error(1, "KIM_API_get_num_free_params failed");

  char const* p = strstr(descriptor_str, "PARAM_FREE_cutoff");

  // don't count cutoff if available
  if (p) {
    numfreeparams--;
    p++;
  } else {
    p = descriptor_str;
  }

  g_kim.freeparams.npar = numfreeparams;
  g_kim.freeparams.name = (char const**)Malloc(numfreeparams * sizeof(char*));
  g_kim.freeparams.value = (void**)Malloc(numfreeparams * sizeof(void*));
  g_kim.freeparams.type = (KIM_PARAM_TYPE*)Malloc(numfreeparams * sizeof(KIM_PARAM_TYPE));
  g_kim.freeparams.size = (int*)Malloc(numfreeparams * sizeof(int));
  g_kim.freeparams.rank = (int*)Malloc(numfreeparams * sizeof(int));
  g_kim.freeparams.shape = (int**)Malloc(numfreeparams * sizeof(int*));

  // read all remaining parameters

  int index = 0;

  while (*p && index < numfreeparams) {
    p = strstr(p + 1,"PARAM_FREE");
    if (!p)
      break;

    char name[255];
    char type[255];

    sscanf(p, "%s %s", name, type);

    if (strncmp(type, "int", 3) == 0) {
      g_kim.freeparams.type[index] = KIM_PARAM_TYPE_INT;
    } else if (strncmp(type, "double", 6) == 0) {
      g_kim.freeparams.type[index] = KIM_PARAM_TYPE_DOUBLE;
    } else {
      g_kim.freeparams.type[index] = KIM_PARAM_TYPE_UNKNOWN;
    }

    g_kim.freeparams.name[index] = (char const*)Malloc(strlen(name) + 1);
    strncpy((char*)g_kim.freeparams.name[index], name, strlen(name) + 1);

    index++;
  }

  free(descriptor_str);

  // read size, rank and shape (if available)

  for (int i = 0; i < numfreeparams; ++i) {
    intptr_t size = KIM_API_get_size(pkim, g_kim.freeparams.name[i], &status);
    if (KIM_STATUS_OK > status)
      error(1, "KIM_API_get_size failed");
    g_kim.freeparams.size[i] = size;

    intptr_t rank = KIM_API_get_rank(pkim, g_kim.freeparams.name[i], &status);
    if (KIM_STATUS_OK > status)
      error(1, "KIM_API_get_rank failed");
    g_kim.freeparams.rank[i] = rank;

    if (rank == 0) {
      g_kim.freeparams.shape[i] = NULL;
    } else {
      g_kim.freeparams.shape[i] = (int*)Malloc(rank * sizeof(int));
      rank = KIM_API_get_shape(pkim, g_kim.freeparams.name[i], g_kim.freeparams.shape[i], &status);
      if (KIM_STATUS_OK > status)
        error(1, "KIM_API_get_shape failed");
    }
  }

  // read and copy values

  for (int i = 0; i < numfreeparams; ++i) {
    void* data = KIM_API_get_data(pkim, g_kim.freeparams.name[i], &status);
    if (KIM_STATUS_OK > status)
      error(1, "KIM_API_get_data failed");
    int data_size = g_kim.freeparams.size[i] * sizeof(intptr_t);
    g_kim.freeparams.value[i] = Malloc(data_size);
    memcpy(g_kim.freeparams.value[i], data, data_size);
  }

  // read neighbor list and boundry conditions string

  char const* NBCstr = "\0";

  status = KIM_API_get_NBC_method(pkim, &NBCstr);
  if (KIM_STATUS_OK > status)
    error(1, "KIM_API_get_NBC_method failed!\n");

  if (strncmp("NEIGH_RVEC", NBCstr, 10) == 0) {
    g_kim.NBC = KIM_NEIGHBOR_TYPE_RVEC;
  } else if (strncmp("NEIGH_PURE", NBCstr, 10) == 0) {
    g_kim.NBC = KIM_NEIGHBOR_TYPE_PURE;
  } else if (strncmp("MI_OPBC", NBCstr, 7) == 0) {
    g_kim.NBC = KIM_NEIGHBOR_TYPE_OPBC;
  } else if (strncmp("CLUSTER", NBCstr, 7) == 0) {
    g_kim.NBC = KIM_NEIGHBOR_TYPE_CLUSTER;
  } else {
    g_kim.NBC = KIM_NEIGHBOR_TYPE_UNKNOWN;
  }

  // read neighbor list type

  g_kim.is_half_neighbors = KIM_API_is_half_neighbors(pkim, &status);
  if (KIM_STATUS_OK > status)
    error(1, "KIM_API_is_half_neighbors failed!\n");

  // determine if the KIM Model can compute the total energy

  g_kim.model_has_energy = 1;
  KIM_API_get_index(pkim, "energy", &status);
  if (status != KIM_STATUS_OK)
    g_kim.model_has_energy = 0;

  // determine if the KIM Model can compute the forces

  g_kim.model_has_forces = 1;
  KIM_API_get_index(pkim, "forces", &status);
  if (status != KIM_STATUS_OK)
    g_kim.model_has_forces = 0;

  // determine if the KIM Model can compute the virial

  g_kim.model_has_virial = 1;
  KIM_API_get_index(pkim, (char*) "virial", &status);
  if (status != KIM_STATUS_OK) {
    KIM_API_get_index(pkim, (char*) "process_dEdr", &status);
    if (status != KIM_STATUS_OK)
      g_kim.model_has_virial = 0;
  }

  // check for cutoff

  g_kim.freeparams.cutoff_is_free_param = 0;
  g_kim.freeparams.cutoff_value = 0.0;
  void* data = KIM_API_get_data(pkim, "PARAM_FREE_cutoff", &status);
  if (KIM_STATUS_OK == status) {
    g_kim.freeparams.cutoff_value = *(double*)data;
    g_kim.freeparams.cutoff_is_free_param = 1;
  } else {
    data = KIM_API_get_data(pkim, "cutoff", &status);
    if (KIM_STATUS_OK == status) {
      g_kim.freeparams.cutoff_value = *(double*)data;
    } else {
      error(1, "Cutoff is not available from the model!\n");
    }
  }

  KIM_API_model_destroy(pkim);
  KIM_API_free(&pkim, &status);
}

void print_model_parameters_and_exit(char const* filename)
{
  printf("\n");
  warning("check_kim_opt_param has been found in the potential file.\n");
  warning("All possible parameters will be listed an potfit will be terminated!\n\n");

  printf("The following potential parameters are available for fitting. Include the\n"
         "name(s) (and the initial value(s) and the corresponding lower and upper\n"
         "boundaries) that you want to optimize in your potential file %s.\n\n", filename);

  printf("Name\t\t\tType\t\tExtent\n");
  printf("\n");

  for(int i = 0; i < g_kim.freeparams.npar; ++i) {
    printf("%s\t\t", g_kim.freeparams.name[i]);
    switch (g_kim.freeparams.type[i]) {
      case KIM_PARAM_TYPE_UNKNOWN:
        printf("Unknown\t\t");
        break;
      case KIM_PARAM_TYPE_INT:
        printf("Integer\t\t");
        break;
      case KIM_PARAM_TYPE_DOUBLE:
        printf("Double\t\t");
        break;
    }

    printf("[ ");
    for(int j = 0; j < g_kim.freeparams.rank[i]; ++j) {
      printf("%d ", g_kim.freeparams.shape[i][j]);
    }
    printf("]\n");
  }

  printf("\nNote that an empty parameter extent (i.e. '[ ]') indicates that the\n"
         "parameter is a scalar.\n");
  printf("KIM array parameters are row based. While listing the initial\n"
         "values for such parameters, you should ensure that the sequence is\n"
         "correct. For example, if the extent of a parameter `PARAM_FREE_A' is\n"
         "[ 2 2 ], then you should list the initial values as: A[0 0], A[0 1],\n"
         "A[1 0], A[1 1].\n");

  exit(1);
}

/******************************************************************************
*
* Read the keywords (`type', `cutoff', etc) and parameter names (beginning
* with `PARAM_FREE') in the potential input file.
*
* In this function, the memory of three global variables: `num_opt_param',
* `name_opt_param', `size_opt_param' will be allocated and the first two will be
* initialized here (based on the infomation from the input), `size_opt_param' will
* be initialized in nest_optimizable_param.
*
* pt: potential table
* filename: potential input file name
* infile: potential input file
* FreeParam: data struct that contains the free parameters info of the KIM Model
*
******************************************************************************/

void read_keywords_and_parameters(char const* filename, FILE* infile)
{
  char buffer[255], name[255];
  fpos_t startpos;

  // save starting position
  fgetpos(infile, &startpos);

  // scan for "model" keyword

  do {
    if (NULL == fgets(buffer, 255, infile))
      error(1, "Error reading 'model' keyword from potential file %s\n", filename);
    if (strlen(buffer) > 1 && 1 != sscanf(buffer, "%s", name))
      error(1, "Error reading 'model' keyword from potential file %s\n", filename);
  } while (strncmp(name, "model", 5) != 0 && !feof(infile));

  if (strncmp(name, "model", 5) != 0)
    error(1, "Keyword 'model' is missing in file: %s\n", filename);

  if (1 != sscanf(buffer, "%*s %s", name))
    error(1, "Error reading 'model' keyword from potential file %s\n", filename);

  // copy name
  g_kim.model_name = (char const*)Malloc(strlen(name) + 1);
  strncpy((char*)g_kim.model_name, name, strlen(name));
  *((char*)&g_kim.model_name[strlen(name)]) = '\0';

  printf("\nKIM Model: %s\n", g_kim.model_name);

  read_kim_model_params();

  // find "kim_opt_param" keyword

  fsetpos(infile, &startpos);

  do {
    if (NULL == fgets(buffer, 255, infile))
      error(1, "Error reading 'kim_opt_param' keyword from potential file %s\n", filename);
    if (strlen(buffer) > 1 && 1 != sscanf(buffer, "%s", name))
      error(1, "Error reading 'kim_opt_param' keyword from potential file %s\n", filename);
  } while (strncmp(name, "kim_opt_param", 13) != 0 && !feof(infile));

  // read "kim_opt_param"

  if (strncmp(buffer,"kim_opt_param", 13) != 0) {
    error(1, "Keyword 'kim_opt_param' is missing in potential file: %s\n", filename);
  } else {
    if (1 != sscanf(buffer, "%*s %s", name))
      error(1, "Argument for 'kim_opt_param' is missing/invalid in potential file: %s\n", filename);
    if (strncmp(name, "list", 4) == 0) {
      print_model_parameters_and_exit(filename);
    } else {
      if(1 != sscanf(name, "%d", &g_kim.num_opt_param)) {
        error(0, "Argument for 'kim_opt_param' (%s) is invalid in file: %s\n", name, filename);
        error(1, "Please provide an integer or the keyword 'list'\n");
      }
    }
  }

  if (g_kim.num_opt_param > g_kim.freeparams.npar)
    error(1, "kim_opt_param in potential file is set to %d, but the "
      "model only has %d parameters.", g_kim.num_opt_param, g_kim.freeparams.npar);

  // allocate storage

  g_kim.idx_opt_param = (int*)Malloc(g_kim.num_opt_param * sizeof(int));
  g_kim.size_opt_param = (int*)Malloc(g_kim.num_opt_param * sizeof(int));

  // count how many optimization parameters we have
  // values are not read here

  fsetpos(infile, &startpos);

  for (int i = 0; i < g_kim.num_opt_param; ++i) {
    do {
      if (NULL == fgets(buffer, 255, infile))
        error(1, "Error reading 'PARAM_FREE' keyword #%d from potential file %s\n", i, filename);
      if (strlen(buffer) > 1 && 1 != sscanf(buffer, "%s", name))
        error(1, "Error reading 'PARAM_FREE' keyword #%d from potential file %s\n", i, filename);
    } while (strncmp(name, "PARAM_FREE", 10) != 0 && !feof(infile));
    if (feof(infile) ) {
        error(0, "Not enough parameter(s) 'PARAM_FREE_*' in file: %s.\n", filename);
        error(1, "You listed %d parameter(s), but required are %d.\n", i + 1, g_kim.num_opt_param);
    }
    int found = 0;
    for (int j = 0; j < g_kim.freeparams.npar; ++j) {
      if (strcmp(g_kim.freeparams.name[j], name) == 0) {
        if (g_kim.freeparams.type[j] != KIM_PARAM_TYPE_DOUBLE)
          error(1, "Only double parameters can currently be optimized!\n");
        g_kim.idx_opt_param[i] = j;
        g_kim.size_opt_param[i] = 1;
        for (int k = 0; k < g_kim.freeparams.rank[j]; ++k)
          g_kim.size_opt_param[i] *= g_kim.freeparams.shape[j][k];
        g_kim.total_num_opt_param += g_kim.size_opt_param[i];
        found = 1;
        break;
      }
    }
    if (!found)
      error(1, "The KIM model does not have a parameter with the name %s\n", name);
  }
}

double get_freeparam_value(int param_idx, int entry_idx)
{
  if (g_kim.freeparams.type[param_idx] != KIM_PARAM_TYPE_DOUBLE)
    error(1, "Only double values may currently be optimized\n");

  if (g_kim.freeparams.rank[param_idx] == 0) {
    if (entry_idx != 0) {
      error(1, "Parameter %d does only have a single entry!\n", param_idx);
    } else {
      return ((double*)g_kim.freeparams.value[param_idx])[0];
    }
  }

  if (g_kim.freeparams.rank[param_idx] == 1) {
    if (entry_idx > g_kim.freeparams.size[param_idx]) {
      error(1, "Parameter %d does only have %d entries!\n", param_idx, g_kim.freeparams.size[param_idx]);
    } else {
      return ((double*)g_kim.freeparams.value[param_idx])[entry_idx];
    }
  }

  if (g_kim.freeparams.rank[param_idx] == 2) {
    int i = entry_idx / g_kim.freeparams.shape[param_idx][1];
    int j = entry_idx % g_kim.freeparams.shape[param_idx][0];

    if (entry_idx > g_kim.freeparams.size[param_idx]) {
      error(1, "Parameter %d does only have %d entries!\n", param_idx, g_kim.freeparams.size[param_idx]);
    } else {
      return ((double**)g_kim.freeparams.value[param_idx])[i][j];
    }
  }

  error(1, "Only up to rank 2 is implemented here!\n");

  return -1.0;
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

  char  buffer[255], name[255];
  fpos_t startpos;

  // save starting position
  fgetpos(pfile, &startpos);

  // read the keywords and the names of parameters that will be optimized
  read_keywords_and_parameters(filename, pfile);

  apt->n_par[0] = g_kim.total_num_opt_param;
  apt->total_par += apt->n_par[0];

  // set begin
  apt->begin[0] = 0.0;

  // allocate memory
  apt->values[0] = (double *)Malloc(apt->n_par[0] * sizeof(double));
  // last slot stores the total number of invariable parameters
  apt->invar_par[0] = (int *)Malloc((apt->n_par[0] + 1) * sizeof(int));
  apt->pmin[0] = (double *)Malloc(apt->n_par[0] * sizeof(double));
  apt->pmax[0] = (double *)Malloc(apt->n_par[0] * sizeof(double));

  // read parameter values

  fsetpos(pfile, &startpos);

  int param_idx = 0;

  for (int j = 0; j < g_kim.num_opt_param; j++) {
    memset(buffer, 0, sizeof(buffer));
    memset(name, 0, sizeof(name));

    do {
      if (NULL == fgets(buffer, 255, pfile))
        error(1, "Error reading 'PARAM_FREE' keyword #%d from potential file %s\n", j, filename);
      if (strlen(buffer) > 1 && 1 != sscanf(buffer, "%s", name))
        error(1, "Error reading 'PARAM_FREE' keyword #%d from potential file %s\n", j, filename);
    } while (strncmp(name, "PARAM_FREE", 10) != 0 && !feof(pfile));

    if (strncmp(name, "PARAM_FREE_cutoff", 17) == 0) {
      j--;
      continue;
    }

    // make sure we have the same order as before
    const int idx = g_kim.idx_opt_param[j];
    if (strcmp(g_kim.freeparams.name[idx], name) != 0)
      error(1, "Parameter order mismatch, this should never happen!\n");

    double min = 0.0;
    double max = 0.0;

    // read all lines for this parameter
    for (int k = 0; k < g_kim.size_opt_param[j]; ++k) {
      if (NULL == fgets(buffer, 255, pfile))
        error(1, "Error reading '%s' entry #%d from potential file %s\n", name, k, filename);
      int ret_val = sscanf(buffer, "%s %lf %lf", name, &min, &max);
      if (ret_val == 3) {
        if (min > max) {
          double temp = 0.0;
          SWAP(min, max, temp);
        }
        apt->pmin[0][param_idx] = min;
        apt->pmax[0][param_idx] = max;
        // fixed param or not
        if (min == max) {
          apt->invar_par[0][param_idx] = 1;
          apt->invar_par[0][apt->n_par[0]] += 1;
        }
      }
      for (int c = 0; name[c]; ++c)
        name[c] = tolower(name[c]);
      if (strncmp(name, "kim", 3) == 0) {
        apt->values[0][param_idx] = get_freeparam_value(idx, k);
      } else if (1 > sscanf(name, "%lf", &apt->values[0][param_idx])) {
        error(1, "First data for parameter '%s' corrupted, it should be a type float or 'KIM'.\n", g_kim.freeparams.name[idx]);
      }
      if(min > apt->values[0][param_idx] || max < apt->values[0][param_idx])
        error(1, "Value %d of '%s' is not within its limits in file '%s'.\n", k + 1, g_kim.freeparams.name[idx], filename);
      param_idx++;
    }
  }

  printf("Successfully read %d potential parameters that will be optimized.\n", g_kim.total_num_opt_param);

  fsetpos(pfile, &startpos);

  do {
    if (NULL == fgets(buffer, 255, pfile))
      error(1, "Error while searching for 'cutoff' keyword in potential file %s\n", filename);
    if (strlen(buffer) > 1 && 1 != sscanf(buffer, "%s", name))
      error(1, "Error while searching for 'cutoff' keyword in potential file %s\n", filename);
  } while (strncmp(name, "cutoff", 6) != 0 && !feof(pfile));

  if (strncmp(name, "cutoff", 6) == 0) {
    char temp_str[255];

    if(2 != sscanf(buffer, "%s %s", name, temp_str))
      error(1,"Error reading in cutoff in file '%s'.\n", filename);

    for (int c = 0; temp_str[c]; ++c)
      temp_str[c] = tolower(temp_str[c]);

    if (strncmp(temp_str, "kim", 3) == 0) {
      printf("Cutoff from KIM model will be used: %f\n", g_kim.freeparams.cutoff_value);
      apt->end[0] = g_kim.freeparams.cutoff_value;
    } else {
      if (g_kim.freeparams.cutoff_is_free_param == 0) {
        warning("Cutoff is not a free parameter. Model value will be used!\n");
      } else {
        if(1 != sscanf(temp_str,"%lf", &apt->end[0]))
          error(1,"Error reading in cutoff in file '%s'.\n", filename);
        g_kim.freeparams.cutoff_value = apt->end[0];
      }
    }
  }

  // initialize optimization table
  pt->begin[0] = 0.0;
  pt->end[0] = apt->end[0];
  g_pot.calc_pot.begin = pt->begin;
  g_pot.calc_pot.end = pt->end;
  pt->step[0] = 0.0;
  pt->invstep[0] = 0.0;
  pt->first[0] = 2;
  pt->last[0] = pt->first[0] + apt->n_par[0] - 1;
  pt->len = pt->first[0] + apt->n_par[0];
  pt->table = (double *)Malloc(pt->len * sizeof(double));
  g_pot.calc_list = (double *)Malloc(pt->len * sizeof(double));
  pt->idx = (int *)Malloc(pt->len * sizeof(int));
  apt->idxpot = (int *)Malloc(apt->total_par * sizeof(int));
  apt->idxparam = (int *)Malloc(apt->total_par * sizeof(int));

  int k = 0;
  int l = 0;
  double * val = pt->table;
  double * list = g_pot.calc_list;

  // loop over potentials
  for (int i = 0; i < apt->number; ++i) {
    // loop over parameters
    for (int j = 0; j < apt->n_par[i]; ++j) {
      *val = apt->values[i][j];
      *list = apt->values[i][j];
      val++;
      list++;
      if (!g_pot.invar_pot[i] && !apt->invar_par[i][j]) {
        // index, the optimiable parameters in pt->table. for example,
        // idx[0] = 2 indicates that the the first variable parameters
        // lies in slot 2 of pt->table
        pt->idx[k] = l;
        // index, the optimizable parameters come from which potential?
        // e.g. idxpot[0] = 0, indicates that the first variable parameter
        // is from the first potential
        apt->idxpot[k] = i;
        // index, the optimizable parameters is the ?th parammeter of the
        // potential. Sould be used together with idxpot. e.g. idxparam[0] = 1
        // indicates that the first variable parameter is the 2nd parameter of
        // potential idxpot[0]
        apt->idxparam[k] = j;
      }
      k++;
      l++;
    }

    if (!g_pot.invar_pot[i])
      pt->idxlen += apt->n_par[i] - apt->invar_par[i][apt->n_par[i]];
    apt->total_par -= apt->invar_par[i][apt->n_par[i]];
  }

#if defined(NOPUNISH)
  if (opt)
    warning("Gauge degrees of freedom are NOT fixed!\n");
#endif // NOPUNISH
}

#endif // !APOT && !KIM
