/****************************************************************
 *
 * param.c: potfit routines for reading the parameter file
 *
 ****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
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

#include <float.h>
#include <limits.h>
#include <strings.h>

#include "potfit.h"

#include "memory.h"
#include "params.h"

void read_parameter_file(char const* param_file);
void get_param_double(char const* param_name, double* value, int line,
                      char const* param_file, double min, double max);
void get_param_int(char const* param_name, int* value, int line,
                   char const* param_file, int min, int max);
void get_param_string(char const* param_name, const char** value, int line,
                      char const* param_file);
void check_parameters_complete(char const* paramfile);

/****************************************************************
  read_parameters
****************************************************************/

void read_parameters(int argc, char** argv)
{
  if (argc != 2)
    error(1, "Usage: %s <paramfile>\n", argv[0]);

  printf("Reading parameter file >> %s << ... \n", argv[1]);

  read_parameter_file(argv[1]);

  check_parameters_complete(argv[1]);

  printf("Reading parameter file >> %s << ... done\n", argv[1]);
}

/****************************************************************
  read_parameter_file
****************************************************************/

// read parameter file line by line and process keywords
// lines beginning with comment characters '#' or blank lines are skipped

void read_parameter_file(char const* param_file)
{
  char buffer[1024];
  char* token = NULL;

  int line = 0;

  // open parameter file, and read it
  FILE* pfile = fopen(param_file, "r");

  if (pfile == NULL)
    error(1, "Could not open parameter file %s.\n", param_file);

  do {
    if (fgets(buffer, 1024, pfile) == NULL)
      break;  // probably EOF reached

    line++;

    token = strtok(buffer, " \t\r\n");

    // skip blank and comment lines
    if (token == NULL || token[0] == '#')
      continue;

    // number of atom types
    if (strcasecmp(token, "ntypes") == 0) {
      get_param_int("ntypes", &g_param.ntypes, line, param_file, 1, INT_MAX);
    }
    // file with start potential
    else if (strcasecmp(token, "startpot") == 0) {
      get_param_string("startpot", &g_files.startpot, line, param_file);
    }
    // file for end potential
    else if (strcasecmp(token, "endpot") == 0) {
      get_param_string("endpot", &g_files.endpot, line, param_file);
    }
    // prefix for general output files
    else if (strcasecmp(token, "output_prefix") == 0) {
      get_param_string("output_prefix", &g_files.output_prefix, line,
                       param_file);
      g_param.write_output_files = 1;
    }
    // prefix for LAMMPS output files
    else if (strcasecmp(token, "output_lammps") == 0) {
      get_param_string("output_lammps", &g_files.output_lammps, line,
                       param_file);
      g_param.write_lammps_files = 1;
    }
    // file for IMD potential
    else if (strcasecmp(token, "imdpot") == 0) {
      get_param_string("imdpot", &g_files.imdpot, line, param_file);
      g_param.writeimd = 1;
    }
    // file for plotting
    else if (strcasecmp(token, "plotfile") == 0) {
      get_param_string("plotfile", &g_files.plotfile, line, param_file);
      g_param.plot = 1;
    }
    // number of steps in IMD potential
    else if (strcasecmp(token, "imdpotsteps") == 0) {
      get_param_int("imdpotsteps", &g_param.imdpotsteps, line, param_file, 1,
                    INT_MAX);
    }
    // how far should the imd pot be extended
    else if (strcasecmp(token, "extend") == 0) {
      get_param_double("extend", &g_param.extend, line, param_file, DBL_MIN,
                       DBL_MAX);
    }
    // file with atom configuration
    else if (strcasecmp(token, "config") == 0) {
      get_param_string("config", &g_files.config, line, param_file);
    }
    // Optimization flag
    else if (strcasecmp(token, "opt") == 0) {
      get_param_int("opt", &g_param.opt, line, param_file, 0, 1);
    }
    // break flagfile
    else if (strcasecmp(token, "flagfile") == 0) {
      get_param_string("flagfile", &g_files.flagfile, line, param_file);
    }
    // write radial pair distribution ?
    else if (strcasecmp(token, "write_pair") == 0) {
      get_param_int("write_pair", &g_param.write_pair, line, param_file, 0, 1);
    }
    // plotpoint file
    else if (strcasecmp(token, "plotpointfile") == 0) {
      get_param_string("plotpointfile", &g_files.plotpointfile, line,
                       param_file);
    }
    // temporary potential file
    else if (strcasecmp(token, "tempfile") == 0) {
      get_param_string("tempfile", &g_files.tempfile, line, param_file);
    }
    // seed for RNG
    else if (strcasecmp(token, "seed") == 0) {
      get_param_int("seed", &g_param.rng_seed, line, param_file, INT_MIN,
                    INT_MAX);
    }
    // Energy Weight
    else if (strcasecmp(token, "eng_weight") == 0) {
      get_param_double("eng_weight", &g_param.eweight, line, param_file, 0,
                       DBL_MAX);
    }
    // error margin delta epsilon
    else if (strcasecmp(token, "d_eps") == 0) {
      get_param_double("d_eps", &g_calc.d_eps, line, param_file, 0, DBL_MAX);
    }

    // write final potential in lammps format
    else if (strcasecmp(token, "write_lammps") == 0) {
      get_param_int("write_lammps", &g_param.write_lammps, line, param_file, 0,
                    1);
    }
    // global scaling parameter
    else if (strcasecmp(token, "cell_scale") == 0) {
      get_param_double("cell_scale", &g_param.global_cell_scale, line,
                       param_file, 0, DBL_MAX);
    }
#if defined(APOT)
    // Scaling Constant for APOT Punishment
    else if (strcasecmp(token, "apot_punish") == 0) {
      get_param_double("apot_punish", &g_param.apot_punish_value, line,
                       param_file, 0, DBL_MAX);
    }
    // minimum for plotfile
    else if (strcasecmp(token, "plotmin") == 0) {
      get_param_double("plotmin", &g_param.plotmin, line, param_file, DBL_MIN,
                       DBL_MAX);
    }
#if defined(PAIR)
    // exclude chemical potential from energy calculations
    else if (strcasecmp(token, "enable_cp") == 0) {
      get_param_int("enable_cp", &g_param.enable_cp, line, param_file, INT_MIN,
                    INT_MAX);
    }
#endif  // PAIR
#else   // APOT
    // file for maximal change
    else if (strcasecmp(token, "maxchfile") == 0) {
      get_param_string("maxchfile", &g_files.maxchfile, line, param_file);
      g_param.usemaxch = 1;
    }
#endif  // APOT

#if defined(COULOMB)
    // cutoff-radius for long-range interactions
    else if (strcasecmp(token, "dp_cut") == 0) {
      get_param_double("dp_cut", &g_config.dp_cut, line, param_file, DBL_MIN,
                       DBL_MAX);
    }
#endif  // COULOMB

#if defined(DIPOLE)
    // dipole iteration precision
    else if (strcasecmp(token, "dp_tol") == 0) {
      get_param_double("dp_tol", &g_config.dp_tol, line, param_file, DBL_MIN,
                       DBL_MAX);
    }
    // mixing parameter for damping dipole iteration loop
    else if (strcasecmp(token, "dp_mix") == 0) {
      get_param_double("dp_mix", &g_config.dp_mix, line, param_file, DBL_MIN,
                       DBL_MAX);
    }
#endif  // DIPOLE

#if defined(EVO)
    // stopping criterion for differential evolution
    else if (strcasecmp(token, "evo_threshold") == 0) {
      get_param_double("evo_threshold", &g_param.evo_threshold, line,
                       param_file, 0, DBL_MAX);
    }
#else
    // starting temperature for annealing
    else if (strcasecmp(token, "anneal_temp") == 0) {
      get_param_string("anneal_temp", &g_param.anneal_temp, line, param_file);
    }
#endif  // EVO

#if defined(PDIST)
    // file for pair distribution
    else if (strcasecmp(token, "distfile") == 0) {
      get_param_string("distfile", &g_files.distfile, line, param_file);
    }
#endif  // PDIST

#if defined(STRESS)
    // Stress Weight
    else if (strcasecmp(token, "stress_weight") == 0) {
      get_param_double("stress_weight", &g_param.sweight, line, param_file, 0,
                       DBL_MAX);
    }
#endif  // STRESS

    // unknown tag
    else
      warning("Unknown tag <%s> in parameter file ignored!\n", token);

  } while (!feof(pfile));

  fclose(pfile);
}

/****************************************************************
  read double value
****************************************************************/

void get_param_int(char const* param_name, int* value, int line,
                   char const* param_file, int min, int max)
{
  char* str = strtok(NULL, " \t\r\n");

  if (str == NULL) {
    error(0, "Missing value in parameter file %s (line %d): %s is <undefined>\n",
          param_file, line, param_name);
  } else
    *value = atoi(str);

  if (*value < min || *value > max) {
    error(0,
          "Illegal value in parameter file %s (line %d): %s is out of bounds! "
          "(min: %i, max: %i, given: %i)\n",
          param_file, line, param_name, min, max, *value);
  }
}

/****************************************************************
  read integer value
****************************************************************/

void get_param_double(char const* param_name, double* value, int line,
                      char const* param_file, double min, double max)
{
  char* str = strtok(NULL, " \t\r\n");

  if (str == NULL) {
    error(0, "Missing value in parameter file %s (line %d): %s is <undefined>\n",
          param_file, line, param_name);
  } else
    *value = atof(str);

  if (*value < min || *value > max) {
    error(0,
          "Illegal value in parameter file %s (line %d): %s is out of bounds! "
          "(min: %d, max: %d, given: %d)\n",
          param_file, line, param_name, min, max, *value);
  }
}

/****************************************************************
  read string value
****************************************************************/

void get_param_string(char const* param_name, const char** value, int line,
                      char const* param_file)
{
  char* str = strtok(NULL, " \t\r\n");

  if (str == NULL) {
    error(0, "Missing value in parameter file %s (line %d): %s is <undefined>\n",
          param_file, line, param_name);
  } else {
    int len = strlen(str) + 1;
    *value = (const char*)Malloc(len * sizeof(char));
    strncpy((char*)*value, str, len);
  }
}

/****************************************************************
  check_parameters_complete
****************************************************************/

void check_parameters_complete(char const* paramfile)
{
  if (g_param.ntypes < 1)
    error(1,
          "Missing parameter or invalid value in %s : g_param.ntypes is \"%d\"\n",
          paramfile, g_param.ntypes);

  if (g_files.startpot == NULL)
    error(0,
          "Missing parameter or invalid value in %s : startpot is <undefined>\n",
          paramfile);

  if (g_files.endpot == NULL) {
    warning("endpot is missing in %s, setting it to %s_end\n", paramfile,
            g_files.startpot);
    g_files.endpot =
        (const char*)Malloc((strlen(g_files.startpot) + 5) * sizeof(char));
    sprintf((char*)g_files.endpot, "%s_end", g_files.startpot);
  }

  if (g_files.config == NULL)
    error(1, "Missing parameter or invalid value in %s : config is <undefined>\n",
          paramfile);

  if (g_files.tempfile == NULL)
    error(1,
          "Missing parameter or invalid value in %s : tempfile is <undefined>\n",
          paramfile);

  if (g_param.eweight < 0)
    error(1, "Missing parameter or invalid value in %s : eng_weight is \"%f\"\n",
          paramfile, g_param.eweight);

#if defined(STRESS)
  if (g_param.sweight < 0)
    error(1,
          "Missing parameter or invalid value in %s : stress_weight is \"%f\"\n",
          paramfile, g_param.sweight);
#endif  // STRESS

  if (g_param.writeimd && g_param.imdpotsteps < 1)
    error(1, "Missing parameter or invalid value in %s : imdpotsteps is \"%d\"\n",
          paramfile, g_param.imdpotsteps);

#if defined(APOT)
  if (g_param.plotmin < 0)
    error(1, "Missing parameter or invalid value in %s : plotmin is \"%f\"\n",
          paramfile, g_param.plotmin);
#if defined(PAIR)
  if (g_param.enable_cp != 0 && g_param.enable_cp != 1)
    error(1, "Missing parameter or invalid value in %s : enable_cp is \"%d\"\n",
          paramfile, g_param.enable_cp);
#endif  // PAIR
#endif  // APOT

#if defined(EVO)
    if (g_param.evo_threshold < 0)
      error(1, "Missing parameter or invalid value in %s : evo_threshold is"
                "\"%f\"\n", paramfile, g_param.evo_threshold);
#else
    if (g_param.anneal_temp == NULL) {
      warning("anneal_temp not provided in %s, setting it to 0!\n", paramfile);
      g_param.anneal_temp = (char*)Malloc(2);
      ((char*)g_param.anneal_temp)[0] = '0';
      ((char*)g_param.anneal_temp)[1] = '\0';
    }
#endif  // EVO

#if defined(PDIST)
  if (g_files.distfile == NULL)
    error(1,
          "Missing parameter or invalid value in %s : distfile is <undefined>\n",
          paramfile, g_files.distfile);
#endif  // PDIST

  if (g_param.global_cell_scale <= 0)
    error(1, "Missing parameter or invalid value in %s : cell_scale is \"%f\"\n",
          paramfile, g_param.global_cell_scale);
}
