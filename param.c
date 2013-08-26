/****************************************************************
 *
 * param.c: Read in parameter files (tag based)
 *
 ****************************************************************
 *
 * Copyright 2002-2013
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

int   curline;			/* number of current line */

/****************************************************************
 *
 *  get parameters from one line
 *
 ****************************************************************/

/* parameters:

   param_name ... parametername (for error messages)
   param ........ address of the variable for that parameter
   ptype ........ type of parameter
                  the following are possible:
                  PARAM_STR : string, declared as char[]
                  PARAM_INT : integer value(s)
                  PARAM_DOUBLE : double value(s)

   pnum_min ..... minimum number of parameters to read
   		  (If there are less parameters than requested an
		  error will be reported.)
   pnum_max ..... maximum number of parameters to read
   		  (Additional values will be ignored. If less than pnum_max
		  values are registered, the reading process is stopped; there
		  will be no error if at least pnum_min values were read.)

   return value . the number of parameters read is returned

*/

int getparam(char *param_name, void *param, param_t ptype, int pnum_min, int pnum_max)
{
  char *str;
  int   i;
  int   numread;

  numread = 0;
  switch (ptype) {
      case PARAM_STR:
	str = strtok(NULL, " \t\r\n");
	if (str == NULL)
	  error(1, "Parameter for %s missing in line %d\nString expected!", param_name, curline);
	else
	  strncpy((char *)param, str, pnum_max);
	numread++;
	break;
      case PARAM_INT:
	for (i = 0; i < pnum_min; i++) {
	  str = strtok(NULL, " \t\r\n");
	  if (str == NULL)
	    error(1,
	      "Parameter for %s missing in line %d!\nInteger vector of length %u expected!",
	      param_name, curline, (unsigned)pnum_min);
	  else
	    ((int *)param)[i] = atoi(str);
	  numread++;
	}
	for (i = pnum_min; i < pnum_max; i++) {
	  if ((str = strtok(NULL, " \t\r\n")) != NULL) {
	    ((int *)param)[i] = atoi(str);
	    numread++;
	  } else
	    break;
	}
	break;
      case PARAM_DOUBLE:
	for (i = 0; i < pnum_min; i++) {
	  str = strtok(NULL, " \t\r\n");
	  if (str == NULL)
	    error(1,
	      "Parameter for %s missing in line %d!\nDouble vector of length %u expected!",
	      param_name, curline, (unsigned)pnum_min);
	  else
	    ((double *)param)[i] = atof(str);
	  numread++;
	}
	for (i = pnum_min; i < pnum_max; i++) {
	  if ((str = strtok(NULL, " \t\r\n")) != NULL) {
	    ((double *)param)[i] = atof(str);
	    numread++;
	  } else
	    break;
	}
	break;
  }

  return numread;
}

/****************************************************************
 *
 *  check all parameters for reasonable values and completeness
 *
 ****************************************************************/

void check_parameters_complete(char *paramfile)
{
  if (ntypes <= 0)
    error(1, "Missing parameter or invalid value in %s : ntypes is \"%d\"", paramfile, ntypes);

  if (strcmp(startpot, "\0") == 0)
    error(1, "Missing parameter or invalid value in %s : startpot is \"%s\"", paramfile, startpot);

  if (strcmp(endpot, "\0") == 0) {
    warning(1, "endpot is missing in %s, setting it to %s_end", paramfile, startpot);
    sprintf(endpot, "%s_end", startpot);
  }

  if (strcmp(config, "\0") == 0)
    error(1, "Missing parameter or invalid value in %s : config is \"%s\"", paramfile, config);

  if (strcmp(tempfile, "\0") == 0)
    error(1, "Missing parameter or invalid value in %s : tempfile is \"%s\"", paramfile, tempfile);

  if (eweight < 0)
    error(1, "Missing parameter or invalid value in %s : eng_weight is \"%f\"", paramfile, eweight);

#ifdef STRESS
  if (sweight < 0)
    error(1, "Missing parameter or invalid value in %s : stress_weight is \"%f\"", paramfile, sweight);
#endif /* STRESS */

  if (writeimd && imdpotsteps <= 0)
    error(1, "Missing parameter or invalid value in %s : imdpotsteps is \"%d\"", paramfile, imdpotsteps);

#ifdef APOT
  if (plotmin < 0)
    error(1, "Missing parameter or invalid value in %s : plotmin is \"%f\"", paramfile, plotmin);
#ifdef PAIR
  if (enable_cp != 0 && enable_cp != 1)
    error(1, "Missing parameter or invalid value in %s : enable_cp is \"%d\"", paramfile, enable_cp);
#endif /* PAIR */
#endif /* APOT */

#ifdef PDIST
  if (strcmp(distfile, "\0") == 0)
    error(1, "Missing parameter or invalid value in %s : distfile is \"%s\"", paramfile, distfile);
#endif /* PDIST */

  return;
}

/****************************************************************
 *
 *  read command line
 *
 ****************************************************************/

void read_parameters(int argc, char **argv)
{
  FILE *pf;

  /* check command line */
  if (argc < 2)
    error(1, "Usage: %s <paramfile>\n", argv[0]);

  /* open parameter file, and read it */
  pf = fopen(argv[1], "r");
  if (NULL == pf)
    error(1, "Could not open parameter file \"%s\".\n", argv[1]);
  read_paramfile(pf);
  fclose(pf);
  printf("Reading parameter file >> %s << ... ", argv[1]);
  check_parameters_complete(argv[1]);
  printf("done\n");

  return;
}

/****************************************************************
 *
 *  read in parameter file
 *  lines beginning with comment characters '#' or blank lines are skipped
 *
 ****************************************************************/

void read_paramfile(FILE *pf)
{
  char  buffer[1024];
  char *token, *res;

  curline = 0;

  do {
    res = fgets(buffer, 1024, pf);
    if (NULL == res)
      break;			/* probably EOF reached */
    curline++;
    token = strtok(buffer, " \t\r\n");
    if (NULL == token)
      continue;			/* skip blank lines */
    if (token[0] == '#')
      continue;			/* skip comments */

    /* number of atom types */
    if (strcasecmp(token, "ntypes") == 0) {
      getparam("ntypes", &ntypes, PARAM_INT, 1, 1);
    }
    /* file with start potential */
    else if (strcasecmp(token, "startpot") == 0) {
      getparam("startpot", startpot, PARAM_STR, 1, 255);
    }
    /* file for end potential */
    else if (strcasecmp(token, "endpot") == 0) {
      getparam("endpot", endpot, PARAM_STR, 1, 255);
    }
    /* prefix for all output files */
    else if (strcasecmp(token, "output_prefix") == 0) {
      getparam("output_prefix", output_prefix, PARAM_STR, 1, 255);
      if (strcmp(output_prefix, "") != 0)
	write_output_files = 1;
    } else if (strcasecmp(token, "output_lammps") == 0) {
      getparam("output_lammps", output_lammps, PARAM_STR, 1, 255);
      if (strcmp(output_lammps, "") != 0)
	write_lammps_files = 1;
    }
    /* file for IMD potential */
    else if (strcasecmp(token, "imdpot") == 0) {
      getparam("imdpot", imdpot, PARAM_STR, 1, 255);
      writeimd = 1;
    }
    /* file for plotting */
    else if (strcasecmp(token, "plotfile") == 0) {
      getparam("plotfile", plotfile, PARAM_STR, 1, 255);
      plot = 1;
    }
#ifndef APOT
    /* file for maximal change */
    else if (strcasecmp(token, "maxchfile") == 0) {
      getparam("maxchfile", maxchfile, PARAM_STR, 1, 255);
      usemaxch = 1;
    }
#endif /* !APOT */
#ifdef PDIST
    /* file for pair distribution */
    else if (strcasecmp(token, "distfile") == 0) {
      getparam("distfile", distfile, PARAM_STR, 1, 255);
    }
#endif /* PDIST */
    /* number of steps in IMD potential */
    else if (strcasecmp(token, "imdpotsteps") == 0) {
      getparam("imdpotsteps", &imdpotsteps, PARAM_INT, 1, 1);
    }
#ifdef APOT
    /* minimum for plotfile */
    else if (strcasecmp(token, "plotmin") == 0) {
      getparam("plotmin", &plotmin, PARAM_DOUBLE, 1, 1);
    }
#ifdef PAIR
    /* exclude chemical potential from energy calculations */
    else if (strcasecmp(token, "enable_cp") == 0) {
      getparam("enable_cp", &enable_cp, PARAM_INT, 1, 1);
    }
#endif /* PAIR */
#endif /* APOT */
    /* how far should the imd pot be extended */
    else if (strcasecmp(token, "extend") == 0) {
      getparam("extend", &extend, PARAM_DOUBLE, 1, 1);
    }
    /* file with atom configuration */
    else if (strcasecmp(token, "config") == 0) {
      getparam("config", config, PARAM_STR, 1, 255);
    }
    /* Optimization flag */
    else if (strcasecmp(token, "opt") == 0) {
      getparam("opt", &opt, PARAM_INT, 1, 1);
    }
    /* break flagfile */
    else if (strcasecmp(token, "flagfile") == 0) {
      getparam("flagfile", flagfile, PARAM_STR, 1, 255);
    }
    /* write radial pair distribution ? */
    else if (strcasecmp(token, "write_pair") == 0) {
      getparam("write_pair", &write_pair, PARAM_INT, 1, 1);
    }
    /* plotpoint file */
    else if (strcasecmp(token, "plotpointfile") == 0) {
      getparam("plotpointfile", plotpointfile, PARAM_STR, 1, 255);
    }
    /* temporary potential file */
    else if (strcasecmp(token, "tempfile") == 0) {
      getparam("tempfile", tempfile, PARAM_STR, 1, 255);
    }
    /* seed for RNG */
    else if (strcasecmp(token, "seed") == 0) {
      getparam("seed", &seed, PARAM_INT, 1, 1);
    }
#ifdef EVO
    /* stopping criterion for differential evolution */
    else if (strcasecmp(token, "evo_threshold") == 0) {
      getparam("evo_threshold", &evo_threshold, PARAM_DOUBLE, 1, 1);
    }
#else /* EVO */
    /* starting temperature for annealing */
    else if (strcasecmp(token, "anneal_temp") == 0) {
      getparam("anneal_temp", &anneal_temp, PARAM_STR, 1, 20);
    }
#endif /* EVO */
#ifdef APOT
    /* Scaling Constant for APOT Punishment */
    else if (strcasecmp(token, "apot_punish") == 0) {
      getparam("apot_punish", &apot_punish_value, PARAM_DOUBLE, 1, 1);
    }
#endif /* APOT */
    /* Energy Weight */
    else if (strcasecmp(token, "eng_weight") == 0) {
      getparam("eng_weight", &eweight, PARAM_DOUBLE, 1, 1);
    }
    /* error margin delta epsilon */
    else if (strcasecmp(token, "d_eps") == 0) {
      getparam("d_eps", &d_eps, PARAM_DOUBLE, 1, 1);
    }
#ifdef STRESS
    /* Energy Weight */
    else if (strcasecmp(token, "stress_weight") == 0) {
      getparam("stress_weight", &sweight, PARAM_DOUBLE, 1, 1);
    }
#endif /* STRESS */
    /* write final potential in lammps format */
    else if (strcasecmp(token, "write_lammps") == 0) {
      getparam("write_lammps", &write_lammps, PARAM_INT, 1, 1);
    }
#ifdef COULOMB
    /* cutoff-radius for long-range interactions */
    else if (strcasecmp(token, "dp_cut") == 0) {
      getparam("dp_cut", &dp_cut, PARAM_DOUBLE, 1, 1);
    }
#endif /* COULOMB */
#ifdef DIPOLE
    /* dipole iteration precision */
    else if (strcasecmp(token, "dp_tol") == 0) {
      getparam("dp_tol", &dp_tol, PARAM_DOUBLE, 1, 1);
    }
    /* mixing parameter for damping dipole iteration loop */
    else if (strcasecmp(token, "dp_mix") == 0) {
      getparam("dp_mix", &dp_mix, PARAM_DOUBLE, 1, 1);
    }
#endif /* DIPOLE */
    /* unknown tag */
    else {
      fprintf(stderr, "Unknown tag <%s> in parameter file ignored!\n", token);
    }
  } while (!feof(pf));

  return;
}
