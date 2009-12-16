/******************************************************************************
*
*  param.c: Read in parameter files (tag based)
*
******************************************************************************/
/*
*   Copyright 2002-2009 Peter Brommer, Franz G"ahler, Daniel Schopf
*             Institute for Theoretical and Applied Physics
*             University of Stuttgart, D-70550 Stuttgart, Germany
*             http://www.itap.physik.uni-stuttgart.de/
*
*****************************************************************/
/*
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
*   along with potfit; if not, write to the Free Software
*   Foundation, Inc., 51 Franklin St, Fifth Floor,
*   Boston, MA  02110-1301  USA
*/
/**************************************************************************
*  $Revision: 1.31 $
*  $Date: 2009/12/16 12:10:56 $
***************************************************************************/

#ifndef POTSCALE
#include "potfit.h"
#else
#include "potscale.h"
#endif


/* the following are needed for gettimeofday */
#include <sys/time.h>
#include <unistd.h>

#if defined(__GNUC__) && defined(__STRICT_ANSI__)
extern char *strdup(char *);
#endif

#ifdef __WATCOMC__
#define strcasecmp strcmpi
#endif

typedef enum ParamType {
  PARAM_STR, PARAM_STRPTR,
  PARAM_INT, PARAM_DOUBLE
} PARAMTYPE;

int   curline;			/* number of current line */

/******************************************************************************
*
*  get parameters from one line
*
******************************************************************************/

/* Parameter:

   param_name ... Parametername (fuer Fehlermeldungen)
   param ........ Adresse der Variable fuer den Parameter
   ptype ........ Parametertyp
                  folgende Werte sind zulaessig:
                  PARAM_STR : String, deklariert als char[]
                  PARAM_STRPTR : String, deklariert als Zeiger auf char*
                  PARAM_INT : Integer-Wert(e)
                  PARAM_DOUBLE : Double-Wert(e)

   pnum_min ..... Minimale Anzahl der einzulesenden Werte
                  (Falls weniger Werte gelesen werden koennen als verlangt,
                  wird ein Fehler gemeldet).
   pnum_max ..... Maximale Anzahl der einzulesenden Werte
                  (Ueberzaehlige Werte werden ignoriert. Falls weniger als
                  pnum_max Werte vorhanden sind, wird das Lesen
                  abgebrochen, es wird aber kein Fehler gemeldet,
                  wenn mindestens pnum_min Werte abgesaettigt wurden.

  Resultat: Die Anzahl der gelesenen Werte wird zurueckgegeben.

*/

int getparam(char *param_name, void *param, PARAMTYPE ptype,
	     int pnum_min, int pnum_max)
{
  static char errmsg[256];
  char *str;
  int   i;
  int   numread;

  numread = 0;
  if (ptype == PARAM_STR) {
    str = strtok(NULL, " \t\r\n");
    if (str == NULL) {
      sprintf(errmsg,
	      "Parameter for %s missing in line %u\nstring expected!\n",
	      param_name, curline);
      error(errmsg);
    } else
      strncpy((char *)param, str, pnum_max);
    numread++;
  } else if (ptype == PARAM_STRPTR) {
    str = strtok(NULL, " \t\r\n");
    if (str == NULL) {
      sprintf(errmsg,
	      "Parameter for %s missing in line %u\nstring expected!\n",
	      param_name, curline);
      error(errmsg);
    } else
      *((char **)param) = strdup(str);
    numread++;
  } else if (ptype == PARAM_INT) {
    for (i = 0; i < pnum_min; i++) {
      str = strtok(NULL, " \t\r\n");
      if (str == NULL) {
	sprintf(errmsg, "Parameter for %s missing in line %u!\n",
		param_name, curline);
	sprintf(errmsg + strlen(errmsg),
		"Integer vector of length %u expected!\n",
		(unsigned)pnum_min);
	error(errmsg);
      } else
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
  } else if (ptype == PARAM_DOUBLE) {
    for (i = 0; i < pnum_min; i++) {
      str = strtok(NULL, " \t\r\n");
      if (str == NULL) {
	sprintf(errmsg, "Parameter for %s missing in line %u!\n",
		param_name, curline);
	sprintf(errmsg + strlen(errmsg),
		"Double vector of length %u expected!\n", (unsigned)pnum_min);
	error(errmsg);
      } else
	((real *)param)[i] = atof(str);
      numread++;
    }
    for (i = pnum_min; i < pnum_max; i++) {
      if ((str = strtok(NULL, " \t\r\n")) != NULL) {
	((real *)param)[i] = atof(str);
	numread++;
      } else
	break;
    }
  }
  return numread;
}

/******************************************************************************
*
*  read in parameter file
*  lines beginning with comment characters '#' or blank lines are skipped
*
******************************************************************************/

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
/*     file for end force */
/*    else if (strcasecmp(token, "endforce") == 0) {*/
/*      getparam("endforce", endforce, PARAM_STR, 1, 255);*/
/*    }*/
/*     file for end energy */
/*    else if (strcasecmp(token, "endenergy") == 0) {*/
/*      getparam("endenergy", endenergy, PARAM_STR, 1, 255);*/
/*    }*/
/*     file for end stress */
/*    else if (strcasecmp(token, "endstress") == 0) {*/
/*      getparam("endstress", endstress, PARAM_STR, 1, 255);*/
/*    }*/
    /* prefix for all output files */
    else if (strcasecmp(token, "output_prefix") == 0) {
      getparam("output_prefix", output_prefix, PARAM_STR, 1, 255);
      if (strcmp(output_prefix, "") != 0)
	write_output_files = 1;
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
    /* file for maximal change */
    else if (strcasecmp(token, "maxchfile") == 0) {
      getparam("maxchfile", maxchfile, PARAM_STR, 1, 255);
      usemaxch = 1;
    }
    /* file for pair distribution */
    else if (strcasecmp(token, "distfile") == 0) {
      getparam("distfile", distfile, PARAM_STR, 1, 255);
    }
    /* number of steps in IMD potential */
    else if (strcasecmp(token, "imdpotsteps") == 0) {
      getparam("imdpotsteps", &imdpotsteps, PARAM_INT, 1, 1);
    }
#ifdef APOT
    /* minimum for plotfile */
    else if (strcasecmp(token, "plotmin") == 0) {
      getparam("plotmin", &plotmin, PARAM_DOUBLE, 1, 1);
    }
#ifndef EAM
    /* exclude chemical potential from energy calculations */
    else if (strcasecmp(token, "disable_cp") == 0) {
      getparam("disable_cp", &disable_cp, PARAM_INT, 1, 1);
    }
#endif
#endif
    /* Energy Weight */
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
    /* starting temperature for annealing */
    else if (strcasecmp(token, "anneal_temp") == 0) {
      getparam("anneal_temp", &anneal_temp, PARAM_DOUBLE, 1, 1);
    }
    /* starting width for normal distribution for evo */
    else if (strcasecmp(token, "evo_width") == 0) {
      getparam("evo_width", &evo_width, PARAM_DOUBLE, 1, 1);
    }
    /* Energy Weight */
    else if (strcasecmp(token, "eng_weight") == 0) {
      getparam("eng_weight", &eweight, PARAM_DOUBLE, 1, 1);
    }
#ifdef STRESS
    /* Energy Weight */
    else if (strcasecmp(token, "stress_weight") == 0) {
      getparam("stress_weight", &sweight, PARAM_DOUBLE, 1, 1);
    }
#endif
    /* unknown tag */
    else {
      fprintf(stderr, "Unknown tag <%s> ignored!\n", token);
      fflush(stderr);
    }
  } while (!feof(pf));
  fclose(pf);

}

/******************************************************************************
*
*  read command line
*
******************************************************************************/

void read_parameters(int argc, char **argv)
{
  char  msg[255];
  FILE *pf;

  /* check command line */
  if (argc < 2) {
    sprintf(msg, "Usage: %s <paramfile>\n", argv[0]);
    error(msg);
  }

  /* open parameter file, and read it */
  pf = fopen(argv[1], "r");
  if (NULL == pf) {
    fprintf(stderr, "ERROR: Could not open parameter file %s!\n", argv[1]);
    fflush(stderr);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 2);
#endif
    exit(2);
  }
  read_paramfile(pf);
  printf("Read parameters from file %s\n", argv[1]);

}
