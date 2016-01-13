/****************************************************************
 *
 * potential_input.c: Routines for reading a potential table
 *
 ****************************************************************
 *
 * Copyright 2002-2014
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
#include "potential.h"
#include "utils.h"

/*added */
#ifdef KIM
#include "kim/kim.h"
#endif
/*added ends*/

/****************************************************************
 *
 * read potential tables
 *
 ****************************************************************/

void read_pot_table(pot_table_t *pt, char *filename)
{
  FILE *infile;
  char  buffer[1024], *res, *str;
  int   have_format = 0, end_header = 0;
  int   size, i, j, k = 0, npots = 0;

#ifdef APOT
  apot_table_t *apt = &apot_table;
#endif /* APOT */

  /* set paircol to the number of pair potentials */
  paircol = (ntypes * (ntypes + 1)) / 2;

  /* open file */
  infile = fopen(filename, "r");
  if (NULL == infile)
    error(1, "Could not open file %s\n", filename);

  printf("Starting to read the potential file:\n");

  /* read the header */
  do {
    /* read one line */
    res = fgets(buffer, 1024, infile);
    if (NULL == res)
      error(1, "Unexpected end of file in %s", filename);
    /* check if it is a header line */
    if (buffer[0] != '#')
      error(1, "Header corrupt in file %s", filename);
    /* stop after last header line */
    if (buffer[1] == 'E') {
      end_header = 1;
    }

    /* see if it is the format line */
    if (buffer[1] == 'F') {
      /* format complete? */
      if (2 != sscanf((const char *)(buffer + 2), "%d %d", &format, &size))
	error(1, "Corrupt format header line in file %s", filename);

/* added */
#ifndef KIM 

#ifndef APOT
      if (format == 0)
	error(1, "potfit binary compiled without analytic potential support.\n");
#else
      if (format > 0)
	error(1, "potfit binary compiled without tabulated potential support.\n");
#endif /* !APOT */

      if (format != 0)
	printf(" - Potential file format %d detected\n", format);
      else
	printf(" - Potential file format %d (analytic potentials) detected\n", format);

      /* only pair potentials for
       * - pair interactions
       * - coulomb interactions
       * - dipole interactions
       */
      npots = paircol;

      /* more potential functions for other interactions */
#ifdef EAM
      npots = paircol + 2 * ntypes;
#ifdef TBEAM			/* TBEAM */
      npots += 2 * ntypes;
#endif /* TBEAM */
#endif /* EAM */

#ifdef ADP
      npots = 3 * paircol + 2 * ntypes;
#endif /* ADP */

#ifdef MEAM
      npots = 2 * paircol + 3 * ntypes;
#endif /* MEAM */

#ifdef STIWEB
      npots = 2 * paircol + 1;
#endif /* STIWEB */

#if defined TERSOFF && !defined TERSOFFMOD
      npots = ntypes * ntypes;
#endif /* TERSOFF && !TERSOFFMOD */

      if (size == npots) {
	printf(" - Using %d %s potentials to calculate forces\n", npots, interaction_name);
	fflush(stdout);
      } else {
	error(0, "Wrong number of data columns in %s potential file \"%s\".\n", interaction_name, filename);
	error(1, "For ntypes=%d there should be %d, but there are %d.", ntypes, npots, size);
      }
      /* recognized format? */
      if ((format != 0) && (format != 3) && (format != 4))
	error(1, "Unrecognized potential format specified for file %s", filename);


#else /* !KIM */
      if(format != 5)
        error(1,"Potential format in file '%s' should be '5' to use KIM potential.\n",filename);
      if(size != 1) {
        warning("The number of potentials should always be '1' when KIM Potential is used.\n"
            "You specified %d in file '%s', and it is reset to 1.\n", size,filename);
        size = 1;
      }
#endif /* !KIM */
/*added ends*/


      gradient = (int *)malloc(size * sizeof(int));
      invar_pot = (int *)malloc(size * sizeof(int));
#ifdef APOT
      smooth_pot = (int *)malloc(size * sizeof(int));
#endif /* APOT */
      for (i = 0; i < size; i++) {
	gradient[i] = 0;
	invar_pot[i] = 0;
#ifdef APOT
	smooth_pot[i] = 0;
#endif /* APOT */
      }
      have_format = 1;
    }

/* added */
/* KIM does not support T, I and G */
#ifndef KIM 

    /* header line with potential type */
    else if (buffer[1] == 'T') {
      if ((str = strchr(buffer + 3, '\n')) != NULL)
	*str = '\0';
      if (strcmp(buffer + 3, interaction_name) != 0) {
	error(0, "Wrong potential type found in potential file!\n");
	error(0, "This binary only supports %s potentials.\n", interaction_name);
	error(1, "Your potential file contains a %s potential.\n", buffer + 3);
      }
    }

    /* header line with invariant potentials */
    else if (buffer[1] == 'I') {
      if (have_format) {
#ifdef APOT
	apot_table.invar_pots = 0;
#endif /* APOT */
	/* gradient complete */
	for (i = 0; i < size; i++) {
	  str = strtok(((i == 0) ? buffer + 2 : NULL), " \t\r\n");
	  if (str == NULL) {
	    error(1, "Not enough items in #I header line.");
	  } else {
	    ((int *)invar_pot)[i] = atoi(str);
#ifdef APOT
	    apot_table.invar_pots++;
#endif /* APOT */
	  }
	}
	have_invar = 1;
      } else
	error(1, "#I needs to be specified after #F in file %s", filename);
    }
#ifndef APOT
    /* header line with gradients */
    else if (buffer[1] == 'G') {
      if (have_format) {
	/* gradient complete */
	for (i = 0; i < size; i++) {
	  str = strtok(((i == 0) ? buffer + 2 : NULL), " \t\r\n");
	  if (str == NULL)
	    error(1, "Not enough items in #G header line.");
	  else
	    ((int *)gradient)[i] = atoi(str);
	}
	have_grad = 1;
      } else
	error(1, "#G needs to be specified after #F in file %s", filename);
    }
#endif /* !APOT */

#endif /* !KIM */
/* added ends */

  } while (!end_header);

  /* do we have a format in the header? */
  if (!have_format)
    error(1, "Format not specified in header of potential file %s", filename);

  /* allocate info block of function table */
  pt->len = 0;
  pt->ncols = size;
  pt->begin = (double *)malloc(size * sizeof(double));
  pt->end = (double *)malloc(size * sizeof(double));
  pt->step = (double *)malloc(size * sizeof(double));
  pt->invstep = (double *)malloc(size * sizeof(double));
  pt->first = (int *)malloc(size * sizeof(int));
  pt->last = (int *)malloc(size * sizeof(int));
  if ((pt->begin == NULL) || (pt->end == NULL) || (pt->step == NULL)
    || (pt->invstep == NULL) || (pt->first == NULL) || (pt->last == NULL))
    error(1, "Cannot allocate info block for potential table %s", filename);
#ifdef APOT
  /* allocate memory for analytic potential table */
  apt->number = size;
  apt->total_par = 0;

  apt->n_par = (int *)malloc(size * sizeof(int));
  apt->begin = (double *)malloc(size * sizeof(double));
  apt->end = (double *)malloc(size * sizeof(double));
  apt->param_name = (char ***)malloc(size * sizeof(char **));
  apt->fvalue = (fvalue_pointer *) malloc(size * sizeof(fvalue_pointer));
#ifdef PAIR
  if (enable_cp) {
    apt->values = (double **)malloc((size + 1) * sizeof(double *));
    apt->values[size] = (double *)malloc(ntypes * sizeof(double));
    apt->invar_par = (int **)malloc(size * sizeof(int *));
    apt->chempot = apt->values[size];
    apt->pmin = (double **)malloc((size + 1) * sizeof(double *));
    apt->pmin[size] = (double *)malloc(ntypes * sizeof(double));
    apt->pmax = (double **)malloc((size + 1) * sizeof(double *));
    apt->pmax[size] = (double *)malloc(ntypes * sizeof(double));
  } else {
#elif defined COULOMB
  if (1) {
    apt->ratio = (double *)malloc(ntypes * sizeof(double));
    apt->values = (double **)malloc((size + 5) * sizeof(double *));
    apt->param_name = (char ***)malloc((size + 5) * sizeof(char **));
    apt->pmin = (double **)malloc((size + 5) * sizeof(double *));
    apt->pmax = (double **)malloc((size + 5) * sizeof(double *));
    apt->invar_par = (int **)malloc((size + 5) * sizeof(int *));

    apt->values[size] = (double *)malloc((ntypes - 1) * sizeof(double));
    apt->param_name[size] = (char **)malloc((ntypes - 1) * sizeof(char *));
    apt->pmin[size] = (double *)malloc((ntypes - 1) * sizeof(double));
    apt->pmax[size] = (double *)malloc((ntypes - 1) * sizeof(double));
    apt->invar_par[size] = (int *)malloc((ntypes - 1) * sizeof(int));

    apt->values[size + 1] = (double *)malloc(sizeof(double));
    apt->param_name[size + 1] = (char **)malloc(sizeof(char *));
    apt->pmin[size + 1] = (double *)malloc(sizeof(double));
    apt->pmax[size + 1] = (double *)malloc(sizeof(double));
    apt->invar_par[size + 1] = (int *)malloc(sizeof(int));

    apt->values[size + 2] = (double *)malloc(ntypes * sizeof(double));
    apt->param_name[size + 2] = (char **)malloc(ntypes * sizeof(char *));
    apt->pmin[size + 2] = (double *)malloc(ntypes * sizeof(double));
    apt->pmax[size + 2] = (double *)malloc(ntypes * sizeof(double));
    apt->invar_par[size + 2] = (int *)malloc(ntypes * sizeof(int));

    for (i = 3; i < 5; i++) {
      apt->values[size + i] = (double *)malloc(paircol * sizeof(double));
      apt->param_name[size + i] = (char **)malloc(paircol * sizeof(char *));
      apt->pmin[size + i] = (double *)malloc(paircol * sizeof(double));
      apt->pmax[size + i] = (double *)malloc(paircol * sizeof(double));
      apt->invar_par[size + i] = (int *)malloc(paircol * sizeof(int));
    }
    apt->charge = apt->values[size];
    apt->dp_kappa = apt->values[size + 1];
#ifdef DIPOLE
    apt->dp_alpha = apt->values[size + 2];
    apt->dp_b = apt->values[size + 3];
    apt->dp_c = apt->values[size + 4];
#endif /* DIPOLE */
  } else {
#endif /* COULOMB */
    apt->values = (double **)malloc(size * sizeof(double *));
    apt->invar_par = (int **)malloc(size * sizeof(int *));
    apt->pmin = (double **)malloc(size * sizeof(double *));
    apt->pmax = (double **)malloc(size * sizeof(double *));
#if defined PAIR || defined COULOMB
  }
#endif /* PAIR || COULOMB */

  apt->names = (char **)malloc(size * sizeof(char *));
  for (i = 0; i < size; i++)

/* added */
#ifndef KIM
    apt->names[i] = (char *)malloc(20 * sizeof(char));
#else
    /* KIM potential has longer name */
    apt->names[i] = (char *)malloc(255 * sizeof(char));
#endif /* !KIM */
/* added ends */

  if ((apt->n_par == NULL) || (apt->begin == NULL) || (apt->end == NULL)
    || (apt->fvalue == NULL) || (apt->names == NULL) || (apt->pmin == NULL)
    || (apt->pmax == NULL) || (apt->param_name == NULL)
    || (apt->values == NULL))
    error(1, "Cannot allocate info block for analytic potential table %s", filename);
#endif /* APOT */

  switch (format) {

/* added */
#ifndef KIM

#ifdef APOT
      case 0:
	read_pot_table0(pt, apt, filename, infile);
	break;
#else
      case 3:
	read_pot_table3(pt, size, filename, infile);
	break;
      case 4:
	read_pot_table4(pt, size, filename, infile);
#endif /* APOT */

#else /* !KIM */
      case 5:
#ifndef NOLIMITS
	read_pot_table5_no_nolimits(pt, apt, filename, infile);
	break;
#else /* !NOLIMITS */
	read_pot_table5_with_nolimits(pt, size, filename, infile);
	break;
#endif /* !NOLIMITS */

#endif /* !KIM */
/* added ends */
  }

  fclose(infile);

  printf("Reading potential file >> %s << ... done\n", filename);


/* added */
#ifndef KIM

  /* compute rcut and rmin */
  rcut = (double *)malloc(ntypes * ntypes * sizeof(double));
  if (NULL == rcut)
    error(1, "Cannot allocate rcut");
  rmin = (double *)malloc(ntypes * ntypes * sizeof(double));
  if (NULL == rmin)
    error(1, "Cannot allocate rmin");
  for (i = 0; i < ntypes; i++)
    for (j = 0; j < ntypes; j++) {
      k = (i <= j) ? i * ntypes + j - ((i * (i + 1)) / 2)
	: j * ntypes + i - ((j * (j + 1)) / 2);
      rmin[i * ntypes + j] = pt->begin[k];
      rcut[i * ntypes + j] = pt->end[k];
    }
#if defined EAM || defined ADP || defined MEAM
  for (i = 0; i < ntypes; i++) {
    for (j = 0; j < ntypes; j++) {
      rcut[i * ntypes + j] = MAX(rcut[i * ntypes + j], pt->end[(ntypes * (ntypes + 1)) / 2 + i]);
      rcut[i * ntypes + j] = MAX(rcut[i * ntypes + j], pt->end[(ntypes * (ntypes + 1)) / 2 + j]);
      rmin[i * ntypes + j] = MAX(rmin[i * ntypes + j], pt->begin[(ntypes * (ntypes + 1)) / 2 + i]);
      rmin[i * ntypes + j] = MAX(rmin[i * ntypes + j], pt->begin[(ntypes * (ntypes + 1)) / 2 + j]);
    }
  }
#endif /* EAM || ADP || MEAM */

#else /* !KIM */

		/* copy cutoff, cutoff equal pt->end[0], because there is only one potential */
		/* compute rcut and rmin */
		rcut = (double *)malloc(ntypes * ntypes * sizeof(double));
		if (NULL == rcut)
			error(1, "Cannot allocate rcut");
		rmin = (double *)malloc(ntypes * ntypes * sizeof(double));
		if (NULL == rmin)
			error(1, "Cannot allocate rmin");
		for (i = 0; i < ntypes; i++)
			for (j = 0; j < ntypes; j++) {
				rmin[i * ntypes + j] = 0.0;
				rcut[i * ntypes + j] = pt->end[0];
			}

#endif /* !KIM */
/* added ends */


#ifndef APOT
  double *val;

  /* read maximal changes file */
  maxchange = (double *)malloc(pt->len * sizeof(double));
  if (usemaxch) {
    /* open file */
    infile = fopen(maxchfile, "r");
    if (NULL == infile)
      error(1, "Could not open file %s\n", maxchfile);
    val = maxchange;
    for (i = 0; i < pt->len; i++) {
      if (1 > fscanf(infile, " %lf\n", val))
	error(1, "Premature end of maxch file %s", maxchfile);
      else
	val++;
    }
    fclose(infile);
  }
#endif /* !APOT */

  for (i = 0; i < ntypes; i++) {
    for (j = 0; j < ntypes; j++) {
      rcutmin = MIN(rcutmin, rcut[i + ntypes * j]);
      rcutmax = MAX(rcutmax, rcut[i + ntypes * j]);
    }
  }

  /* clean up locals and mark globals for later */
  reg_for_free(gradient, "gradient");
  reg_for_free(invar_pot, "invar_pot");
#ifdef APOT
  reg_for_free(smooth_pot, "smooth_pot");
  reg_for_free(apt->n_par, "apt->n_par");
  reg_for_free(apt->begin, "apt->begin");
  reg_for_free(apt->end, "apt->end");
  reg_for_free(apt->param_name, "apt->param_name");
  reg_for_free(apt->fvalue, "apt->fvalue");
  reg_for_free(apt->values, "apt->values");
  reg_for_free(apt->invar_par, "apt->invar_par");
  reg_for_free(apt->pmin, "apt->pmin");
  reg_for_free(apt->pmax, "apt->pmax");
  reg_for_free(apt->names, "apt->names");
  for (i = 0; i < size; i++) {
    reg_for_free(apt->names[i], "apt->names[%d]", i);
  }
#else /* APOT */
  reg_for_free(maxchange, "maxchange");
#endif /* APOT */
  reg_for_free(pt->begin, "pt->begin");
  reg_for_free(pt->end, "pt->end");
  reg_for_free(pt->step, "pt->step");
  reg_for_free(pt->invstep, "pt->invstep");
  reg_for_free(pt->first, "pt->first");
  reg_for_free(pt->last, "pt->last");
#if defined PAIR && defined APOT
  if (enable_cp) {
    reg_for_free(apt->chempot, "apt->chempot");
    reg_for_free(apt->pmin[size], "apt->pmin[%d]", size);
    reg_for_free(apt->pmax[size], "apt->pmax[%d]", size);
  }
#endif /* PAIR && APOT */
#ifdef COULOMB
  reg_for_free(apt->ratio, "apt->ratio");
  reg_for_free(apt->charge, "apt->charge");
  reg_for_free(apt->dp_kappa, "apt->dp_kappa");
  for (i = 0; i < 5; i++) {
    reg_for_free(apt->pmin[size + i], "apt->pmin[%d]", size + i);
    reg_for_free(apt->pmax[size + i], "apt->pmax[%d]", size + i);
    reg_for_free(apt->invar_par[size + i], "apt->invar_par[%d]", size + i);
  }
#endif /* COULOMB */
#ifdef DIPOLE
  reg_for_free(apt->dp_alpha, "apt->dp_alpha");
  reg_for_free(apt->dp_b, "apt->dp_b");
  reg_for_free(apt->dp_c, "apt->dp_c");
#endif /* DIPOLE */
  reg_for_free(rcut, "rcut");
  reg_for_free(rmin, "rmin");

  return;
}

/* added */
#ifndef KIM 

#ifdef APOT

/****************************************************************
 *
 *  read potential in analytic format:
 *  	for more information an how to specify an analytic potential
 *  	please check the documentation
 *
 *  parameters:
 *  	pot_table_t * ... pointer to the potential table
 *  	apot_table_t * ... pointer to the analytic potential table
 *  	char * ... name of the potential file (for error messages)
 *  	FILE * ... open file handle of the potential file
 *
 ****************************************************************/

void read_pot_table0(pot_table_t *pt, apot_table_t *apt, char *filename, FILE *infile)
{
  int   i, j, k, l, ret_val;
  char  buffer[255], name[255];
  char *token;
  double *val, *list, temp;
  fpos_t filepos, startpos;

  /* initialize the function table for analytic potentials */
  apot_init();

  /* save starting position */
  fgetpos(infile, &startpos);

#ifdef PAIR
  /* read cp */
  if (enable_cp) {

    /* search for cp */
    do {
      fgetpos(infile, &filepos);
      fscanf(infile, "%s", buffer);
    } while (strncmp(buffer, "cp", 2) != 0 && !feof(infile));
    /* and save the position */
    fsetpos(infile, &filepos);

    /* shortcut for apt->number */
    i = apt->number;

    /* allocate memory for global parameters */
    apt->names = (char **)realloc(apt->names, (i + 1) * sizeof(char *));
    apt->names[i] = (char *)malloc(20 * sizeof(char));
    strcpy(apt->names[i], "chemical potentials");
    reg_for_free(apt->names[i], "apt->names[%d] (chem_pot)", i);

    apt->invar_par = (int **)realloc(apt->invar_par, (i + 1) * sizeof(int *));
    apt->invar_par[i] = (int *)malloc((ntypes + 1) * sizeof(int));
    reg_for_free(apt->invar_par[i], "apt->invar_par[%d]", i);

    apt->param_name = (char ***)realloc(apt->param_name, (i + 1) * sizeof(char **));
    apt->param_name[i] = (char **)malloc(ntypes * sizeof(char *));
    reg_for_free(apt->param_name[i], "apt->param_name[%d]", i);

    /* check if the allocation was successfull */
    if (apt->names[i] == NULL || apt->invar_par[i] == NULL || apt->param_name[i] == NULL)
      error(1, "Cannot allocate memory for chemical potentials.");

    /* loop over all atom types */
    for (j = 0; j < ntypes; j++) {

      /* allocate memory for parameter name */
      apt->param_name[i][j] = (char *)malloc(30 * sizeof(char));
      reg_for_free(apt->param_name[i][j], "apt->param_name[%d][%d]", i, j);
      if (apt->param_name[i][j] == NULL)
	error(1, "Cannot allocate memory for chemical potential names.");

      /* read one line */
      if (4 > fscanf(infile, "%s %lf %lf %lf", buffer, &apt->chempot[j], &apt->pmin[i][j], &apt->pmax[i][j]))
	error(1, "Could not read chemical potential for %d. atomtype.", j);

      /* split cp and _# */
      token = strchr(buffer, '_');
      if (token != NULL) {
	strncpy(name, buffer, strlen(buffer) - strlen(token));
	name[strlen(buffer) - strlen(token)] = '\0';
      }
      if (strcmp("cp", name) != 0) {
	fprintf(stderr, "Found \"%s\" instead of \"cp\"\n", name);
	error(1, "No chemical potentials found in %s.\n", filename);
      }

      /* check for invariance and proper value (respect boundaries) */
      apt->invar_par[i][j] = 0.0;
      /* parameter will not be optimized if min==max */
      if (apt->pmin[i][j] == apt->pmax[i][j]) {
	apt->invar_par[i][j] = 1;
	apt->invar_par[i][ntypes]++;
	/* swap min and max if max<min */
      } else if (apt->pmin[i][j] > apt->pmax[i][j]) {
	temp = apt->pmin[i][j];
	apt->pmin[i][j] = apt->pmax[i][j];
	apt->pmax[i][j] = temp;
	/* reset value if >max or <min */
      } else if ((apt->values[i][j] < apt->pmin[i][j])
	|| (apt->values[i][j] > apt->pmax[i][j])) {
	/* Only print warning if we are optimizing */
	if (opt) {
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
    printf(" - Enabled %d chemical potential(s)\n", ntypes);

    /* disable composition nodes for now */
#ifdef CN
    /* read composition nodes */
    if (2 > fscanf(infile, "%s %d", buffer, &compnodes)) {
      if (strcmp("type", buffer) == 0)
	compnodes = -1;
      else
	error(1, "Could not read number of composition nodes from potential file.\n");
    }
    if (strcmp(buffer, "cn") != 0 && ntypes > 1 && compnodes != -1)
      error(1, "No composition nodes found in %s.\nUse \"cn 0\" for none.\n", filename);
    if (ntypes == 1) {
      compnodes = 0;
    }
    if (compnodes != -1) {
      apt->values[apt->number] =
	(double *)realloc(apt->values[apt->number], (ntypes + compnodes) * sizeof(double));
      apt->pmin[apt->number] =
	(double *)realloc(apt->pmin[apt->number], (ntypes + compnodes) * sizeof(double));
      apt->pmax[apt->number] =
	(double *)realloc(apt->pmax[apt->number], (ntypes + compnodes) * sizeof(double));
      apt->chempot = apt->values[apt->number];
      compnodelist = (double *)malloc((ntypes + compnodes) * sizeof(double));

      for (j = 0; j < compnodes; j++) {
	if (4 > fscanf(infile, "%lf %lf %lf %lf", &compnodelist[j],
	    &apt->chempot[ntypes + j], &apt->pmin[apt->number][ntypes + j],
	    &apt->pmax[apt->number][ntypes + j]))
	  error(1, "Could not read composition node %d\n", j + 1);
	if (apt->pmin[apt->number][ntypes + j] > apt->chempot[ntypes + j]
	  || apt->pmax[apt->number][ntypes + j] < apt->chempot[ntypes + j])
	  error(1, "composition node %d is out of bounds.\n", j + 1);
      }

      /* check compnodes for valid values */
      if (ntypes == 2) {
	for (j = 0; j < compnodes; j++)
	  if (compnodelist[j] > 1 || compnodelist[j] < 0)
	    error(1, "Composition node %d is %f but should be inside [0,1].\n", j + 1, compnodelist[j]);
      }
    }
    if (compnodes != -1)
      printf("Enabled chemical potentials with %d extra composition node(s).\n", compnodes);
    if (compnodes == -1)
      compnodes = 0;
#endif /* CN */
  }
#endif /* PAIR */

#ifdef COULOMB
  fsetpos(infile, &startpos);
  /* skip to electrostatic section */
  do {
    fgetpos(infile, &filepos);
    fscanf(infile, "%s", buffer);
  } while (strcmp(buffer, "elstat") != 0 && !feof(infile));

  /* check for elstat keyword */
  if (strcmp("elstat", buffer) != 0) {
    error(1, "No elstat option found in %s.\n", filename);
  }

  /* read electrostatic parameters */
  fscanf(infile, " %s", buffer);
  if (strcmp("ratio", buffer) != 0) {
    error(1, "Could not read ratio");
  }
  for (i = 0; i < ntypes; i++) {
    if (1 > fscanf(infile, "%lf", &apt->ratio[i])) {
      error(1, "Could not read ratio for atomtype #%d\n", i);
    }
  }
  for (i = 0; i < ntypes - 1; i++) {
    apt->param_name[apt->number][i] = (char *)malloc(30 * sizeof(char));
    if (4 > fscanf(infile, "%s %lf %lf %lf", apt->param_name[apt->number][i],
	&apt->charge[i], &apt->pmin[apt->number][i], &apt->pmax[apt->number][i])) {
      error(1, "Could not read charge for atomtype #%d\n", i);
    }
    apt->invar_par[apt->number][i] = 0;
    if (apt->pmin[apt->number][i] == apt->pmax[apt->number][i]) {
      apt->invar_par[apt->number][i]++;
    }
    reg_for_free(apt->param_name[apt->number][i], "apt->param_name[%d][%d]", apt->number, i);
  }
  apt->param_name[apt->number + 1][0] = (char *)malloc(30 * sizeof(char));
  if (4 > fscanf(infile, "%s %lf %lf %lf", apt->param_name[apt->number + 1][0],
      &apt->dp_kappa[0], &apt->pmin[apt->number + 1][0], &apt->pmax[apt->number + 1][0])) {
    error(1, "Could not read kappa");
  }
  apt->invar_par[apt->number + 1][0] = 0;
  if (apt->pmin[apt->number + 1][0] == apt->pmax[apt->number + 1][0]) {
    apt->invar_par[apt->number + 1][0]++;
  }
  reg_for_free(apt->param_name[apt->number + 1][0], "apt->param_name[%d][%d]", apt->number + 1, 0);
  apt->sw_kappa = apt->invar_par[apt->number + 1][0];
#ifndef DIPOLE
  printf(" - Read elstat table\n");
#endif /* !DIPOLE */
#endif /* COULOMB */

#ifdef DIPOLE
  int   ncols = ntypes * (ntypes + 1) / 2;

  for (i = 0; i < ntypes; i++) {
    apt->param_name[apt->number + 2][i] = (char *)malloc(30 * sizeof(char));
    if (4 > fscanf(infile, "%s %lf %lf %lf",
	apt->param_name[apt->number + 2][i], &apt->dp_alpha[i],
	&apt->pmin[apt->number + 2][i], &apt->pmax[apt->number + 2][i])) {
      error(1, "Could not read polarisability for atomtype #%d\n", i);
    }
    apt->invar_par[apt->number + 2][i] = 0;
    if (apt->pmin[apt->number + 2][i] == apt->pmax[apt->number + 2][i]) {
      apt->invar_par[apt->number + 2][i]++;
    }
    reg_for_free(apt->param_name[apt->number + 2][i], "apt->param_name[%d][%d]", apt->number + 2, i);
  }
  for (i = 0; i < ncols; i++) {
    apt->param_name[apt->number + 3][i] = (char *)malloc(30 * sizeof(char));
    if (4 > fscanf(infile, "%s %lf %lf %lf",
	apt->param_name[apt->number + 3][i], &apt->dp_b[i],
	&apt->pmin[apt->number + 3][i], &apt->pmax[apt->number + 3][i])) {
      error(1, "Could not read parameter dp_b for potential #%d\n", i);
    }
    apt->invar_par[apt->number + 3][i] = 0;
    if (apt->pmin[apt->number + 3][i] == apt->pmax[apt->number + 3][i]) {
      apt->invar_par[apt->number + 3][i]++;
    }
    reg_for_free(apt->param_name[apt->number + 3][i], "apt->param_name[%d][%d]", apt->number + 3, i);
  }
  for (i = 0; i < ncols; i++) {
    apt->param_name[apt->number + 4][i] = (char *)malloc(30 * sizeof(char));
    if (4 > fscanf(infile, "%s %lf %lf %lf",
	apt->param_name[apt->number + 4][i], &apt->dp_c[i],
	&apt->pmin[apt->number + 4][i], &apt->pmax[apt->number + 4][i])) {
      error(1, "Could not read parameter dp_c for potential #%d\n", i);
    }
    apt->invar_par[apt->number + 4][i] = 0;
    if (apt->pmin[apt->number + 4][i] == apt->pmax[apt->number + 4][i]) {
      apt->invar_par[apt->number + 4][i]++;
    }
    reg_for_free(apt->param_name[apt->number + 4][i], "apt->param_name[%d][%d]", apt->number + 4, i);
  }

  printf(" - Read elstat table\n");
#endif /* DIPOLE */

  /* skip to global section */
  fsetpos(infile, &startpos);
  do {
    fgetpos(infile, &filepos);
    fscanf(infile, "%s", buffer);
  } while (strcmp(buffer, "global") != 0 && !feof(infile));
  fsetpos(infile, &filepos);

  /* check for global keyword */
  if (strcmp(buffer, "global") == 0) {
    if (2 > fscanf(infile, "%s %d", buffer, &apt->globals))
      error(1, "Premature end of potential file %s", filename);
    have_globals = 1;
    apt->total_par += apt->globals;

    i = apt->number + enable_cp;
    j = apt->globals;
    global_pot = i;

    /* allocate memory for global parameters */
    apt->names = (char **)realloc(apt->names, (global_pot + 1) * sizeof(char *));
    apt->names[global_pot] = (char *)malloc(20 * sizeof(char));
    strcpy(apt->names[global_pot], "global parameters");
    reg_for_free(apt->names[global_pot], "apt->names[%d] (global_pot)", global_pot);

    apt->n_glob = (int *)malloc(apt->globals * sizeof(int));
    reg_for_free(apt->n_glob, "apt->n_glob");

    apt->global_idx = (int ***)malloc(apt->globals * sizeof(int **));
    reg_for_free(apt->global_idx, "apt->global_idx");

    apt->values = (double **)realloc(apt->values, (global_pot + 1) * sizeof(double *));
    apt->values[global_pot] = (double *)malloc(j * sizeof(double));
    reg_for_free(apt->values[global_pot], "apt->values[%d]", global_pot);

    apt->invar_par = (int **)realloc(apt->invar_par, (global_pot + 1) * sizeof(int *));
    apt->invar_par[global_pot] = (int *)malloc((j + 1) * sizeof(int));
    reg_for_free(apt->invar_par[global_pot], "apt->invar_par[%d]", global_pot);

    apt->pmin = (double **)realloc(apt->pmin, (global_pot + 1) * sizeof(double *));
    apt->pmin[global_pot] = (double *)malloc(j * sizeof(double));
    reg_for_free(apt->pmin[global_pot], "apt->pmin[%d]", global_pot);

    apt->pmax = (double **)realloc(apt->pmax, (global_pot + 1) * sizeof(double *));
    apt->pmax[global_pot] = (double *)malloc(j * sizeof(double));
    reg_for_free(apt->pmax[global_pot], "apt->pmax[%d]", global_pot);

    apt->param_name = (char ***)realloc(apt->param_name, (global_pot + 1) * sizeof(char **));
    apt->param_name[global_pot] = (char **)malloc(j * sizeof(char *));
    reg_for_free(apt->param_name[global_pot], "apt->param_name[%d]", global_pot);

    pt->first = (int *)realloc(pt->first, (global_pot + 1) * sizeof(int));

    if (NULL == apt->values[global_pot] || NULL == apt->pmin[global_pot]
      || NULL == apt->n_glob || NULL == apt->global_idx || NULL == apt->pmax[global_pot]
      || NULL == apt->param_name[global_pot])
      error(1, "Cannot allocate memory for global paramters.");

    /* initialize properly */
    apt->invar_par[i][j] = 0;

    /* read the global parameters */
    for (j = 0; j < apt->globals; j++) {
      apt->param_name[global_pot][j] = (char *)malloc(30 * sizeof(char));
      reg_for_free(apt->param_name[global_pot][j], "apt->param_name[%d][%d]", global_pot, j);

      if (NULL == apt->param_name[global_pot][j])
	error(1, "Error in allocating memory for global parameter name");

      strcpy(apt->param_name[global_pot][j], "\0");
      ret_val =
	fscanf(infile, "%s %lf %lf %lf", apt->param_name[global_pot][j],
	&apt->values[global_pot][j], &apt->pmin[global_pot][j], &apt->pmax[global_pot][j]);
      if (4 > ret_val)
	if (strcmp(apt->param_name[global_pot][j], "type") == 0) {
	  error(0, "Not enough global parameters!\n");
	  error(1, "You specified %d parameter(s), but needed are %d.\nAborting", j, apt->globals);
	}

      /* check for duplicate names */
      for (k = j - 1; k >= 0; k--)
	if (strcmp(apt->param_name[global_pot][j], apt->param_name[global_pot][k]) == 0) {
	  error(0, "\nFound duplicate global parameter name!\n");
	  error(1, "Parameter #%d (%s) is the same as #%d (%s)\n", j + 1,
	    apt->param_name[global_pot][j], k + 1, apt->param_name[global_pot][k]);
	}
      apt->n_glob[j] = 0;

      /* check for invariance and proper value (respect boundaries) */
      /* parameter will not be optimized if min==max */
      apt->invar_par[i][j] = 0;
      if (apt->pmin[i][j] == apt->pmax[i][j]) {
	apt->invar_par[i][j] = 1;
	apt->invar_par[i][apt->globals]++;
      } else if (apt->pmin[i][j] > apt->pmax[i][j]) {
	temp = apt->pmin[i][j];
	apt->pmin[i][j] = apt->pmax[i][j];
	apt->pmax[i][j] = temp;
      } else if ((apt->values[i][j] < apt->pmin[i][j])
	|| (apt->values[i][j] > apt->pmax[i][j])) {
	/* Only print warning if we are optimizing */
	if (opt) {
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

  /* skip to actual potentials */
  fsetpos(infile, &startpos);
  do {
    fgetpos(infile, &filepos);
    fscanf(infile, "%s", buffer);
  } while (strcmp(buffer, "type") != 0 && !feof(infile));
  fsetpos(infile, &filepos);

  for (i = 0; i < apt->number; i++) {
    /* scan for "type" keyword */
    do {
      fgetpos(infile, &filepos);
      fscanf(infile, "%s", buffer);
    } while (strcmp(buffer, "type") != 0 && !feof(infile));
    fsetpos(infile, &filepos);
    /* read type */
    if (2 > fscanf(infile, "%s %s", buffer, name))
      error(1, "Premature end of potential file %s", filename);
    if (strcmp(buffer, "type") != 0)
      error(1, "Unknown keyword in file %s, expected \"type\" but found \"%s\".", filename, buffer);

    /* split name and _sc */
    token = strrchr(name, '_');
    if (token != NULL && strcmp(token + 1, "sc") == 0) {
      strncpy(buffer, name, strlen(name) - 3);
      buffer[strlen(name) - 3] = '\0';
      strcpy(name, buffer);
      smooth_pot[i] = 1;
    }

    /* check if potential is "pohlong" and change it to bjs */
    if (strcmp(name, "pohlong") == 0)
      strcpy(name, "bjs\0");

    if (apot_parameters(name) == -1)
      error(1, "Unknown function type in file %s, please define \"%s\" in functions.c.", filename, name);

    strcpy(apt->names[i], name);
    apt->n_par[i] = apot_parameters(name);

    /* add one parameter for cutoff function if _sc is found */
    if (smooth_pot[i] == 1)
      apt->n_par[i]++;
    apt->total_par += apt->n_par[i];

    /* read cutoff */
    if (2 > fscanf(infile, "%s %lf", buffer, &apt->end[i]))
      error(1, "Could not read cutoff for potential #%d in file %s\nAborting", i, filename);
    if (strcmp(buffer, "cutoff") != 0)
      error(1,
	"No cutoff found for the %d. potential (%s) after \"type\" in file %s.\nAborting",
	i + 1, apt->names[i], filename);

    /* set very small begin, needed for EAM embedding function */
    apt->begin[i] = 0.0001;

    /* allocate memory for this parameter */
    apt->values[i] = (double *)malloc(apt->n_par[i] * sizeof(double));
    reg_for_free(apt->values[i], "apt->values[%d]", i);

    apt->invar_par[i] = (int *)malloc((apt->n_par[i] + 1) * sizeof(int));
    reg_for_free(apt->invar_par[i], "apt->invar_par[%d]", i);

    apt->pmin[i] = (double *)malloc(apt->n_par[i] * sizeof(double));
    reg_for_free(apt->pmin[i], "apt->pmin[%d]", i);

    apt->pmax[i] = (double *)malloc(apt->n_par[i] * sizeof(double));
    reg_for_free(apt->pmax[i], "apt->pmax[%d]", i);

    apt->param_name[i] = (char **)malloc(apt->n_par[i] * sizeof(char *));
    reg_for_free(apt->param_name[i], "apt->param_name[%d]", i);

    if (NULL == apt->values[i] || NULL == apt->pmin[i]
      || NULL == apt->pmax[i] || NULL == apt->param_name[i])
      error(1, "Cannot allocate memory for potential paramters.\nAborting");

    /* check for comments */
    do {
      j = fgetc(infile);
    } while (j != 10);

    fgetpos(infile, &filepos);
    fgets(buffer, 255, infile);
    while (buffer[0] == '#') {
      fgetpos(infile, &filepos);
      fgets(buffer, 255, infile);
    }
    fsetpos(infile, &filepos);

    /* read parameters */
    apt->invar_par[i][apt->n_par[i]] = 0;
    for (j = 0; j < apt->n_par[i]; j++) {
      apt->param_name[i][j] = (char *)malloc(30 * sizeof(char));
      if (NULL == apt->param_name[i][j])
	error(1, "Error in allocating memory for parameter name");
      reg_for_free(apt->param_name[i][j], "apt->param_name[%d][%d]", i, j);
      strcpy(apt->param_name[i][j], "empty");
      fgetpos(infile, &filepos);
      fgets(name, 255, infile);
      while (name[0] == '#' && !feof(infile) && (j != apt->n_par[i] - 1)) {
	fgets(name, 255, infile);
      }
      if ((j != (apt->n_par[i] - 1)) && (feof(infile) || name[0] == '\0')) {
	error(0, "Premature end of potential definition or file.\n");
	error(1, "Probably your potential definition is missing some parameters.\n");
      }
      if (feof(infile))
	name[0] = '\0';
      buffer[0] = '\0';
      ret_val =
	sscanf(name, "%s %lf %lf %lf", buffer, &apt->values[i][j], &apt->pmin[i][j], &apt->pmax[i][j]);
      if (buffer[0] != '\0')
	strncpy(apt->param_name[i][j], buffer, 30);

      /* if last char of name is "!" we have a global parameter */
      if (strrchr(apt->param_name[i][j], '!') != NULL) {
	apt->param_name[i][j][strlen(apt->param_name[i][j]) - 1] = '\0';
	l = -1;
	for (k = 0; k < apt->globals; k++)
	  if (strcmp(apt->param_name[i][j], apt->param_name[global_pot][k]) == 0)
	    l = k;
	if (-1 == l)
	  error(1, "Could not find global parameter %s!\n", apt->param_name[i][j]);
	sprintf(apt->param_name[i][j], "%s!", apt->param_name[i][j]);

	/* write index array for global parameters */
	if (++apt->n_glob[l] > 1) {
	  apt->global_idx[l] = (int **)realloc(apt->global_idx[l], apt->n_glob[l] * sizeof(int *));
	} else {
	  apt->global_idx[l] = (int **)malloc(1 * sizeof(int *));
	}
	apt->global_idx[l][apt->n_glob[l] - 1] = (int *)malloc(2 * sizeof(int));
	apt->global_idx[l][apt->n_glob[l] - 1][0] = i;
	apt->global_idx[l][apt->n_glob[l] - 1][1] = j;

	apt->values[i][j] = apt->values[global_pot][l];
	apt->pmin[i][j] = apt->pmin[global_pot][l];
	apt->pmax[i][j] = apt->pmax[global_pot][l];
	apt->invar_par[i][j] = 1;
	apt->invar_par[i][apt->n_par[i]]++;
      } else {
	/* this is no global parameter */
	if (4 > ret_val) {
	  if (smooth_pot[i] && j == apot_parameters(apt->names[i])) {
	    if (0 == strcmp(apt->param_name[i][j], "type") ||
	      0 == strcmp(apt->param_name[i][j], "empty") || feof(infile)) {
	      warning("No cutoff parameter given for potential #%d: adding one parameter.\n", i);
	      strcpy(apt->param_name[i][j], "h");
	      apt->values[i][j] = 1;
	      apt->pmin[i][j] = 0.5;
	      apt->pmax[i][j] = 2;
	      fsetpos(infile, &filepos);
	    }
	  } else {
	    if (strcmp(apt->param_name[i][j], "type") == 0 || ret_val == EOF) {
	      error(0, "Not enough parameters for potential #%d (%s) in file %s!\n", i + 1, apt->names[i],
		filename);
	      error(1, "You specified %d parameter(s), but required are %d.\n", j, apt->n_par[i]);
	    }
	    error(1, "Could not read parameter #%d of potential #%d in file %s", j + 1, i + 1, filename);
	  }
	}

	/* check for invariance and proper value (respect boundaries) */
	/* parameter will not be optimized if min==max */
	apt->invar_par[i][j] = 0;
	if (apt->pmin[i][j] == apt->pmax[i][j]) {
	  apt->invar_par[i][j] = 1;
	  apt->invar_par[i][apt->n_par[i]]++;
	} else if (apt->pmin[i][j] > apt->pmax[i][j]) {
	  temp = apt->pmin[i][j];
	  apt->pmin[i][j] = apt->pmax[i][j];
	  apt->pmax[i][j] = temp;
	} else if ((apt->values[i][j] < apt->pmin[i][j])
	  || (apt->values[i][j] > apt->pmax[i][j])) {
	  /* Only print warning if we are optimizing */
	  if (opt) {
	    if (apt->values[i][j] < apt->pmin[i][j])
	      apt->values[i][j] = apt->pmin[i][j];
	    if (apt->values[i][j] > apt->pmax[i][j])
	      apt->values[i][j] = apt->pmax[i][j];
	    warning("Starting value for parameter #%d in potential #%d is ", j + 1, i + 1);
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

  /* clean up globals table */
  if (have_globals) {
    for (i = 0; i < apt->globals; i++) {
      reg_for_free(apt->global_idx[i], "apt->global_idx[%d]", i);
      for (j = 0; j < apt->n_glob[i]; j++)
	reg_for_free(apt->global_idx[i][j], "apt->global_idx[%d][%d]", i, j);
    }
  }
#ifdef COULOMB
  apt->total_ne_par = apt->total_par;
#endif /* COULOMB */

  /* if we have global parameters, are they actually used ? */
  if (have_globals) {
    j = 0;
    for (i = 0; i < apt->globals; i++)
      j += apt->n_glob[i];
    if (j == 0) {
      have_globals = 0;
      printf("You defined global parameters but did not use them.\n");
      printf("Disabling global parameters.\n\n");
    }
  }

  /* assign the potential functions to the function pointers */
  if (apot_assign_functions(apt) == -1)
    error(1, "Could not assign the function pointers.\n");
#ifdef PAIR
  if (enable_cp) {
    cp_start = apt->total_par - apt->globals + ntypes * (ntypes + 1);
    apt->total_par += (ntypes + compnodes - apt->invar_par[apt->number][ntypes]);
  }
#endif /* PAIR */

#ifdef COULOMB
  apt->total_par += ntypes;
#endif /* COULOMB */
#ifdef DIPOLE
  apt->total_par += ntypes;
  apt->total_par += (2 * ncols);
#endif /* DIPOLE */

  /* initialize function table and write indirect index */
  for (i = 0; i < apt->number; i++) {
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
  if (have_globals)
    pt->len += apt->globals;

#ifdef PAIR
  if (enable_cp) {
    pt->len += (ntypes + compnodes);
  }
#endif /* PAIR */
#ifdef COULOMB
  pt->len += ntypes;
  pt->len += ntypes - 1;
#endif /* COULOMB */
#ifdef DIPOLE
  pt->len += ntypes;
  pt->len += (2 * ncols);
#endif /* DIPOLE */

  pt->table = (double *)malloc(pt->len * sizeof(double));
  reg_for_free(pt->table, "pt->table");

  calc_list = (double *)malloc(pt->len * sizeof(double));
  reg_for_free(calc_list, "calc_list");

  pt->idx = (int *)malloc(pt->len * sizeof(int));
  reg_for_free(pt->idx, "pt->idx");

  apt->idxpot = (int *)malloc(apt->total_par * sizeof(int));
  reg_for_free(apt->idxpot, "apt->idxpot");

  apt->idxparam = (int *)malloc(apt->total_par * sizeof(int));
  reg_for_free(apt->idxparam, "apt->idxparam");

  if ((NULL == pt->table) || (NULL == pt->idx) || (apt->idxpot == NULL)
    || (apt->idxparam == NULL))
    error(1, "Cannot allocate memory for potential table.\n");
  for (i = 0; i < pt->len; i++) {
    pt->table[i] = 0;
    calc_list[i] = 0;
    pt->idx[i] = 0;
  }

  /* this is the indirect index */
  k = 0;
  l = 0;
  val = pt->table;
  list = calc_list;
  for (i = 0; i < apt->number; i++) {	/* loop over potentials */
    val += 2;
    list += 2;
    l += 2;
    for (j = 0; j < apt->n_par[i]; j++) {	/* loop over parameters */
      *val = apt->values[i][j];
      *list = apt->values[i][j];
      val++;
      list++;
      if (!invar_pot[i] && !apt->invar_par[i][j]) {
	pt->idx[k] = l++;
	apt->idxpot[k] = i;
	apt->idxparam[k++] = j;
      } else
	l++;
    }
    if (!invar_pot[i])
      pt->idxlen += apt->n_par[i] - apt->invar_par[i][apt->n_par[i]];
    apt->total_par -= apt->invar_par[i][apt->n_par[i]];
  }

#ifdef PAIR
  if (enable_cp) {
    init_chemical_potential(ntypes);
    i = apt->number;
    for (j = 0; j < (ntypes + compnodes); j++) {
      *val = apt->values[i][j];
      val++;
      if (!apt->invar_par[i][j]) {
	pt->idx[k] = l++;
	apt->idxpot[k] = i;
	apt->idxparam[k++] = j;
      }
    }
    pt->idxlen += (ntypes + compnodes - apt->invar_par[apt->number][ntypes]);
    global_idx += (ntypes + compnodes - apt->invar_par[apt->number][ntypes]);
  }
#endif /* PAIR */

#ifdef COULOMB
  i = apt->number;
  for (j = 0; j < (ntypes - 1); j++) {
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
  pt->idxlen += ntypes;
#endif /* COULOMB */

#ifdef DIPOLE
  i = apt->number + 2;
  for (j = 0; j < (ntypes); j++) {
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
    for (j = 0; j < (ncols); j++) {
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
  pt->idxlen += ntypes;
  pt->idxlen += (2 * ncols);
#endif /* DIPOLE */

  if (have_globals) {
    i = global_pot;
    for (j = 0; j < apt->globals; j++) {
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
  global_idx += pt->last[apt->number - 1] + 1;

#ifdef NOPUNISH
  if (opt)
    warning("Gauge degrees of freedom are NOT fixed!\n");
#endif /* NOPUNISH */

  check_apot_functions();

  init_calc_table(pt, &calc_pot);

  return;
}

#else /* APOT */

/****************************************************************
 *
 *  read potential in third format:
 *
 *  Sampling points are equidistant.
 *
 *  Header:  one line for each function with
 *           rbegin rstart npoints
 *
 *  Table: Function values at sampling points,
 *         functions separated by blank lines
 *
 *  parameters:
 * 	pot_table_t * ... pointer to the potential table
 *  	int ... number of potential functions
 *  	char * ... name of the potential file (for error messages)
 *  	FILE * ... open file handle of the potential file
 *
 ****************************************************************/

void read_pot_table3(pot_table_t *pt, int size, char *filename, FILE *infile)
{
  int   i, j, k, l;
  int   nvals[size];
  double *val;

  /* read the info block of the function table */
  for (i = 0; i < size; i++) {
    if (3 > fscanf(infile, "%lf %lf %d", &pt->begin[i], &pt->end[i], &nvals[i]))
      error(1, "Premature end of potential file %s", filename);
    pt->step[i] = (pt->end[i] - pt->begin[i]) / (nvals[i] - 1);
    pt->invstep[i] = 1.0 / pt->step[i];
    /* in the two slots between last[i-1] and first[i] the gradients
       of the respective functions are stored */
    if (i == 0)
      pt->first[i] = 2;
    else
      pt->first[i] = pt->last[i - 1] + 3;
    pt->last[i] = pt->first[i] + nvals[i] - 1;
    pt->len = pt->first[i] + nvals[i];
  }
  /* allocate the function table */
  pt->table = (double *)malloc(pt->len * sizeof(double));
  reg_for_free(pt->table, "pt->table");
  pt->xcoord = (double *)malloc(pt->len * sizeof(double));
  reg_for_free(pt->xcoord, "pt->xcoord");
  pt->d2tab = (double *)malloc(pt->len * sizeof(double));
  reg_for_free(pt->d2tab, "pt->d2tab");
  pt->idx = (int *)malloc(pt->len * sizeof(int));
  reg_for_free(pt->idx, "pt->idx");
  for (i = 0; i < pt->len; i++) {
    pt->table[i] = 0.0;
    pt->xcoord[i] = 0.0;
    pt->d2tab[i] = 0.0;
    pt->idx[i] = 0.0;
  }
  if ((NULL == pt->table) || (NULL == pt->idx) || (NULL == pt->d2tab))
    error(1, "Cannot allocate memory for potential table");

  /* input loop */
  val = pt->table;
  k = 0;
  l = 0;

  /* read pair potentials */
  for (i = 0; i < paircol; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.0;
    }
    val += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(infile, "%lf\n", val))
	error(1, "Premature end of potential file %s", filename);
      else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;
    }
  }

#if defined EAM || defined ADP || defined MEAM
  /* read EAM transfer function rho(r) */
  for (i = paircol; i < paircol + ntypes; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.0;
    }
    val += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(infile, "%lf\n", val))
	error(1, "Premature end of potential file %s", filename);
      else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;
    }
  }

  /* read EAM embedding function F(n) */
  for (i = paircol + ntypes; i < paircol + 2 * ntypes; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);
    } else {
      *val = 1.e30;
      *(val + 1) = 1.e30;
    }
    val += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(infile, "%lf\n", val))
	error(1, "Premature end of potential file %s", filename);
      else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if (!invar_pot[i])
	pt->idx[k++] = l++;
      else
	l++;
    }
  }

#ifdef TBEAM
  /* read TBEAM transfer function rho(r) for the s-band */
  for (i = paircol + 2 * ntypes; i < paircol + 3 * ntypes; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.0;
    }
    val += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(infile, "%lf\n", val))
	error(1, "Premature end of potential file %s", filename);
      else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;
    }
  }

  /* read TBEAM embedding function F(n) for the s-band */
  for (i = paircol + 3 * ntypes; i < paircol + 4 * ntypes; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);
    } else {
      *val = 1.e30;
      *(val + 1) = 1.e30;
    }
    val += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(infile, "%lf\n", val))
	error(1, "Premature end of potential file %s", filename);
      else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if (!invar_pot[i])
	pt->idx[k++] = l++;
      else
	l++;
    }
  }
#endif /* TBEAM */
#endif /* EAM || ADP || MEAM */

#ifdef ADP
  /* read ADP dipole function u(r) */
  for (i = paircol + 2 * ntypes; i < 2 * (paircol + ntypes); i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.0;
    }
    val += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(infile, "%lf\n", val))
	error(1, "Premature end of potential file %s", filename);
      else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;
    }
  }

  /* read adp quadrupole function w(r) */
  for (i = 2 * (paircol + ntypes); i < 3 * paircol + 2 * ntypes; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);
    } else {
      *val = 1.e30;
      *(val + 1) = 1.e30;
    }
    val += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(infile, "%lf\n", val))
	error(1, "Premature end of potential file %s", filename);
      else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;
    }
  }
#endif /* ADP */

#ifdef MEAM
  for (i = paircol + 2 * ntypes; i < 2 * paircol + 2 * ntypes; i++) {	/* read in second pair pot    f */
    if (have_grad) {		/* read gradient */
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);;
    } else {
      *val = 1e30;
      *(val + 1) = 0.0;
    }
    val += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    for (j = 0; j < nvals[i]; j++) {	/* read values */
      if (1 > fscanf(infile, "%lf\n", val)) {
	error(1, "Premature end of potential file %s", filename);;
      } else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      // Clamp first spline knot in first f_ij potential only
      // to remove degeneracy of f*f*g where f' = f/b and g' = b^2*g
#ifndef MEAMf
      if ((!invar_pot[i]) && (j < nvals[i] - 1 && (j != 0 || i != paircol + 2 * ntypes)))
#else
      if (!invar_pot[i])
#endif /* MEAMf */
	pt->idx[k++] = l++;
      else
	l++;
    }
  }
  for (i = 2 * paircol + 2 * ntypes; i < 2 * paircol + 3 * ntypes; i++) {	/* read in angl part */
    if (have_grad) {		/* read gradient */
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);;
    } else {
      *val = 0;
      *(val + 1) = 0;
    }
    val += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    for (j = 0; j < nvals[i]; j++) {	/* read values */
      if (1 > fscanf(infile, "%lf\n", val)) {
	error(1, "Premature end of potential file %s", filename);;
      } else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if (!invar_pot[i])
	pt->idx[k++] = l++;
      else
	l++;
    }
  }
#endif /* MEAM */

  pt->idxlen = k;
  init_calc_table(pt, &calc_pot);

  return;
}

/****************************************************************
 *
 *  read potential in fourth format:
 *
 *  Sampling points are NON-equidistant.
 *
 *  Header:  one line for each function with
 *           npoints
 *
 *  Table: Sampling points, function values
 *            r f(r)
 *         functions separated by blank lines
 *
 ****************************************************************/

void read_pot_table4(pot_table_t *pt, int size, char *filename, FILE *infile)
{
  int   i, k, l, j;
  int   den_count;
  int   emb_count;
  int   nvals[size];
  double *val, *ord;

  /* read the info block of the function table */
  for (i = 0; i < size; i++) {
    if (1 > fscanf(infile, "%d", &nvals[i]))
      error(1, "Premature end of potential file %s", filename);
    pt->step[i] = 0.0;
    pt->invstep[i] = 0.0;
    if (i == 0)
      pt->first[i] = 2;
    else
      pt->first[i] = pt->last[i - 1] + 3;
    pt->last[i] = pt->first[i] + nvals[i] - 1;
    pt->len = pt->first[i] + nvals[i];
  }
  /* allocate the function table */
  pt->table = (double *)malloc(pt->len * sizeof(double));
  pt->xcoord = (double *)malloc(pt->len * sizeof(double));
  pt->d2tab = (double *)malloc(pt->len * sizeof(double));
  pt->idx = (int *)malloc(pt->len * sizeof(int));
  if ((NULL == pt->table) || (NULL == pt->xcoord) || (NULL == pt->idx)
    || (NULL == pt->d2tab))
    error(1, "Cannot allocate memory for potential table");

  for (i = 0; i < pt->len; i++) {
    pt->table[i] = 0.0;
    pt->xcoord[i] = 0.0;
    pt->d2tab[i] = 0.0;
    pt->idx[i] = 0;
  }

  /* input loop */
  val = pt->table;
  ord = pt->xcoord;
  k = 0;
  l = 0;

  /* read pair potentials */
  for (i = 0; i < paircol; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.0;
    }
    val += 2;
    ord += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (j = 0; j < nvals[i]; j++) {
      if (2 > fscanf(infile, "%lf %lf\n", ord, val))
	error(1, "Premature end of potential file %s", filename);
      else {
	val++;
	ord++;
      }
      if ((j > 0) && (*(ord - 1) <= *(ord - 2)))
	error(1, "Abscissa not monotonous in potential %d.", i);
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;

    }
    pt->begin[i] = pt->xcoord[pt->first[i]];
    pt->end[i] = pt->xcoord[pt->last[i]];
    /* pt->step is average step length.. */
    pt->step[i] = (pt->end[i] - pt->begin[i]) / ((double)nvals[i] - 1);
    pt->invstep[i] = 1.0 / pt->step[i];

  }
#if defined EAM || defined ADP
#ifndef TBEAM			/* EAM ADP MEAM */
  den_count = ntypes;
  emb_count = ntypes;
#else /* TBEAM */
  if (ntypes == 1) {
    den_count = ntypes + 1;
  } else {
    den_count = ntypes * (ntypes + 1) / 2;
  }
  emb_count = 2 * ntypes;
#endif /* END EAM or TBEAM */
  /* read EAM transfer function rho(r) */
  for (i = paircol; i < paircol + den_count; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.0;
    }
    val += 2;
    ord += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (j = 0; j < nvals[i]; j++) {
      if (2 > fscanf(infile, "%lf %lf\n", ord, val))
	error(1, "Premature end of potential file %s", filename);
      else {
	ord++;
	val++;
      }
      if ((j > 0) && (*(ord - 1) <= *(ord - 2)))
	error(1, "Abscissa not monotonous in potential %d.", i);
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;
    }
    pt->begin[i] = pt->xcoord[pt->first[i]];
    pt->end[i] = pt->xcoord[pt->last[i]];
    /* pt->step is average step length.. */
    pt->step[i] = (pt->end[i] - pt->begin[i]) / ((double)nvals[i] - 1);
    pt->invstep[i] = 1.0 / pt->step[i];
  }

  /* read EAM embedding function F(n) */
  for (i = paircol + den_count; i < paircol + den_count + emb_count; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 1.e30;
    }
    val += 2;
    ord += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(infile, "%lf %lf\n", ord, val))
	error(1, "Premature end of potential file %s", filename);
      else {
	ord++;
	val++;
      }
      if ((j > 0) && (*(ord - 1) <= *(ord - 2)))
	error(1, "Abscissa not monotonous in potential %d.", i);
      if (!invar_pot[i])
	pt->idx[k++] = l++;
      else
	l++;
    }
    pt->begin[i] = pt->xcoord[pt->first[i]];
    pt->end[i] = pt->xcoord[pt->last[i]];
    /* pt->step is average step length.. */
    pt->step[i] = (pt->end[i] - pt->begin[i]) / ((double)nvals[i] - 1);
    pt->invstep[i] = 1.0 / pt->step[i];
  }
#endif /* EAM || ADP */

#ifdef ADP
  /* read ADP dipole function u(r) */
  for (i = paircol + 2 * ntypes; i < 2 * (paircol + ntypes); i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.0;
    }
    val += 2;
    ord += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (j = 0; j < nvals[i]; j++) {
      if (2 > fscanf(infile, "%lf %lf\n", ord, val))
	error(1, "Premature end of potential file %s", filename);
      else {
	ord++;
	val++;
      }
      if ((j > 0) && (*(ord - 1) <= *(ord - 2)))
	error(1, "Abscissa not monotonous in potential %d.", i);
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;
    }
    pt->begin[i] = pt->xcoord[pt->first[i]];
    pt->end[i] = pt->xcoord[pt->last[i]];
    /* pt->step is average step length.. */
    pt->step[i] = (pt->end[i] - pt->begin[i]) / ((double)nvals[i] - 1);
    pt->invstep[i] = 1.0 / pt->step[i];
  }

  /* read adp quadrupole function w(r) */
  for (i = 2 * (paircol + ntypes); i < 3 * paircol + 2 * ntypes; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error(1, "Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 1.e30;
    }
    val += 2;
    ord += 2;
    if ((!invar_pot[i]) && (gradient[i] >> 1))
      pt->idx[k++] = l++;
    else
      l++;
    if ((!invar_pot[i]) && (gradient[i] % 2))
      pt->idx[k++] = l++;
    else
      l++;
    /* read values */
    for (j = 0; j < nvals[i]; j++) {
      if (1 > fscanf(infile, "%lf %lf\n", ord, val))
	error(1, "Premature end of potential file %s", filename);
      else {
	ord++;
	val++;
      }
      if ((j > 0) && (*(ord - 1) <= *(ord - 2)))
	error(1, "Abscissa not monotonous in potential %d.", i);
      if (!invar_pot[i])
	pt->idx[k++] = l++;
      else
	l++;
    }
    pt->begin[i] = pt->xcoord[pt->first[i]];
    pt->end[i] = pt->xcoord[pt->last[i]];
    /* pt->step is average step length.. */
    pt->step[i] = (pt->end[i] - pt->begin[i]) / ((double)nvals[i] - 1);
    pt->invstep[i] = 1.0 / pt->step[i];
  }
#endif /* ADP */

  pt->idxlen = k;
  init_calc_table(pt, &calc_pot);

  return;
}

#endif /* APOT */

#else /* !KIM */

/****************************************************************
 *
 *  read KIM potential with lower and upper limits for each parameter. 
 *
 *  parameters:
 *  	pot_table_t * ... pointer to the potential table
 *  	apot_table_t * ... pointer to the analytic potential table
 *  	char * ... name of the potential file (for error messages)
 *  	FILE * ... open file handle of the potential file
 *
 ****************************************************************/
#ifndef NOLIMITS

void read_pot_table5_no_nolimits(pot_table_t *pt, apot_table_t *apt, char *filename, FILE *infile)
{
	int   i, j, k, l, ret_val;
	char  buffer[255], name[255];
	char *token;
	double *val, *list, temp;
	fpos_t filepos, startpos;
	double  tmp_pmin, tmp_pmax;
	char tmp_value[255];
	int jj, kk, tmp_size;
	void* pkim;
	int status;
	FreeParamType FreeParamSet;


	/* save starting position */
	fgetpos(infile, &startpos);

#ifdef PAIR
	/* read cp */
	if (enable_cp) {

		/* search for cp */
		do {
			fgetpos(infile, &filepos);
			fscanf(infile, "%s", buffer);
		} while (strncmp(buffer, "cp", 2) != 0 && !feof(infile));
		/* and save the position */
		fsetpos(infile, &filepos);

		/* shortcut for apt->number */
		i = apt->number;

		/* allocate memory for global parameters */
		apt->names = (char **)realloc(apt->names, (i + 1) * sizeof(char *));
		apt->names[i] = (char *)malloc(20 * sizeof(char));
		strcpy(apt->names[i], "chemical potentials");
		reg_for_free(apt->names[i], "apt->names[%d] (chem_pot)", i);

		apt->invar_par = (int **)realloc(apt->invar_par, (i + 1) * sizeof(int *));
		apt->invar_par[i] = (int *)malloc((ntypes + 1) * sizeof(int));
		reg_for_free(apt->invar_par[i], "apt->invar_par[%d]", i);

		apt->param_name = (char ***)realloc(apt->param_name, (i + 1) * sizeof(char **));
		apt->param_name[i] = (char **)malloc(ntypes * sizeof(char *));
		reg_for_free(apt->param_name[i], "apt->param_name[%d]", i);

		/* check if the allocation was successfull */
		if (apt->names[i] == NULL || apt->invar_par[i] == NULL || apt->param_name[i] == NULL)
			error(1, "Cannot allocate memory for chemical potentials.");

		/* loop over all atom types */
		for (j = 0; j < ntypes; j++) {
			apt->param_name[i][j] = (char *)malloc(30 * sizeof(char));
			reg_for_free(apt->param_name[i][j], "apt->param_name[%d][%d]", i, j);
			if (apt->param_name[i][j] == NULL)
				error(1, "Cannot allocate memory for chemical potential names.");

			/* read one line */
			if (4 > fscanf(infile, "%s %lf %lf %lf", buffer, &apt->chempot[j], 
						&apt->pmin[i][j], &apt->pmax[i][j]))
				error(1, "Could not read chemical potential for %d. atomtype.", j);

			/* split cp and _# */
			token = strchr(buffer, '_');
			if (token != NULL) {
				strncpy(name, buffer, strlen(buffer) - strlen(token));
				name[strlen(buffer) - strlen(token)] = '\0';
			}
			if (strcmp("cp", name) != 0) {
				fprintf(stderr, "Found \"%s\" instead of \"cp\"\n", name);
				error(1, "No chemical potentials found in %s.\n", filename);
			}

			/* check for invariance and proper value (respect boundaries) */
			apt->invar_par[i][j] = 0.0;
			/* parameter will not be optimized if min==max */
			if (apt->pmin[i][j] == apt->pmax[i][j]) {
				apt->invar_par[i][j] = 1;
				apt->invar_par[i][ntypes]++;
				/* swap min and max if max<min */
			} else if (apt->pmin[i][j] > apt->pmax[i][j]) {
				temp = apt->pmin[i][j];
				apt->pmin[i][j] = apt->pmax[i][j];
				apt->pmax[i][j] = temp;
				/* reset value if >max or <min */
			} else if ((apt->values[i][j] < apt->pmin[i][j])
					|| (apt->values[i][j] > apt->pmax[i][j])) {
				/* Only print warning if we are optimizing */
				if (opt) {
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
		printf(" - Enabled %d chemical potential(s)\n", ntypes);

		/* disable composition nodes for now */
#ifdef CN
		/* read composition nodes */
		if (2 > fscanf(infile, "%s %d", buffer, &compnodes)) {
			if (strcmp("type", buffer) == 0)
				compnodes = -1;
			else
				error(1, "Could not read number of composition nodes from potential file.\n");
		}
		if (strcmp(buffer, "cn") != 0 && ntypes > 1 && compnodes != -1)
			error(1, "No composition nodes found in %s.\nUse \"cn 0\" for none.\n", filename);
		if (ntypes == 1) {
			compnodes = 0;
		}
		if (compnodes != -1) {
			apt->values[apt->number] =
				(double *)realloc(apt->values[apt->number], (ntypes + compnodes) * sizeof(double));
			apt->pmin[apt->number] =
				(double *)realloc(apt->pmin[apt->number], (ntypes + compnodes) * sizeof(double));
			apt->pmax[apt->number] =
				(double *)realloc(apt->pmax[apt->number], (ntypes + compnodes) * sizeof(double));
			apt->chempot = apt->values[apt->number];
			compnodelist = (double *)malloc((ntypes + compnodes) * sizeof(double));

			for (j = 0; j < compnodes; j++) {
				if (4 > fscanf(infile, "%lf %lf %lf %lf", &compnodelist[j],
							&apt->chempot[ntypes + j], &apt->pmin[apt->number][ntypes + j],
							&apt->pmax[apt->number][ntypes + j]))
					error(1, "Could not read composition node %d\n", j + 1);
				if (apt->pmin[apt->number][ntypes + j] > apt->chempot[ntypes + j]
						|| apt->pmax[apt->number][ntypes + j] < apt->chempot[ntypes + j])
					error(1, "composition node %d is out of bounds.\n", j + 1);
			}

			/* check compnodes for valid values */
			if (ntypes == 2) {
				for (j = 0; j < compnodes; j++)
					if (compnodelist[j] > 1 || compnodelist[j] < 0)
						error(1, "Composition node %d is %f but should be inside [0,1].\n", j + 1, compnodelist[j]);
			}
		}
		if (compnodes != -1)
			printf("Enabled chemical potentials with %d extra composition node(s).\n", compnodes);
		if (compnodes == -1)
			compnodes = 0;
#endif /* CN */
	}
#endif /* PAIR */


	/* read the keywords and the names of parameters that will be optimized */
	read_potential_keyword(pt, filename, infile, &FreeParamSet);

	/* write the descriptor.kim file for the test */
	write_temporary_descriptor_file(kim_model_name);

	/* create KIM object with 1 atom and 1 species */
	setup_KIM_API_object(&pkim, 1, 1, kim_model_name);
	
  get_compute_const(pkim);

  /* initialze the data struct for the free parameters with type double */
  get_free_param_double(pkim, &FreeParamSet);

  /* nest the optimizable params */
  nest_optimizable_param(pkim, &FreeParamSet, name_opt_param, num_opt_param);

	/* There is only 1 potential for KIM Model. */
	i = 0;
	apt->n_par[i] = FreeParamSet.Nnestedvalue; 
	apt->total_par += apt->n_par[i];
	
	/* copy name */
	strcpy(apt->names[i], kim_model_name);

	/* set begin */
	apt->begin[i] = 0.0;

	/* allocate memory */
	apt->values[i] = (double *)malloc(apt->n_par[i] * sizeof(double));
	reg_for_free(apt->values[i], "apt->values[%d]", i);
	/* last slot stores the total number of invariable parameters */
	apt->invar_par[i] = (int *)malloc((apt->n_par[i] + 1) * sizeof(int));
	reg_for_free(apt->invar_par[i], "apt->invar_par[%d]", i);
	apt->pmin[i] = (double *)malloc(apt->n_par[i] * sizeof(double));
	reg_for_free(apt->pmin[i], "apt->pmin[%d]", i);
	apt->pmax[i] = (double *)malloc(apt->n_par[i] * sizeof(double));
	reg_for_free(apt->pmax[i], "apt->pmax[%d]", i);
	if(NULL == apt->pmin[i] || NULL == apt->pmax[i]
	   || NULL == apt->values[i] || NULL == apt->invar_par[i])
		error(1, "Cannot allocate memory for potential paramters.\nAborting");


	/* read parameter values */
  /* j: jth optimizable parameter 
     k: kth free parameter, index of jth optimizable param in the free param struct 
    jj: the current slot in the nested apot->value
    kk: the current slot of the jth param
  */
	jj = 0;  
	fsetpos(infile, &startpos);
	for (j = 0; j < num_opt_param; j++) {
		buffer[0] = '\0';
		name[0] = '\0';
		/* find the keyword `PARAM_FREE_*' */
		do {
			fgets(buffer, 255, infile);
			sscanf(buffer, "%s", name);
		} while (strncmp(name, "PARAM_FREE", 10) != 0 && !feof(infile));
		/* k is the palce(slot) that the param in the FreeParam struct */
		for ( k = 0; k < FreeParamSet.Nparam; k++) {
			if (strcmp(FreeParamSet.name[k], name_opt_param[j]) == 0)
				break;
		}
		/* read in the value and pmin and pmax */
		for (kk = 0; kk < size_opt_param[j]; kk++) {
			fgets(buffer, 255, infile);
			ret_val = sscanf(buffer, "%s %lf %lf", tmp_value, &tmp_pmin, &tmp_pmax);
			if (ret_val == 3) {
				apt->pmin[i][jj] = tmp_pmin;
				apt->pmax[i][jj] = tmp_pmax;
				/*if give `KIM', use KIM parameter, otherwise use readin */
				if (strncmp(tmp_value, "KIM", 3) == 0) {
					apt->values[i][jj] = *FreeParamSet.nestedvalue[jj];
				} else if (1 > sscanf(tmp_value, "%lf", &apt->values[i][jj])) {
					error(1, "First data for parameter '%s' corrupted, it should be a type float or 'KIM'.\n",
                    FreeParamSet.name[k]);
				}
				if(apt->pmin[i][j]>apt->values[i][j] || apt->pmax[i][j] < apt->values[i][j])
				{
					error(1, "Value of parameter '%s' is not within its limits in file '%s'.\n",
								name_opt_param[j], filename);
				}
			} else {
				error(0, "Data for parameter '%s' corrupted in file '%s'.\n",
                 FreeParamSet.name[k], filename);
        error(1, "Or you may want to try compiling with nolimits enabled.\n", kk+1,FreeParamSet.name[k]);
			}
			jj++;
		}  
	}


	/* invarable parameters. Here all are optimizable. */
	/* the last slot is the number of invar */
	for (j = 0; j < apt->n_par[i]; j++) {
		apt->invar_par[i][j] = 0;
	}
	apt->invar_par[i][apt->n_par[i]] = 0;


	/* find cutoff */
	int have_cutoff;
  int have_PARAM_FREE_cutoff;
          have_cutoff = 0;	
  /* If it is anaytic potential (every PARAM_FREE* is scalar), keyword `cutoff' has
	 * to be specified in the input file. You can give value directly, or give `KIM' 
   * to use the cutoff from KIM model.
   * If it is tabualted potential, cutoff will read in from KIM model, even
	 * if it is given in the input file. */
	
  /* whether we have PARAM_FREE_cutoff in the KIM Model. 
  If not, cannot publish cutoff, then have_cutoff needs to be 0. */

  KIM_API_get_data(pkim, "PARAM_FREE_cutoff", &status);
   if (KIM_STATUS_OK == status)  have_PARAM_FREE_cutoff = 1;


	fsetpos(infile, &startpos);
	do {
		fgets(buffer, 255, infile);
		sscanf(buffer, "%s", name);
	} while (strcmp(name, "cutoff") != 0 && !feof(infile));

	if (apt->n_par[i] == num_opt_param) {
		if (strcmp(name, "cutoff") == 0) {
			if(2 > sscanf(buffer, "%s %s", name, tmp_value)) {
				error(1,"Error reading in cutoff in file '%s'.\n", filename);
			} else {
				if (strcmp(tmp_value, "KIM") == 0) {
					have_cutoff = 0;
				} else if (have_PARAM_FREE_cutoff){
					if(1 > sscanf(tmp_value,"%lf", &apt->end[i])){
				    error(1,"Error reading in cutoff in file '%s'.\n", filename);
					}
					have_cutoff = 1; 
				} else {
          have_cutoff = 0;
        }
			}
		} else {
			error(1,"'cutoff' is missing in file: %s.\n", filename);
		}
	} else {
		if (strcmp(name, "cutoff") == 0) {
			printf("'cutoff' read in from file '%s' is deprecated. Will use the 'cutoff' in "
					" the KIM model throughout the fitting.\n", filename);
		}
		have_cutoff = 0;
	}

	/* get cutoff from KIM if cutoff has not been determined yet. */
	if (!have_cutoff) {
    int status;
    double* pcutoff;
    
    pcutoff = KIM_API_get_data(pkim, "cutoff", &status);
    if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", status);
      return(status);
    }
      apt->end[i] = *pcutoff;
  }


#ifdef PAIR
	if (enable_cp) {
		cp_start = apt->total_par - apt->globals + ntypes * (ntypes + 1);
		apt->total_par += (ntypes + compnodes - apt->invar_par[apt->number][ntypes]);
	}
#endif /* PAIR */


	/* no derivative slots are needed */  

	/* initialize function table */
	for (i = 0; i < apt->number; i++) {
		pt->begin[i] = apt->begin[i];
		pt->end[i] = apt->end[i];
		pt->step[i] = 0;
		pt->invstep[i] = 0;
		if (i == 0)
			pt->first[i] = 0;
		else
			pt->first[i] = pt->last[i - 1] + 1;
		pt->last[i] = pt->first[i] + apt->n_par[i] - 1;
	}
	pt->len = pt->first[apt->number - 1] + apt->n_par[apt->number - 1];
	if (have_globals)
		pt->len += apt->globals;

#ifdef PAIR
	if (enable_cp) {
		pt->len += (ntypes + compnodes);
	}
#endif /* PAIR */

	pt->table = (double *)malloc(pt->len * sizeof(double));
	reg_for_free(pt->table, "pt->table");

	calc_list = (double *)malloc(pt->len * sizeof(double));
	reg_for_free(calc_list, "calc_list");

	pt->idx = (int *)malloc(pt->len * sizeof(int));
	reg_for_free(pt->idx, "pt->idx");

	apt->idxpot = (int *)malloc(apt->total_par * sizeof(int));
	reg_for_free(apt->idxpot, "apt->idxpot");

	apt->idxparam = (int *)malloc(apt->total_par * sizeof(int));
	reg_for_free(apt->idxparam, "apt->idxparam");

	if ((NULL == pt->table) || (NULL == pt->idx) || (apt->idxpot == NULL)
			|| (apt->idxparam == NULL))
		error(1, "Cannot allocate memory for potential table.\n");
	for (i = 0; i < pt->len; i++) {
		pt->table[i] = 0;
		calc_list[i] = 0;
		pt->idx[i] = 0;
	}


	/* no derivative slots are needed, val += 2 etc. is disabled )  */  
	/* this is the indirect index */
	k = 0;
	l = 0;
	val = pt->table;
	list = calc_list;
	for (i = 0; i < apt->number; i++) {	/* loop over potentials */
		for (j = 0; j < apt->n_par[i]; j++) {	/* loop over parameters */
			*val = apt->values[i][j];
			*list = apt->values[i][j];
			val++;
			list++;
			if (!invar_pot[i] && !apt->invar_par[i][j]) {
				pt->idx[k] = l++;      /*index, the optimiable parameters in pt->table. for
																 example, idx[0] = 2 indicates that the the first
																 variable parameters lies in slot 2 of pt->table */
				apt->idxpot[k] = i;   /* index, the optimizable parameters come from which
																 potential? e.g. idxpot[0] = 0, indicates that the first 
																 variable parameter is from the first potential*/
				apt->idxparam[k++] = j;/*index, the optimizable parameters is the ?th
																 parammeter of the potential. Sould be used together
																 with idxpot. e.g. idxparam[0] = 1 indicates that the
																 first variable parameter is the 2nd parameter 
																 of potential idxpot[0] */
			} else
				l++;
		}
		if (!invar_pot[i])
			pt->idxlen += apt->n_par[i] - apt->invar_par[i][apt->n_par[i]];
		apt->total_par -= apt->invar_par[i][apt->n_par[i]];
	}

#ifdef PAIR
	if (enable_cp) {
		init_chemical_potential(ntypes);
		i = apt->number;
		for (j = 0; j < (ntypes + compnodes); j++) {
			*val = apt->values[i][j];
			val++;
			if (!apt->invar_par[i][j]) {
				pt->idx[k] = l++;
				apt->idxpot[k] = i;
				apt->idxparam[k++] = j;
			}
		}
		pt->idxlen += (ntypes + compnodes - apt->invar_par[apt->number][ntypes]);
		global_idx += (ntypes + compnodes - apt->invar_par[apt->number][ntypes]);
	}
#endif /* PAIR */

#ifdef NOPUNISH
	if (opt)
		warning("Gauge degrees of freedom are NOT fixed!\n");
#endif /* NOPUNISH */

	/* free memory */
	free_model_object(&pkim);

	printf(" - Successfully read potential parameters that will be optimized.\n");
	return;
}

#else /* NOLIMITS */

/****************************************************************
 *
 *  read KIM potential without specifying the lower and upper limits
 *  for each parameter. 
 *
 *  parameters:
 * 	  pot_table_t * ... pointer to the potential table
 *  	int ... number of potential functions
 *  	char * ... name of the potential file (for error messages)
 *  	FILE * ... open file handle of the potential file
 *
 ****************************************************************/

void read_pot_table5_with_nolimits(pot_table_t *pt, int size, char *filename, FILE *infile)
{
	int   i, j, k, ret_val;
	char  buffer[255], name[255];
	fpos_t filepos, startpos; 
	void* pkim;
  int status;
	FreeParamType FreeParamSet;

	/* save starting position */
	fgetpos(infile, &startpos);

	/* read the keywords and the names of parameters that will be optimized */
	read_potential_keyword(pt, filename, infile, &FreeParamSet);

	/* write the descriptor.kim file for the test */
	write_temporary_descriptor_file(kim_model_name);

	/* create KIM object with 1 atom and 1 species */
	setup_KIM_API_object(&pkim, 1, 1, kim_model_name);
  
  get_compute_const(pkim);

  /* initialze the data struct for the free parameters with type double */
  get_free_param_double(pkim, &FreeParamSet);

  /* nest the optimizable params */
  nest_optimizable_param(pkim, &FreeParamSet, name_opt_param, num_opt_param);

	/* some potential table value */
	pt->first[0] = 0;
	pt->last[0] = pt->first[0] + FreeParamSet.Nnestedvalue - 1;
	pt->len = FreeParamSet.Nnestedvalue;
	pt->idxlen = FreeParamSet.Nnestedvalue;

	/* allocate the function table */
	pt->table = (double *)malloc(pt->len * sizeof(double));
	reg_for_free(pt->table, "pt->table");
	pt->idx = (int *)malloc(pt->len * sizeof(int));
	reg_for_free(pt->idx, "pt->idx");
	if ((NULL == pt->table) || (NULL == pt->idx) )
		error(1, "Cannot allocate memory for potential table");

	/* copy parameters values to potfit */
	for (i = 0; i < FreeParamSet.Nnestedvalue; i++) {
		pt->table[i] = *FreeParamSet.nestedvalue[i];
		pt->idx[i] = i;
	}

	/* set end (end is cutoff) */
  double* pcutoff;
  
  pcutoff = KIM_API_get_data(pkim, "cutoff", &status);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", status);
    return(status);
  }
    pt->end[0] = *pcutoff;


	/* free memory */
	free_model_object(&pkim);
  
	printf(" - Successfully read potential parameters that will be optimized.\n");
	return;
}

#endif /* NOLIMITS */
#endif /* !KIM */
/* added ends */


/****************************************************************
 *
 *  init_calc_table: Initialize table used for calculation.
 *
 *  *  Header:  one line for each function with
 *           rbegin rstart npoints
 *
 *  Table: Center, width, amplitude of Gaussians,
 *         functions separated by blank lines
 *
 ****************************************************************/

void init_calc_table(pot_table_t *optt, pot_table_t *calct)
{
#ifdef APOT
  int   i, j, index, x = 0, size;
  double *val, f, h;
#endif /* APOT */

  switch (format) {
#ifdef APOT
      case 0:
	{
	  /* allocate memory for calc_pot potential table */
	  size = apot_table.number;
	  calct->len = size * APOT_STEPS + 2 * optt->ncols + ntypes + compnodes;
	  calct->idxlen = APOT_STEPS;
	  calct->ncols = optt->ncols;
	  calct->begin = optt->begin;
	  calct->end = optt->end;
	  calct->first = (int *)malloc(size * sizeof(int));
	  reg_for_free(calct->first, "calct->first");
	  calct->last = (int *)malloc(size * sizeof(int));
	  reg_for_free(calct->last, "calct->last");
	  calct->step = (double *)malloc(size * sizeof(double));
	  reg_for_free(calct->step, "calct->step");
	  calct->invstep = (double *)malloc(size * sizeof(double));
	  reg_for_free(calct->invstep, "calct->invstep");
	  calct->xcoord = (double *)malloc(calct->len * sizeof(double));
	  reg_for_free(calct->xcoord, "calct->xcoord");
	  calct->table = (double *)malloc(calct->len * sizeof(double));
	  reg_for_free(calct->table, "calct->table");
	  calct->d2tab = (double *)malloc(calct->len * sizeof(double));
	  reg_for_free(calct->d2tab, "calct->d2tab");
	  calct->idx = (int *)malloc(calct->len * sizeof(int));
	  reg_for_free(calct->idx, "calct->idx");
	  if (calct->first == NULL || calct->last == NULL || calct->step == NULL
	    || calct->invstep == NULL || calct->xcoord == NULL
	    || calct->table == NULL || calct->d2tab == NULL || calct->idx == NULL)
	    error(1, "Cannot allocate info block for calc potential table\n");

	  for (i = 0; i < calct->len; i++) {
	    calct->xcoord[i] = 0.0;
	    calct->table[i] = 0.0;
	    calct->d2tab[i] = 0.0;
	  }

	  /* initialize the calc_pot table */
	  for (i = 0; i < size; i++) {
	    val = apot_table.values[i];
	    h = apot_table.values[i][apot_table.n_par[i] - 1];
	    calct->table[i * APOT_STEPS + i * 2] = 10e30;
	    calct->table[i * APOT_STEPS + i * 2 + 1] = 0;
	    calct->first[i] = (x += 2);
	    calct->last[i] = (x += APOT_STEPS - 1);
	    x++;
	    calct->step[i] = (calct->end[i] - calct->begin[i]) / (APOT_STEPS - 1);
	    calct->invstep[i] = 1.0 / calct->step[i];
	    for (j = 0; j < APOT_STEPS; j++) {
	      index = i * APOT_STEPS + (i + 1) * 2 + j;
	      calct->xcoord[index] = calct->begin[i] + j * calct->step[i];

	      apot_table.fvalue[i] (calct->xcoord[index], val, &f);
	      calct->table[index] = smooth_pot[i] ? f * cutoff(calct->xcoord[index], calct->begin[i], h) : f;
	      calct->idx[i * APOT_STEPS + j] = index;
	    }
	  }
	}
	break;
#else
      case 3:			/* fall through */
      case 4:
	calct->len = optt->len;
	calct->idxlen = optt->idxlen;
	calct->ncols = optt->ncols;
	calct->begin = optt->begin;
	calct->end = optt->end;
	calct->step = optt->step;
	calct->invstep = optt->invstep;
	calct->first = optt->first;
	calct->last = optt->last;
	calct->xcoord = optt->xcoord;
	calct->table = optt->table;
	calct->d2tab = optt->d2tab;
	calct->idx = optt->idx;
#endif /* APOT */
  }

  return;
}

#ifdef APOT

/****************************************************************
 *
 * update apot_table from opt_pot_table, including globals
 *
 ****************************************************************/

void update_apot_table(double *xi)
{
  int   i = 0, j = 0, m = 0, n = 0;

  for (i = 0; i < ndim; i++)
    apot_table.values[apot_table.idxpot[i]][apot_table.idxparam[i]] = xi[idx[i]];
  if (have_globals) {
    for (i = 0; i < apot_table.globals; i++) {
      for (j = 0; j < apot_table.n_glob[i]; j++) {
	m = apot_table.global_idx[i][j][0];
	n = apot_table.global_idx[i][j][1];
	apot_table.values[m][n] = *(xi + global_idx + i);
      }
    }
  }

  return;
}

/****************************************************************
 *
 * update calc_pot.table from opt_pot.table, including globals
 *
 ****************************************************************/

void update_calc_table(double *xi_opt, double *xi_calc, int do_all)
{
  int   i, j, k, m, n, change;
  double f, h = 0;
  double *list, *val;

  val = xi_opt;
  list = calc_list + 2;
  /* copy global parameters to the right positions */
  if (have_globals) {
    for (i = 0; i < apot_table.globals; i++) {
      for (j = 0; j < apot_table.n_glob[i]; j++) {
	m = apot_table.global_idx[i][j][0];
	n = apot_table.global_idx[i][j][1];
	*(val + opt_pot.first[m] + n) = *(val + global_idx + i);
      }
    }
  }
  for (i = 0; i < calc_pot.ncols; i++) {
    if (smooth_pot[i] && (do_all || !invar_pot[i])) {
      h = *(val + 1 + apot_table.n_par[i]);
      if (h == 0)
	error(1, "The cutoff parameter for potential %d is 0!", i);
    }

    (*val) = apot_grad(calc_pot.begin[i], val + 2, apot_table.fvalue[i]);
    val += 2;
    /* check if something has changed */
    change = 0;
    for (j = 0; j < apot_table.n_par[i]; j++) {
      if (list[j] != val[j]) {
	change = 1;
	list[j] = val[j];
      }
    }
    if (do_all || (change && !invar_pot[i])) {
      for (j = 0; j < APOT_STEPS; j++) {
	k = i * APOT_STEPS + (i + 1) * 2 + j;
	apot_table.fvalue[i] (calc_pot.xcoord[k], val, &f);
	*(xi_calc + k) = smooth_pot[i] ? f * cutoff(calc_pot.xcoord[k], apot_table.end[i], h) : f;
	if (isnan(f) || isnan(*(xi_calc + k))) {
#ifdef DEBUG
	  error(0, "Potential value was nan or inf. Aborting.\n");
	  error(0, "This occured in potential %d (%s)\n", i, apot_table.names[i]);
	  error(0, "at distance r=%f with the parameters:\n", calc_pot.xcoord[k]);
	  for (m = 0; m < apot_table.n_par[i]; m++)
	    error(0, "%s %f\n", apot_table.param_name[i][m], *(val + m));
	  if (smooth_pot[i])
	    error(0, "h %f\n", h);
#endif /* DEBUG */
	  error(1, "Potential value was nan or inf. Aborting.\n");
	}
      }
    }
    val += apot_table.n_par[i];
    list += apot_table.n_par[i] + 2;
  }

  return;
}

#endif /* APOT */
