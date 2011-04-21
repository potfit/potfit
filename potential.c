/****************************************************************
 *
 * potential.c: Routines for reading, writing and interpolating a
 *	potential table
 *
 ****************************************************************
 *
 * Copyright 2002-2011 Peter Brommer, Franz G"ahler, Daniel Schopf,
 *                     Philipp Beck
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://www.itap.physik.uni-stuttgart.de/
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
#include "splines.h"
#include "utils.h"

#define NPLOT 1000

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
  int   size, i, j, k = 0, *nvals, ncols, npots = 0;
#ifdef APOT
  apot_table_t *apt = &apot_table;
#else
  real *val;
#endif /* APOT */

  /* open file */
  infile = fopen(filename, "r");
  if (NULL == infile)
    error("Could not open file %s\n", filename);

  /* read the header */
  do {
    /* read one line */
    res = fgets(buffer, 1024, infile);
    if (NULL == res)
      error("Unexpected end of file in %s", filename);
    /* check if it is a header line */
    if (buffer[0] != '#')
      error("Header corrupt in file %s", filename);
    /* stop after last header line */
    if (buffer[1] == 'E') {
      end_header = 1;
    }
    if (buffer[1] == 'T') {
      if ((str = strchr(buffer + 3, '\n')) != NULL)
	*str = '\0';
      if (strcmp(buffer + 3, interaction_name) != 0)
	error
	  ("Wrong potential type!\nThis binary only supports %s-potentials.\nYour potential file contains a %s-potential.\n",
	  interaction_name, buffer + 3);
    }
    /* invariant potentials */
    else if (buffer[1] == 'I') {
      if (have_format) {
#ifdef APOT
	apot_table.invar_pots = 0;
#endif /* APOT */
	/* gradient complete */
	for (i = 0; i < size; i++) {
	  str = strtok(((i == 0) ? buffer + 2 : NULL), " \t\r\n");
	  if (str == NULL) {
	    error("Not enough items in #I header line.");
	  } else {
	    ((int *)invar_pot)[i] = atoi(str);
#ifdef APOT
	    apot_table.invar_pots++;
#endif /* APOT */
	  }
	}
	have_invar = 1;
      } else
	error("#I needs to be specified after #F in file %s", filename);
    }
#ifndef APOT
    else if (buffer[1] == 'G') {
      if (have_format) {
	/* gradient complete */
	for (i = 0; i < size; i++) {
	  str = strtok(((i == 0) ? buffer + 2 : NULL), " \t\r\n");
	  if (str == NULL)
	    error("Not enough items in #G header line.");
	  else
	    ((int *)gradient)[i] = atoi(str);
	}
	have_grad = 1;
      } else
	error("#G needs to be specified after #F in file %s", filename);
    }
#endif /* !APOT */

    /* see if it is the format line */
    else if (buffer[1] == 'F') {
      /* format complete? */
      if (2 != sscanf((const char *)(buffer + 2), "%d %d", &format, &size))
	error("Corrupt format header line in file %s", filename);
#ifndef APOT
      if (format == 0)
	error("potfit binary compiled without analytic potential support.\n");
#else
      if (format > 0)
	error("potfit binary compiled without tabulated potential support.\n");
#endif /* !APOT */

      ncols = ntypes * (ntypes + 1) / 2;
      /* right number of columns? */
      switch (interaction) {
	  case I_EAM:
	    npots = ncols + 2 * ntypes;
	    break;
	  case I_ADP:
	    npots = 3 * ncols + 2 * ntypes;
	    break;
	  default:
	    npots = ncols;
      }

      if (size == npots) {
	printf("Using %s potentials from file \"%s\".\n", interaction_name,
	  filename);
      } else
	error
	  ("Wrong number of data columns in file \"%s\",\n should be %d for %s, but are %d.",
	  filename, npots, interaction_name, size);
      /* recognized format? */
      if ((format != 0) && (format != 3) && (format != 4))
	error("Unrecognized format specified for file %s", filename);
      gradient = (int *)malloc(size * sizeof(int));
      invar_pot = (int *)malloc(size * sizeof(int));
#ifdef APOT
      smooth_pot = (int *)malloc(size * sizeof(int));
#endif
      for (i = 0; i < size; i++) {
	gradient[i] = 0;
	invar_pot[i] = 0;
#ifdef APOT
	smooth_pot[i] = 0;
#endif
      }
      have_format = 1;
    }
  } while (!end_header);

  /* do we have a format in the header? */
  if (!have_format)
    error("Format not specified in header of file %s", filename);
  else if (format != 0)
    printf("Potential file format %d detected.\n", format);
  else
    printf("Potential file format %d (analytic potentials) detected.\n",
      format);

  /* allocate info block of function table */
  pt->len = 0;
  pt->ncols = size;
  pt->begin = (real *)malloc(size * sizeof(real));
  pt->end = (real *)malloc(size * sizeof(real));
  pt->step = (real *)malloc(size * sizeof(real));
  pt->invstep = (real *)malloc(size * sizeof(real));
  pt->first = (int *)malloc(size * sizeof(int));
  pt->last = (int *)malloc(size * sizeof(int));
  nvals = (int *)malloc(size * sizeof(int));
  if ((pt->begin == NULL) || (pt->end == NULL) || (pt->step == NULL)
    || (pt->invstep == NULL) || (pt->first == NULL) || (pt->last == NULL)
    || (nvals == NULL))
    error("Cannot allocate info block for potential table %s", filename);
#ifdef APOT
  /* allocate memory for analytic potential table */
  apt->number = size;
  apt->total_par = 0;

  apt->n_par = (int *)malloc(size * sizeof(int));
  apt->begin = (real *)malloc(size * sizeof(real));
  apt->end = (real *)malloc(size * sizeof(real));
  apt->param_name = (char ***)malloc(size * sizeof(char **));
  apt->fvalue = (fvalue_pointer *) malloc(size * sizeof(fvalue_pointer));
#ifdef PAIR
  if (enable_cp) {
    apt->values = (real **)malloc((size + 1) * sizeof(real *));
    apt->values[size] = (real *)malloc(ntypes * sizeof(real));
    apt->invar_par = (int **)malloc(size * sizeof(int *));
    apt->chempot = apt->values[size];
    apt->pmin = (real **)malloc((size + 1) * sizeof(real *));
    apt->pmin[size] = (real *)malloc(ntypes * sizeof(real));
    apt->pmax = (real **)malloc((size + 1) * sizeof(real *));
    apt->pmax[size] = (real *)malloc(ntypes * sizeof(real));
  } else {
#elif defined COULOMB
  if (1) {
    apt->ratio = (real *)malloc(ntypes * sizeof(real));
    apt->values = (real **)malloc((size + 4) * sizeof(real *));
    apt->param_name = (char ***)malloc((size + 4) * sizeof(char **));
    apt->pmin = (real **)malloc((size + 4) * sizeof(real *));
    apt->pmax = (real **)malloc((size + 4) * sizeof(real *));
    apt->invar_par = (int **)malloc((size + 4) * sizeof(int *));

    apt->values[size] = (real *)malloc((ntypes - 1) * sizeof(real));
    apt->param_name[size] = (char **)malloc((ntypes - 1) * sizeof(char *));
    apt->pmin[size] = (real *)malloc((ntypes - 1) * sizeof(real));
    apt->pmax[size] = (real *)malloc((ntypes - 1) * sizeof(real));
    apt->invar_par[size] = (int *)malloc((ntypes - 1) * sizeof(int));

    apt->values[size + 1] = (real *)malloc(ntypes * sizeof(real));
    apt->param_name[size + 1] = (char **)malloc(ntypes * sizeof(char *));
    apt->pmin[size + 1] = (real *)malloc(ntypes * sizeof(real));
    apt->pmax[size + 1] = (real *)malloc(ntypes * sizeof(real));
    apt->invar_par[size + 1] = (int *)malloc(ntypes * sizeof(int));

    for (i = 2; i < 4; i++) {
      apt->values[size + i] = (real *)malloc(size * sizeof(real));
      apt->param_name[size + i] = (char **)malloc(size * sizeof(char *));
      apt->pmin[size + i] = (real *)malloc(size * sizeof(real));
      apt->pmax[size + i] = (real *)malloc(size * sizeof(real));
      apt->invar_par[size + i] = (int *)malloc(size * sizeof(int));
    }
    apt->charge = apt->values[size];
#ifdef DIPOLE
    apt->dp_alpha = apt->values[size + 1];
    apt->dp_b = apt->values[size + 2];
    apt->dp_c = apt->values[size + 3];
#endif
  } else {
#endif /* COULOMB */
    apt->values = (real **)malloc(size * sizeof(real *));
    apt->invar_par = (int **)malloc(size * sizeof(int *));
    apt->pmin = (real **)malloc(size * sizeof(real *));
    apt->pmax = (real **)malloc(size * sizeof(real *));
#if defined PAIR || defined COULOMB
  }
#endif
  apt->names = (char **)malloc(size * sizeof(char *));
  for (i = 0; i < size; i++) {
    apt->names[i] = (char *)malloc(20 * sizeof(char));
  }
  if ((apt->n_par == NULL) || (apt->begin == NULL) || (apt->end == NULL)
    || (apt->fvalue == NULL) || (apt->names == NULL) || (apt->pmin == NULL)
    || (apt->pmax == NULL) || (apt->param_name == NULL)
    || (apt->values == NULL))
    error("Cannot allocate info block for analytic potential table %s",
      filename);
#endif /* APOT */
  switch (format) {
#ifdef APOT
      case 0:
	read_pot_table0(pt, apt, filename, infile);
	break;
#else
      case 3:
	read_pot_table3(pt, size, ncols, nvals, filename, infile);
	break;
      case 4:
	read_pot_table4(pt, size, ncols, nvals, filename, infile);
#endif /* APOT */
  }
  fclose(infile);

  /* compute rcut and rmin */
  rcut = (real *)malloc(ntypes * ntypes * sizeof(real));
  if (NULL == rcut)
    error("Cannot allocate rcut");
  rmin = (real *)malloc(ntypes * ntypes * sizeof(real));
  if (NULL == rmin)
    error("Cannot allocate rmin");
  for (i = 0; i < ntypes; i++)
    for (j = 0; j < ntypes; j++) {
      k = (i <= j) ? i * ntypes + j - ((i * (i + 1)) / 2)
	: j * ntypes + i - ((j * (j + 1)) / 2);
      rmin[i * ntypes + j] = pt->begin[k];
      rcut[i * ntypes + j] = pt->end[k];
    }
#if defined EAM || defined ADP
  for (i = 0; i < ntypes; i++) {
    for (j = 0; j < ntypes; j++) {
      rcut[i * ntypes + j] =
	MAX(rcut[i * ntypes + j], pt->end[(ntypes * (ntypes + 1)) / 2 + i]);
      rcut[i * ntypes + j] =
	MAX(rcut[i * ntypes + j], pt->end[(ntypes * (ntypes + 1)) / 2 + j]);
      rmin[i * ntypes + j] =
	MIN(rmin[i * ntypes + j], pt->begin[(ntypes * (ntypes + 1)) / 2 + i]);
      rmin[i * ntypes + j] =
	MIN(rmin[i * ntypes + j], pt->begin[(ntypes * (ntypes + 1)) / 2 + j]);
    }
  }
#endif /* EAM || ADP */

  paircol = (ntypes * (ntypes + 1)) / 2;

#ifndef APOT
  /* read maximal changes file */
  maxchange = (real *)malloc(pt->len * sizeof(real));
  if (usemaxch) {
    /* open file */
    infile = fopen(maxchfile, "r");
    if (NULL == infile)
      error("Could not open file %s\n", maxchfile);
    val = maxchange;
    for (i = 0; i < pt->len; i++) {
      if (1 > fscanf(infile, " %lf\n", val))
	error("Premature end of maxch file %s", maxchfile);
      else
	val++;
    }
    fclose(infile);
  }
#endif /* APOT */
  for (i = 0; i < ntypes; i++) {
    for (j = 0; j < ntypes; j++) {
      rcutmax = MAX(rcutmax, rcut[i + ntypes * j]);
    }
  }

  /* clean up locals and mark globals for later */
  free(nvals);
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
  for (i = 0; i < 4; i++) {
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

#ifdef APOT

/****************************************************************
 *
 *  read potential in analytic format:
 *  	for more information an how to specify an analytic potential
 *  	please check the documentation
 *
 ****************************************************************/

void read_pot_table0(pot_table_t *pt, apot_table_t *apt, char *filename,
  FILE *infile)
{
  int   i, j, k, l, ret_val;
  char  buffer[255], name[255];
  char *token;
  real *val, *list, temp;
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
    fsetpos(infile, &filepos);

    for (i = 0; i < ntypes; i++) {
      if (4 > fscanf(infile, "%s %lf %lf %lf", buffer, &apt->chempot[i],
	  &apt->pmin[apt->number][i], &apt->pmax[apt->number][i]))
	error("Could not read chemical potential for atomtype #%d.", i);

      /* split cp and _# */
      token = strchr(buffer, '_');
      if (token != NULL) {
	strncpy(name, buffer, strlen(buffer) - strlen(token));
	name[strlen(buffer) - strlen(token)] = '\0';
      }
      if (strcmp("cp", name) != 0) {
	fprintf(stderr, "Found \"%s\" instead of \"cp\"\n", name);
	error("No chemical potentials found in %s.\n", filename);
      }
    }
    printf("Enabled chemical potentials.\n");

    /* disable composition nodes for now */
#ifdef CN
    /* read composition nodes */
    if (2 > fscanf(infile, "%s %d", buffer, &compnodes)) {
      if (strcmp("type", buffer) == 0)
	compnodes = -1;
      else
	error
	  ("Could not read number of composition nodes from potential file.\n");
    }
    if (strcmp(buffer, "cn") != 0 && ntypes > 1 && compnodes != -1)
      error("No composition nodes found in %s.\nUse \"cn 0\" for none.\n",
	filename);
    if (ntypes == 1) {
      compnodes = 0;
    }
    if (compnodes != -1) {
      apt->values[apt->number] =
	(real *)realloc(apt->values[apt->number],
	(ntypes + compnodes) * sizeof(real));
      apt->pmin[apt->number] =
	(real *)realloc(apt->pmin[apt->number],
	(ntypes + compnodes) * sizeof(real));
      apt->pmax[apt->number] =
	(real *)realloc(apt->pmax[apt->number],
	(ntypes + compnodes) * sizeof(real));
      apt->chempot = apt->values[apt->number];
      compnodelist = (real *)malloc((ntypes + compnodes) * sizeof(real));

      for (j = 0; j < compnodes; j++) {
	if (4 > fscanf(infile, "%lf %lf %lf %lf", &compnodelist[j],
	    &apt->chempot[ntypes + j], &apt->pmin[apt->number][ntypes + j],
	    &apt->pmax[apt->number][ntypes + j]))
	  error("Could not read composition node %d\n", j + 1);
	if (apt->pmin[apt->number][ntypes + j] > apt->chempot[ntypes + j]
	  || apt->pmax[apt->number][ntypes + j] < apt->chempot[ntypes + j])
	  error("composition node %d is out of bounds.\n", j + 1);
      }

      /* check compnodes for valid values */
      if (ntypes == 2) {
	for (j = 0; j < compnodes; j++)
	  if (compnodelist[j] > 1 || compnodelist[j] < 0)
	    error("Composition node %d is %f but should be inside [0,1].\n",
	      j + 1, compnodelist[j]);
      }
    }
    if (compnodes != -1)
      printf("Enabled chemical potentials with %d extra composition node(s).\n",
	compnodes);
    if (compnodes == -1)
      compnodes = 0;
#endif
  }
#endif

#ifdef COULOMB
  fsetpos(infile, &startpos);
  /* skip to electrostatic section */
  do {
    fgetpos(infile, &filepos);
    fscanf(infile, "%s", buffer);
  } while (strcmp(buffer, "elstat") != 0 && !feof(infile));

  /* check for elstat keyword */
  if (strcmp("elstat", buffer) != 0) {
    error("No elstat option found in %s.\n", filename);
  }

  /* read electrostatic parameters */
  fscanf(infile, " %s", buffer);
  if (strcmp("ratio", buffer) != 0) {
    error("Could not read ratio");
  }
  for (i = 0; i < ntypes; i++) {
    if (1 > fscanf(infile, "%lf", &apt->ratio[i])) {
      error("Could not read ratio for atomtype #%d\n", i);
    }
  }
  for (i = 0; i < ntypes - 1; i++) {
    apt->param_name[apt->number][i] = (char *)malloc(30 * sizeof(char));
    if (4 > fscanf(infile, "%s %lf %lf %lf", apt->param_name[apt->number][i],
	&apt->charge[i], &apt->pmin[apt->number][i],
	&apt->pmax[apt->number][i])) {
      error("Could not read charge for atomtype #%d\n", i);
    }
    apt->invar_par[apt->number][i] = 0;
    if (apt->pmin[apt->number][i] == apt->pmax[apt->number][i]) {
      apt->invar_par[apt->number][i]++;
    }
    reg_for_free(apt->param_name[apt->number][i], "apt->param_name[%d][%d]",
      apt->number, i);
  }
#endif
#ifdef DIPOLE
  for (i = 0; i < ntypes; i++) {
    apt->param_name[apt->number + 1][i] = (char *)malloc(30 * sizeof(char));
    if (4 > fscanf(infile, "%s %lf %lf %lf",
	apt->param_name[apt->number + 1][i], &apt->dp_alpha[i],
	&apt->pmin[apt->number + 1][i], &apt->pmax[apt->number + 1][i])) {
      error("Could not read polarisability for atomtype #%d\n", i);
    }
    apt->invar_par[apt->number + 1][i] = 0;
    if (apt->pmin[apt->number + 1][i] == apt->pmax[apt->number + 1][i]) {
      apt->invar_par[apt->number + 1][i]++;
    }
    reg_for_free(apt->param_name[apt->number + 1][i], "apt->param_name[%d][%d]",
      apt->number + 1, i);
  }
  for (i = 0; i < apt->number; i++) {
    apt->param_name[apt->number + 2][i] = (char *)malloc(30 * sizeof(char));
    if (4 > fscanf(infile, "%s %lf %lf %lf",
	apt->param_name[apt->number + 2][i], &apt->dp_b[i],
	&apt->pmin[apt->number + 2][i], &apt->pmax[apt->number + 2][i])) {
      error("Could not read parameter dp_b for potential #%d\n", i);
    }
    apt->invar_par[apt->number + 2][i] = 0;
    if (apt->pmin[apt->number + 2][i] == apt->pmax[apt->number + 2][i]) {
      apt->invar_par[apt->number + 2][i]++;
    }
    reg_for_free(apt->param_name[apt->number + 2][i], "apt->param_name[%d][%d]",
      apt->number + 2, i);
  }
  for (i = 0; i < apt->number; i++) {
    apt->param_name[apt->number + 3][i] = (char *)malloc(30 * sizeof(char));
    if (4 > fscanf(infile, "%s %lf %lf %lf",
	apt->param_name[apt->number + 3][i], &apt->dp_c[i],
	&apt->pmin[apt->number + 3][i], &apt->pmax[apt->number + 3][i])) {
      error("Could not read parameter dp_c for potential #%d\n", i);
    }
    apt->invar_par[apt->number + 3][i] = 0;
    if (apt->pmin[apt->number + 3][i] == apt->pmax[apt->number + 3][i]) {
      apt->invar_par[apt->number + 3][i]++;
    }
    reg_for_free(apt->param_name[apt->number + 3][i], "apt->param_name[%d][%d]",
      apt->number + 3, i);
  }
#endif

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
      error("Premature end of potential file %s", filename);
    have_globals = 1;
    apt->total_par += apt->globals;

    i = apt->number + enable_cp;
    j = apt->globals;
    global_pot = i;

    /* allocate memory for global parameters */
    apt->names =
      (char **)realloc(apt->names, (global_pot + 1) * sizeof(char *));
    apt->names[global_pot] = (char *)malloc(20 * sizeof(char));
    strcpy(apt->names[global_pot], "global parameters");

    apt->n_glob = (int *)malloc(apt->globals * sizeof(int));
    reg_for_free(apt->n_glob, "apt->n_glob");

    apt->global_idx = (int ***)malloc(apt->globals * sizeof(int **));
    reg_for_free(apt->global_idx, "apt->global_idx");

    apt->values =
      (real **)realloc(apt->values, (global_pot + 1) * sizeof(real *));
    apt->values[global_pot] = (real *)malloc(j * sizeof(real));
    reg_for_free(apt->values[global_pot], "apt->values[%d]", global_pot);

    apt->invar_par =
      (int **)realloc(apt->invar_par, (global_pot + 1) * sizeof(int *));
    apt->invar_par[global_pot] = (int *)malloc((j + 1) * sizeof(int));
    reg_for_free(apt->invar_par[global_pot], "apt->invar_par[%d]", global_pot);

    apt->pmin = (real **)realloc(apt->pmin, (global_pot + 1) * sizeof(real *));
    apt->pmin[global_pot] = (real *)malloc(j * sizeof(real));
    reg_for_free(apt->pmin[global_pot], "apt->pmin[%d]", global_pot);

    apt->pmax = (real **)realloc(apt->pmax, (global_pot + 1) * sizeof(real *));
    apt->pmax[global_pot] = (real *)malloc(j * sizeof(real));
    reg_for_free(apt->pmax[global_pot], "apt->pmax[%d]", global_pot);

    apt->param_name =
      (char ***)realloc(apt->param_name, (global_pot + 1) * sizeof(char **));
    apt->param_name[global_pot] = (char **)malloc(j * sizeof(char *));
    reg_for_free(apt->param_name[global_pot], "apt->param_name[%d]",
      global_pot);

    pt->first = (int *)realloc(pt->first, (global_pot + 1) * sizeof(int));

    if (NULL == apt->values[global_pot] || NULL == apt->pmin[global_pot]
      || NULL == apt->n_glob || NULL == apt->global_idx
      || NULL == apt->pmax[global_pot]
      || NULL == apt->param_name[global_pot])
      error("Cannot allocate memory for global paramters.");

    /* initialize properly */
    apt->invar_par[i][j] = 0;

    /* read the global parameters */
    for (j = 0; j < apt->globals; j++) {
      apt->param_name[global_pot][j] = (char *)malloc(30 * sizeof(char));
      reg_for_free(apt->param_name[global_pot][j], "apt->param_name[%d][%d]",
	global_pot, j);

      if (NULL == apt->param_name[global_pot][j])
	error("Error in allocating memory for global parameter name");

      strcpy(apt->param_name[global_pot][j], "\0");
      ret_val =
	fscanf(infile, "%s %lf %lf %lf", apt->param_name[global_pot][j],
	&apt->values[global_pot][j], &apt->pmin[global_pot][j],
	&apt->pmax[global_pot][j]);
      if (4 > ret_val)
	if (strcmp(apt->param_name[global_pot][j], "type") == 0)
	  error
	    ("Not enough global parameters!\nYou specified %d parameter(s), but needed are %d.\nAborting",
	    j, apt->globals);

      /* check for duplicate names */
      for (k = j - 1; k >= 0; k--)
	if (strcmp(apt->param_name[global_pot][j],
	    apt->param_name[global_pot][k]) == 0) {
	  fprintf(stderr, "\nFound duplicate global parameter name!\n");
	  fprintf(stderr, "Parameter #%d (%s) is the same as #%d (%s)\n", j + 1,
	    apt->param_name[global_pot][j], k + 1,
	    apt->param_name[global_pot][k]);
	  error("Aborting");
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
	  warning
	    ("Starting value for gloabl paramter #%d is outside of specified adjustment range.\nResetting it to %f.\n",
	    j + 1, apt->values[i][j]);
	  if (apt->values[i][j] == 0)
	    fprintf(stderr,
	      "New value is >> 0 << ! Please be careful about this.\n");
	}
      }
    }
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
      error("Premature end of potential file %s", filename);
    if (strcmp(buffer, "type") != 0)
      error("Unknown keyword in file %s, expected \"type\" but found \"%s\".",
	filename, buffer);

    /* split name and _sc */
    token = strrchr(name, '_');
    if (token != NULL && strcmp(token + 1, "sc") == 0) {
      strncpy(buffer, name, strlen(name) - 3);
      buffer[strlen(name) - 3] = '\0';
      strcpy(name, buffer);
      smooth_pot[i] = 1;
    }

    if (apot_parameters(name) == -1)
      error
	("Unknown function type in file %s, please define \"%s\" in functions.c.",
	filename, name);

    strcpy(apt->names[i], name);
    apt->n_par[i] = apot_parameters(name);

    /* add one parameter for cutoff function if _sc is found */
    if (smooth_pot[i] == 1)
      apt->n_par[i]++;
    apt->total_par += apt->n_par[i];

    /* read cutoff */
    if (2 > fscanf(infile, "%s %lf", buffer, &apt->end[i]))
      error("Could not read cutoff for potential #%d in file %s\nAborting", i,
	filename);
    if (strcmp(buffer, "cutoff") != 0)
      error
	("No cutoff found for the %d. potential (%s) after \"type\" in file %s.\nAborting",
	i + 1, apt->names[i], filename);

    /* set very small begin, needed for EAM embedding function */
    apt->begin[i] = .0001;

    /* allocate memory for this parameter */
    apt->values[i] = (real *)malloc(apt->n_par[i] * sizeof(real));
    reg_for_free(apt->values[i], "apt->values[%d]", i);

    apt->invar_par[i] = (int *)malloc((apt->n_par[i] + 1) * sizeof(int));
    reg_for_free(apt->invar_par[i], "apt->invar_par[%d]", i);

    apt->pmin[i] = (real *)malloc(apt->n_par[i] * sizeof(real));
    reg_for_free(apt->pmin[i], "apt->pmin[%d]", i);

    apt->pmax[i] = (real *)malloc(apt->n_par[i] * sizeof(real));
    reg_for_free(apt->pmax[i], "apt->pmax[%d]", i);

    apt->param_name[i] = (char **)malloc(apt->n_par[i] * sizeof(char *));
    reg_for_free(apt->param_name[i], "apt->param_name[%d]", i);

    if (NULL == apt->values[i] || NULL == apt->pmin[i]
      || NULL == apt->pmax[i] || NULL == apt->param_name[i])
      error("Cannot allocate memory for potential paramters.\nAborting");

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
	error("Error in allocating memory for parameter name");
      reg_for_free(apt->param_name[i][j], "apt->param_name[%d][%d]", i, j);
      strcpy(apt->param_name[i][j], "\0");
      fgetpos(infile, &filepos);
      ret_val =
	fscanf(infile, "%s %lf %lf %lf", buffer, &apt->values[i][j],
	&apt->pmin[i][j], &apt->pmax[i][j]);
      strncpy(apt->param_name[i][j], buffer, 30);

      /* if last char of name is "!" we have a global parameter */
      if (strrchr(apt->param_name[i][j], '!') != NULL) {
	apt->param_name[i][j][strlen(apt->param_name[i][j]) - 1] = '\0';
	l = -1;
	for (k = 0; k < apt->globals; k++) {
	  if (strcmp(apt->param_name[i][j], apt->param_name[global_pot][k])
	    == 0)
	    l = k;
	}
	if (l == -1)
	  error("Could not find global parameter %s!\n", apt->param_name[i][j]);
	sprintf(apt->param_name[i][j], "%s!", apt->param_name[i][j]);

	/* write index array for global parameters */
	if (++apt->n_glob[l] > 1) {
	  apt->global_idx[l] =
	    (int **)realloc(apt->global_idx[l], apt->n_glob[l] * sizeof(int *));
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
	    if (strcmp(apt->param_name[i][j], "type") == 0 || feof(infile)) {
	      warning
		("No cutoff parameter given for potential #%d: adding one parameter.",
		i);
	      strcpy(apt->param_name[i][j], "h");
	      apt->values[i][j] = 1;
	      apt->pmin[i][j] = 0.5;
	      apt->pmax[i][j] = 2;
	      fsetpos(infile, &filepos);
	    }
	  } else {
	    if (strcmp(apt->param_name[i][j], "type") == 0)
	      error
		("Not enough parameters for potential #%d (%s) in file %s!\nYou specified %d parameters, but needed are %d.",
		i + 1, apt->names[i], filename, j, apt->n_par[i]);
	    error("Could not read parameter #%d of potential #%d in file %s",
	      j + 1, i + 1, filename);
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
	    warning
	      ("Starting value for paramter #%d in potential #%d is outside of specified adjustment range.\nResetting it to %f.\n",
	      j + 1, i + 1, apt->values[i][j]);
	    if (apt->values[i][j] == 0)
	      fprintf(stderr,
		"New value is >> 0 << ! Please be careful about this.\n");
	  }
	}
      }
    }
  }

#ifdef COULOMB
  apt->total_ne_par = apt->total_par;
#endif

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
    error("Could not assign the function pointers.\n");
#ifdef PAIR
  if (enable_cp) {
    cp_start = apt->total_par - apt->globals + ntypes * (ntypes + 1);
    apt->total_par += (ntypes + compnodes);
  }
#endif

#ifdef COULOMB
  apt->total_par += ntypes - 1;
#endif
#ifdef DIPOLE
  apt->total_par += ntypes;
  apt->total_par += (2 * apt->number);
#endif

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
#endif
#ifdef COULOMB
  pt->len += ntypes - 1;
#endif
#ifdef DIPOLE
  pt->len += ntypes;
  pt->len += (2 * apt->number);
#endif

  pt->table = (real *)malloc(pt->len * sizeof(real));
  reg_for_free(pt->table, "pt->table");

  calc_list = (real *)malloc(pt->len * sizeof(real));
  reg_for_free(calc_list, "calc_list");

  pt->idx = (int *)malloc(pt->len * sizeof(int));
  reg_for_free(pt->idx, "pt->idx");

  apt->idxpot = (int *)malloc(apt->total_par * sizeof(int));
  reg_for_free(apt->idxpot, "apt->idxpot");

  apt->idxparam = (int *)malloc(apt->total_par * sizeof(int));
  reg_for_free(apt->idxparam, "apt->idxparam");

  if ((NULL == pt->table) || (NULL == pt->idx) || (apt->idxpot == NULL)
    || (apt->idxparam == NULL))
    error("Cannot allocate memory for potential table.\n");
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
      pt->idx[k] = l++;
      apt->idxpot[k] = i;
      apt->idxparam[k++] = j;
      val++;
    }
    pt->idxlen += (ntypes + compnodes);
    global_idx += (ntypes + compnodes);
  }
#endif
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
  pt->idxlen += ntypes - 1;
#endif
#ifdef DIPOLE
  i = apt->number + 1;
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
  for (i = apt->number + 2; i < apt->number + 4; i++) {
    for (j = 0; j < (apt->number); j++) {
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
  pt->idxlen += (2 * apt->number);
#endif

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
    warning("Gauge degrees of freedom are NOT fixed!");
#endif /* NOPUNISH */

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
 ****************************************************************/

void read_pot_table3(pot_table_t *pt, int size, int ncols, int *nvals,
  char *filename, FILE *infile)
{
  int   i, j, k, l;
  real *val;

  /* read the info block of the function table */
  for (i = 0; i < size; i++) {
    if (3 > fscanf(infile, "%lf %lf %d", &pt->begin[i], &pt->end[i], &nvals[i]))
      error("Premature end of potential file %s", filename);
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
  pt->table = (real *)malloc(pt->len * sizeof(real));
  reg_for_free(pt->table, "pt->table");
  pt->xcoord = (real *)malloc(pt->len * sizeof(real));
  reg_for_free(pt->xcoord, "pt->xcoord");
  pt->d2tab = (real *)malloc(pt->len * sizeof(real));
  reg_for_free(pt->d2tab, "pt->d2tab");
  pt->idx = (int *)malloc(pt->len * sizeof(int));
  reg_for_free(pt->idx, "pt->idx");
  for (i = 0; i < pt->len; i++) {
    pt->table[i] = 0.;
    pt->xcoord[i] = 0.;
    pt->d2tab[i] = 0.;
    pt->idx[i] = 0.;
  }
  if ((NULL == pt->table) || (NULL == pt->idx) || (NULL == pt->d2tab))
    error("Cannot allocate memory for potential table");

  /* input loop */
  val = pt->table;
  k = 0;
  l = 0;

  /* read pair potentials */
  for (i = 0; i < ncols; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error("Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.;
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
	error("Premature end of potential file %s", filename);
      else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;
    }
  }

#if defined EAM || defined ADP
  /* read EAM transfer function rho(r) */
  for (i = ncols; i < ncols + ntypes; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error("Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.;
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
	error("Premature end of potential file %s", filename);
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
  for (i = ncols + ntypes; i < ncols + 2 * ntypes; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error("Premature end of potential file %s", filename);
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
	error("Premature end of potential file %s", filename);
      else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if (!invar_pot[i])
	pt->idx[k++] = l++;
      else
	l++;
    }
  }
#endif /* EAM || ADP */

#ifdef ADP
  /* read ADP dipole function u(r) */
  for (i = ncols + 2 * ntypes; i < 2 * (ncols + ntypes); i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error("Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.;
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
	error("Premature end of potential file %s", filename);
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
  for (i = 2 * (ncols + ntypes); i < 3 * ncols + 2 * ntypes; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error("Premature end of potential file %s", filename);
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
	error("Premature end of potential file %s", filename);
      else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;
    }
  }
#endif /* EAM || ADP */

  pt->idxlen = k;
  init_calc_table(pt, &calc_pot);

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

void read_pot_table4(pot_table_t *pt, int size, int ncols, int *nvals,
  char *filename, FILE *infile)
{
  int   i, k, l, j;
  real *val, *ord;

  /* read the info block of the function table */
  for (i = 0; i < size; i++) {
    if (1 > fscanf(infile, "%d", &nvals[i]))
      error("Premature end of potential file %s", filename);
    pt->step[i] = 0.;
    pt->invstep[i] = 0.;
    if (i == 0)
      pt->first[i] = 2;
    else
      pt->first[i] = pt->last[i - 1] + 3;
    pt->last[i] = pt->first[i] + nvals[i] - 1;
    pt->len = pt->first[i] + nvals[i];
  }
  /* allocate the function table */
  pt->table = (real *)malloc(pt->len * sizeof(real));
  pt->xcoord = (real *)malloc(pt->len * sizeof(real));
  pt->d2tab = (real *)malloc(pt->len * sizeof(real));
  pt->idx = (int *)malloc(pt->len * sizeof(int));
  if ((NULL == pt->table) || (NULL == pt->xcoord) || (NULL == pt->idx)
    || (NULL == pt->d2tab))
    error("Cannot allocate memory for potential table");

  /* input loop */
  val = pt->table;
  ord = pt->xcoord;
  k = 0;
  l = 0;

  /* read pair potentials */
  for (i = 0; i < ncols; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error("Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.;
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
	error("Premature end of potential file %s", filename);
      else {
	val++;
	ord++;
      }
      if ((j > 0) && (*(ord - 1) <= *(ord - 2)))
	error("Abscissa not monotonous in potential %d.", i);
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;

    }
    pt->begin[i] = pt->xcoord[pt->first[i]];
    pt->end[i] = pt->xcoord[pt->last[i]];
    /* pt->step is average step length.. */
    pt->step[i] = (pt->end[i] - pt->begin[i]) / ((real)nvals[i] - 1);
    pt->invstep[i] = 1. / pt->step[i];

  }
#if defined EAM || defined ADP
  /* read EAM transfer function rho(r) */
  for (i = ncols; i < ncols + ntypes; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error("Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.;
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
	error("Premature end of potential file %s", filename);
      else {
	ord++;
	val++;
      }
      if ((j > 0) && (*(ord - 1) <= *(ord - 2)))
	error("Abscissa not monotonous in potential %d.", i);
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;
    }
    pt->begin[i] = pt->xcoord[pt->first[i]];
    pt->end[i] = pt->xcoord[pt->last[i]];
    /* pt->step is average step length.. */
    pt->step[i] = (pt->end[i] - pt->begin[i]) / ((real)nvals[i] - 1);
    pt->invstep[i] = 1. / pt->step[i];
  }

  /* read EAM embedding function F(n) */
  for (i = ncols + ntypes; i < ncols + 2 * ntypes; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error("Premature end of potential file %s", filename);
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
	error("Premature end of potential file %s", filename);
      else {
	ord++;
	val++;
      }
      if ((j > 0) && (*(ord - 1) <= *(ord - 2)))
	error("Abscissa not monotonous in potential %d.", i);
      if (!invar_pot[i])
	pt->idx[k++] = l++;
      else
	l++;
    }
    pt->begin[i] = pt->xcoord[pt->first[i]];
    pt->end[i] = pt->xcoord[pt->last[i]];
    /* pt->step is average step length.. */
    pt->step[i] = (pt->end[i] - pt->begin[i]) / ((real)nvals[i] - 1);
    pt->invstep[i] = 1. / pt->step[i];
  }
#endif /* EAM || ADP */

#ifdef ADP
  /* read ADP dipole function u(r) */
  for (i = ncols + 2 * ntypes; i < 2 * (ncols + ntypes); i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error("Premature end of potential file %s", filename);
    } else {
      *val = 1e30;
      *(val + 1) = 0.;
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
	error("Premature end of potential file %s", filename);
      else {
	ord++;
	val++;
      }
      if ((j > 0) && (*(ord - 1) <= *(ord - 2)))
	error("Abscissa not monotonous in potential %d.", i);
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;
    }
    pt->begin[i] = pt->xcoord[pt->first[i]];
    pt->end[i] = pt->xcoord[pt->last[i]];
    /* pt->step is average step length.. */
    pt->step[i] = (pt->end[i] - pt->begin[i]) / ((real)nvals[i] - 1);
    pt->invstep[i] = 1. / pt->step[i];
  }

  /* read adp quadrupole function w(r) */
  for (i = 2 * (ncols + ntypes); i < 3 * ncols + 2 * ntypes; i++) {
    if (have_grad) {
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1))
	error("Premature end of potential file %s", filename);
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
	error("Premature end of potential file %s", filename);
      else {
	ord++;
	val++;
      }
      if ((j > 0) && (*(ord - 1) <= *(ord - 2)))
	error("Abscissa not monotonous in potential %d.", i);
      if (!invar_pot[i])
	pt->idx[k++] = l++;
      else
	l++;
    }
    pt->begin[i] = pt->xcoord[pt->first[i]];
    pt->end[i] = pt->xcoord[pt->last[i]];
    /* pt->step is average step length.. */
    pt->step[i] = (pt->end[i] - pt->begin[i]) / ((real)nvals[i] - 1);
    pt->invstep[i] = 1. / pt->step[i];
  }
#endif /* ADP */

  pt->idxlen = k;
  init_calc_table(pt, &calc_pot);
}

#endif /* APOT */

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
  real *val, f, h;
#endif

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
	  calct->step = (real *)malloc(size * sizeof(real));
	  reg_for_free(calct->step, "calct->step");
	  calct->invstep = (real *)malloc(size * sizeof(real));
	  reg_for_free(calct->invstep, "calct->invstep");
	  calct->xcoord = (real *)malloc(calct->len * sizeof(real));
	  reg_for_free(calct->xcoord, "calct->xcoord");
	  calct->table = (real *)malloc(calct->len * sizeof(real));
	  reg_for_free(calct->table, "calct->table");
	  calct->d2tab = (real *)malloc(calct->len * sizeof(real));
	  reg_for_free(calct->d2tab, "calct->d2tab");
	  calct->idx = (int *)malloc(calct->len * sizeof(int));
	  reg_for_free(calct->idx, "calct->idx");
	  if (calct->first == NULL || calct->last == NULL || calct->step == NULL
	    || calct->invstep == NULL || calct->xcoord == NULL
	    || calct->table == NULL || calct->d2tab == NULL
	    || calct->idx == NULL)
	    error("Cannot allocate info block for calc potential table\n");

	  /* initialize the calc_pot table */
	  for (i = 0; i < size; i++) {
	    val = apot_table.values[i];
	    h = apot_table.values[i][apot_table.n_par[i] - 1];
	    calct->table[i * APOT_STEPS + i * 2] = 10e30;
	    calct->table[i * APOT_STEPS + i * 2 + 1] = 0;
	    calct->first[i] = (x += 2);
	    calct->last[i] = (x += APOT_STEPS - 1);
	    x++;
	    calct->step[i] =
	      (calct->end[i] - calct->begin[i]) / (APOT_STEPS - 1);
	    calct->invstep[i] = 1. / calct->step[i];
	    for (j = 0; j < APOT_STEPS; j++) {
	      index = i * APOT_STEPS + (i + 1) * 2 + j;
	      calct->xcoord[index] = calct->begin[i] + j * calct->step[i];

	      apot_table.fvalue[i] (calct->xcoord[index], val, &f);
	      calct->table[index] =
		smooth_pot[i] ? f * cutoff(calct->xcoord[index],
		calct->begin[i], h) : f;
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
}

void update_calc_table(real *xi_opt, real *xi_calc, int do_all)
{
#ifdef APOT
  char  msg[255];
  int   i, j, k, m, n, change;
  real  f, h = 0;
  real *list, *val;
#endif

  switch (format) {
#ifdef APOT
      case 0:
	{
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
	    if (smooth_pot[i] && !invar_pot[i]) {
	      h = *(val + 1 + apot_table.n_par[i]);
	      if (h == 0)
		error("The cutoff parameter for potential %d is 0!", i);
	    }

	    (*val) =
	      apot_grad(calc_pot.begin[i], val + 2, apot_table.fvalue[i]);
	    val += 2;
	    /* check if something has changed */
	    change = 0;
	    for (j = 0; j < apot_table.n_par[i]; j++) {
	      if (list[j] != val[j]) {
		change = 1;
		list[j] = val[j];
	      }
	    }
	    if ((change || do_all) && !invar_pot[i]) {
	      for (j = 0; j < APOT_STEPS; j++) {
		k = i * APOT_STEPS + (i + 1) * 2 + j;
		apot_table.fvalue[i] (calc_pot.xcoord[k], val, &f);
		*(xi_calc + k) =
		  smooth_pot[i] ? f * cutoff(calc_pot.xcoord[k],
		  apot_table.end[i], h) : f;
#ifndef COULOMB
		if (isnan(f) || isnan(*(xi_calc + k))) {
		  sprintf(msg, "Potential value was nan or inf. Aborting.\n");
		  fprintf(stderr, "This occured in potential %d (%s)\n", i,
		    apot_table.names[i]);
		  fprintf(stderr, "at distance r=%f with the parameters:\n",
		    calc_pot.xcoord[k]);
		  for (m = 0; m < apot_table.n_par[i]; m++)
		    fprintf(stderr, "%s %f\n", apot_table.param_name[i][m],
		      *(val + m));
		  if (smooth_pot[i])
		    fprintf(stderr, "h %f\n", h);
		  error(msg);
		}
#endif
	      }
	    }
	    val += apot_table.n_par[i];
	    list += apot_table.n_par[i] + 2;
	  }
	}

	return;
#endif
      case 3:			/* fall through */
      case 4:
	return;
  }
}


#ifdef PARABEL

/****************************************************************
 *
 *  Evaluate value from parabole through three points.
 *  Extrapolates for all k.
 *
 ****************************************************************/

real parab_ed(pot_table_t *pt, real *xi, int col, real r)
{
  real  rr, istep, chi, p0, p1, p2, dv, d2v;
  int   k;

  /* renorm to beginning of table */
  rr = r - pt->begin[col];

  /* indices into potential table */
  istep = pt->invstep[col];
  k = 0;
  chi = rr * istep;
  k = pt->first[col];

  /* intermediate values */
  p0 = xi[k++];
  p1 = xi[k++];
  p2 = xi[k];
  dv = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* return the potential value */
  return p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;
}

/****************************************************************
 *
 *  Evaluate value from parabole through three points.
 *  Extrapolates for all k. Nonequidistant points.
 *
 ****************************************************************/

real parab_ne(pot_table_t *pt, real *xi, int col, real r)
{
  real  x0, x1, x2, chi0, chi1, chi2, p0, p1, p2;
  int   k;

  /* renorm to beginning of table */
  /* rr = r - pt->begin[col]; */
  k = pt->first[col];
  x0 = pt->xcoord[k];
  p0 = xi[k++];
  x1 = pt->xcoord[k];
  p1 = xi[k++];
  x2 = pt->xcoord[k];
  p2 = xi[k];

  /* indices into potential table */
  chi0 = (r - x0) / (x2 - x1);
  chi1 = (r - x1) / (x2 - x0);
  chi2 = (r - x2) / (x1 - x0);

  /* intermediate values */
  /* dv  = p1 - p0; */
  /* d2v = p2 - 2 * p1 + p0; */

  /* return the potential value */
  return chi1 * chi2 * p0 - chi0 * chi2 * p1 + chi0 * chi1 * p2;

}

/****************************************************************
 *
 *  Evaluate deritvative from parabole through three points.
 *  Extrapolates for all k.
 *
 ****************************************************************/

real parab_grad_ed(pot_table_t *pt, real *xi, int col, real r)
{
  real  rr, istep, chi, p0, p1, p2, dv, d2v;
  int   k;

  /* renorm to beginning of table */
  rr = r - pt->begin[col];

  /* indices into potential table */
  istep = pt->invstep[col];
  k = 0;
  chi = rr * istep;
  k = pt->first[col];

  /* intermediate values */
  p0 = xi[k++];
  p1 = xi[k++];
  p2 = xi[k];
  dv = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* return the derivative */
  return istep * (dv + (chi - 0.5) * d2v);
}

/****************************************************************
 *
 *  Evaluate deritvative from parabole through three points.
 *  Extrapolates for all k.
 *
 ****************************************************************/

real parab_grad_ne(pot_table_t *pt, real *xi, int col, real r)
{
  real  h0, h1, h2, x0, x1, x2, chi0, chi1, chi2, p0, p1, p2;
  int   k;

  /* renorm to beginning of table */
  k = pt->first[col];
  x0 = pt->xcoord[k];
  p0 = xi[k++];
  x1 = pt->xcoord[k];
  p1 = xi[k++];
  x2 = pt->xcoord[k];
  p2 = xi[k];

  h0 = x2 - x1;
  h1 = x2 - x0;
  h2 = x1 - x0;

  chi0 = (r - x0) / h0;
  chi1 = (r - x1) / h1;
  chi2 = (r - x2) / h2;

  /* return the potential value */
  return (chi2 / h1 + chi1 / h2) * p0 - (chi0 / h2 + chi2 / h0) * p1 +
    (chi0 / h1 + chi1 / h0) * p2;

}

/****************************************************************
 *
 *  Evaluate value and deritvative from parabole through three points.
 *  Extrapolates for all k.
 *
 ****************************************************************/

real parab_comb_ed(pot_table_t *pt, real *xi, int col, real r, real *grad)
{
  real  rr, istep, chi, p0, p1, p2, dv, d2v;
  int   k;

  /* renorm to beginning of table */
  rr = r - pt->begin[col];

  /* indices into potential table */
  istep = pt->invstep[col];
  k = 0;
  chi = rr * istep;
  k = pt->first[col];

  /* intermediate values */
  p0 = xi[k++];
  p1 = xi[k++];
  p2 = xi[k];
  dv = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* set the derivative */
  *grad = istep * (dv + (chi - 0.5) * d2v);
  /* return the potential value */
  return p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;
}

/***************************************************************
 *
 *  Evaluate value and deritvative from parabole through three points.
 *  Extrapolates for all k.
 *
 ***************************************************************/

real parab_comb_ne(pot_table_t *pt, real *xi, int col, real r, real *grad)
{
  real  h0, h1, h2, x0, x1, x2, chi0, chi1, chi2, p0, p1, p2;
  int   k;

  /* renorm to beginning of table */
  k = pt->first[col];
  x0 = pt->xcoord[k];
  p0 = xi[k++];
  x1 = pt->xcoord[k];
  p1 = xi[k++];
  x2 = pt->xcoord[k];
  p2 = xi[k];

  h0 = x2 - x1;
  h1 = x2 - x0;
  h2 = x1 - x0;

  chi0 = (r - x0) / h0;
  chi1 = (r - x1) / h1;
  chi2 = (r - x2) / h2;

  /* return the potential value */
  *grad =
    (chi2 / h1 + chi1 / h2) * p0 - (chi0 / h2 + chi2 / h0) * p1 + (chi0 / h1 +
    chi1 / h0) * p2;

  return chi1 * chi2 * p0 - chi0 * chi2 * p1 + chi0 * chi1 * p2;
}

#endif /* PARABEL */

#ifdef COULOMB

/****************************************************************
 *
 *  calculate tail of coulomb-potential and its first derivative
 *
 ****************************************************************/

void init_tails()
{
  int   i, j;

  for (i = 0; i < natoms; i++)
    for (j = 0; j < atoms[i].n_neigh; j++)
      elstat_shift(atoms[i].neigh[j].r, &atoms[i].neigh[j].fnval_el,
	&atoms[i].neigh[j].grad_el, &atoms[i].neigh[j].ggrad_el);
}

/****************************************************************
 *
 * write coulomb-potential
 *
 ****************************************************************/

void write_coulomb_table()
{
  apot_table_t *apt = &apot_table;
  if (ntypes == 2) {
    int   i, j;
    real  value, c1;
    FILE *outfile;
    char *filename1 = "Coulomb_00";
    char *filename2 = "Coulomb_01";
    char *filename3 = "Coulomb_11";

    c1 = -apt->charge[0] / 2;

    outfile = fopen(filename1, "a");
    for (i = 0; i < natoms; i++) {
      for (j = 0; j < atoms[i].n_neigh; j++) {
	value = apt->charge[0] * apt->charge[0] * atoms[i].neigh[j].fnval_el;
	fprintf(outfile, "%f\t%f\n", atoms[i].neigh[j].r, value);
      }
    }
    fclose(outfile);

    outfile = fopen(filename2, "a");
    for (i = 0; i < natoms; i++) {
      for (j = 0; j < atoms[i].n_neigh; j++) {
	value = apt->charge[0] * c1 * atoms[i].neigh[j].fnval_el;
	fprintf(outfile, "%f\t%f\n", atoms[i].neigh[j].r, value);
      }
    }
    fclose(outfile);

    outfile = fopen(filename3, "a");
    for (i = 0; i < natoms; i++) {
      for (j = 0; j < atoms[i].n_neigh; j++) {
	value = c1 * c1 * atoms[i].neigh[j].fnval_el;
	fprintf(outfile, "%f\t%f\n", atoms[i].neigh[j].r, value);
      }
    }
    fclose(outfile);
  } else {
    printf("Coulomb-outfiles are only available in case of two atom types.");
  }
}

#endif /* COULOMB */

#ifdef APOT

/****************************************************************
 *
 *  write potential table (format 0)
 *
 ****************************************************************/

void write_pot_table0(apot_table_t *apt, char *filename)
{
  int   i, j;
  FILE *outfile;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error("Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 0 %d", apt->number);
  fprintf(outfile, "\n#T %s", interaction_name);
  if (have_elements) {
    fprintf(outfile, "\n#C");
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    fprintf(outfile, "\n##");
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);
#if defined EAM || defined ADP
    /* transfer functions */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    /* embedding functions */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
#endif
#if defined ADP
    /* dipole terms */
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);
    /* quadrupole terms */
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);
#endif
  }
  if (have_invar) {
    fprintf(outfile, "\n#I");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, " %d", invar_pot[i]);
  }
  fprintf(outfile, "\n#E\n\n");

#ifdef PAIR
  if (enable_cp) {
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, "cp_%s %.10f %.2f %.2f\n", elements[i], apt->chempot[i],
	apt->pmin[apt->number][i], apt->pmax[apt->number][i]);
    if (compnodes > 0)
      fprintf(outfile, "cn %d\n", compnodes);
    for (j = 0; j < compnodes; j++)
      fprintf(outfile, "%.2f %.10f %.2f %.2f\n", compnodelist[j],
	apt->chempot[ntypes + j], apt->pmin[apt->number][ntypes + j],
	apt->pmax[apt->number][ntypes + j]);
    fprintf(outfile, "\n");
  }
#endif

#ifdef COULOMB
  fprintf(outfile, "elstat\n");
  for (i = 0; i < ntypes - 1; i++)
    fprintf(outfile, "%s\t %f\t %f\t %f\n", apt->param_name[apt->number][i],
      apt->charge[i], apt->pmin[apt->number][i], apt->pmax[apt->number][i]);
  fprintf(outfile, "charge_%s\t %f\n", elements[ntypes - 1], apt->last_charge);
#ifdef DIPOLE
  for (i = 0; i < ntypes; i++)
    fprintf(outfile, "%s\t %f\t %f\t %f\n", apt->param_name[apt->number + 1][i],
      apt->dp_alpha[i], apt->pmin[apt->number + 1][i],
      apt->pmax[apt->number + 1][i]);
  for (i = 0; i < apt->number; i++) {
    fprintf(outfile, "%s\t %f\t %f\t %f\n", apt->param_name[apt->number + 2][i],
      apt->dp_b[i], apt->pmin[apt->number + 2][i],
      apt->pmax[apt->number + 2][i]);
  }
  for (i = 0; i < apt->number; i++) {
    fprintf(outfile, "%s\t %f\t %f\t %f\n", apt->param_name[apt->number + 3][i],
      apt->dp_c[i], apt->pmin[apt->number + 3][i],
      apt->pmax[apt->number + 3][i]);
  }
#endif
  fprintf(outfile, "\n");
#endif

  if (have_globals) {
    fprintf(outfile, "global %d\n", apt->globals);
    for (i = 0; i < apt->globals; i++)
      fprintf(outfile, "%s %.10f %.2f %.2f\n", apt->param_name[global_pot][i],
	apt->values[global_pot][i], apt->pmin[global_pot][i],
	apt->pmax[global_pot][i]);
    fprintf(outfile, "\n");
  }

  /* write data */
  for (i = 0; i < apt->number; i++) {
    if (smooth_pot[i]) {
      fprintf(outfile, "type %s_sc\n", apt->names[i]);
    } else {
      fprintf(outfile, "type %s\n", apt->names[i]);
    }
    fprintf(outfile, "cutoff\t %f\n", apot_table.end[i]);
    fprintf(outfile, "# rmin\t %f\n", apt->begin[i]);
    for (j = 0; j < apt->n_par[i]; j++) {
      if (apt->param_name[i][j][strlen(apt->param_name[i][j]) - 1] != '!') {
	fprintf(outfile, "%s\t %.10f\t %.2f\t %.2f\n", apt->param_name[i][j],
	  apt->values[i][j], apt->pmin[i][j], apt->pmax[i][j]);
      } else {
	fprintf(outfile, "%s\n", apt->param_name[i][j]);
      }
    }
    if (i != (apt->number - 1))
      fprintf(outfile, "\n");
  }
  fclose(outfile);
}

#else /* APOT */

/****************************************************************
 *
 *  write potential table (format 3)
 *
 ****************************************************************/

void write_pot_table3(pot_table_t *pt, char *filename)
{
  FILE *outfile = NULL, *outfile2 = NULL;
  int   i, j, flag = 0;
  real  r;

  if (*plotpointfile != '\0')
    flag = 1;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error("Could not open file %s\n", filename);

  /* if needed: open file for plotpoints */
  if (flag) {
    outfile2 = fopen(plotpointfile, "w");
    if (NULL == outfile)
      error("Could not open file %s\n", filename);
  }

  /* write header */
  fprintf(outfile, "#F 3 %d", pt->ncols);
  fprintf(outfile, "\n#T %s", interaction_name);
  if (have_elements) {
    fprintf(outfile, "\n#C");
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    fprintf(outfile, "\n##");
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);
#ifdef EAM
    /* transfer functions */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    /* embedding functions */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
#endif
  }
  if (have_invar) {
    fprintf(outfile, "\n#I");
    for (i = 0; i < pt->ncols; i++)
      fprintf(outfile, " %d", invar_pot[i]);
  }
  fprintf(outfile, "\n#G");
  for (i = 0; i < pt->ncols; i++)
    fprintf(outfile, " %d", gradient[i]);
  fprintf(outfile, "\n#E\n");

  /* write info block */
  for (i = 0; i < pt->ncols; i++) {
    fprintf(outfile, "%.16e %.16e %d\n", pt->begin[i], pt->end[i],
      pt->last[i] - pt->first[i] + 1);
  }
  fprintf(outfile, "\n");

  /* write data */
  for (i = 0; i < pt->ncols; i++) {
    r = pt->begin[i];
    /* write gradient */
    fprintf(outfile, "%.16e %.16e\n", pt->table[pt->first[i] - 2],
      pt->table[pt->first[i] - 1]);
    for (j = pt->first[i]; j <= pt->last[i]; j++) {
      fprintf(outfile, "%.16e\n", pt->table[j]);
      if (flag)
	fprintf(outfile2, "%.6e %.6e %d\n", r, pt->table[j], j);
      r += pt->step[i];
    }
    fprintf(outfile, "\n");
    if (flag)
      fprintf(outfile2, "\n\n");
  }
  fclose(outfile);
  if (flag)
    fclose(outfile2);
}

/****************************************************************
 *
 *  write potential table (format 4)
 *
 ****************************************************************/

void write_pot_table4(pot_table_t *pt, char *filename)
{
  FILE *outfile = NULL, *outfile2 = NULL;
  int   i, j, flag = 0;

  if (*plotpointfile != '\0')
    flag = 1;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error("Could not open file %s\n", filename);

  /* if needed: open file for plotpoints */
  if (flag) {
    outfile2 = fopen(plotpointfile, "w");
    if (NULL == outfile)
      error("Could not open file %s\n", filename);
  }

  /* write header */
  fprintf(outfile, "#F 4 %d", pt->ncols);
  fprintf(outfile, "\n#T %s", interaction_name);
  if (have_elements) {
    fprintf(outfile, "\n#C");
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    fprintf(outfile, "\n##");
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);
#ifdef EAM
    /* transfer functions */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    /* embedding functions */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
#endif
  }
  if (have_invar) {
    fprintf(outfile, "\n#I");
    for (i = 0; i < pt->ncols; i++)
      fprintf(outfile, " %d", invar_pot[i]);
  }
  fprintf(outfile, "\n#G");
  for (i = 0; i < pt->ncols; i++)
    fprintf(outfile, " %d", gradient[i]);
  fprintf(outfile, "\n#E\n");

  /* write info block */
  for (i = 0; i < pt->ncols; i++) {
    fprintf(outfile, "%d\n", pt->last[i] - pt->first[i] + 1);
  }
  fprintf(outfile, "\n");

  /* write data */
  for (i = 0; i < pt->ncols; i++) {
    fprintf(outfile, "%.16e %.16e\n", pt->table[pt->first[i] - 2],
      pt->table[pt->first[i] - 1]);
    for (j = pt->first[i]; j <= pt->last[i]; j++) {
      fprintf(outfile, "%.16e %.16e\n", pt->xcoord[j], pt->table[j]);
      if (flag)
	fprintf(outfile2, "%.6e %.6e %d\n", pt->xcoord[j], pt->table[j], j);
    }
    fprintf(outfile, "\n");
    if (flag)
      fprintf(outfile2, "\n\n");
  }
  fclose(outfile);
  if (flag)
    fclose(outfile2);
}

#endif /* APOT */

/****************************************************************
 *
 *  write potential table for IMD (format 2)
 *
 ****************************************************************/

void write_pot_table_imd(pot_table_t *pt, char *prefix)
{
  int   i, j, k, m, m2, col1, col2;
  real  r2, temp;
  real *r2begin, *r2end, *r2step;
#ifndef APOT
  real  temp2;
#endif
#if defined EAM || defined ADP
  real  root;
#endif /* EAM */
  FILE *outfile;
  char  filename[255];

  /* allocate memory */
  r2begin = (real *)malloc(ntypes * ntypes * sizeof(real));
  r2end = (real *)malloc(ntypes * ntypes * sizeof(real));
  r2step = (real *)malloc(ntypes * ntypes * sizeof(real));
  if ((r2begin == NULL) || (r2end == NULL) || (r2step == NULL))
    error("Cannot allocate memory in  write_pot_table_imd");

  /* pair potential part (over r^2) */
  sprintf(filename, "%s_phi.imd.pt", prefix);
  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error("Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes * ntypes);

  /* write info block */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col2 = i * ntypes + j;
      /* Extrapolation possible  */
#ifdef APOT
      r2begin[col2] = SQR((plotmin == 0 ? 0.1 : plotmin));
#else
      r2begin[col2] = SQR(MAX(pt->begin[col1] - extend * pt->step[col1], 0));
#endif
      r2end[col2] = SQR(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2],
	r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col2 = i * ntypes + j;
      r2 = r2begin[col2];
      for (k = 0; k < imdpotsteps; k++) {
#ifdef APOT
	apot_table.fvalue[col1] (sqrt(r2), apot_table.values[col1], &temp);
	temp =
	  smooth_pot[col1] ? temp * cutoff(sqrt(r2), apot_table.end[col1],
	  apot_table.values[col1][apot_table.n_par[col1] - 1]) : temp;
	fprintf(outfile, "%.16e\n", temp);
#else
#ifdef NEWSCALE
	/* Pair potentials corrected so that U'(1)   =0 with NORESCALE */
	/*                               and U'(n_av)=0 without */
	fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1,
	    sqrt(r2)) + (sqrt(r2) <=
	    pt->end[paircol + j] ? lambda[i] * splint_ne(pt, pt->table,
	      paircol + j,
	      sqrt(r2)) : 0.) + (sqrt(r2) <=
	    pt->end[paircol + i] ? lambda[j] * splint_ne(pt, pt->table,
	      paircol + i, sqrt(r2)) : 0.));
#else
	fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, sqrt(r2)));
#endif /* NEWSCALE */
#endif /* APOT */
	r2 += r2step[col2];
      }
      fprintf(outfile, "%.16e\n", 0.0);
      fprintf(outfile, "\n");
    }
  }
  fclose(outfile);
  printf("IMD pair potential written to \t\t%s\n", filename);

#if defined EAM || defined ADP
  /* write transfer function (over r^2) */
  sprintf(filename, "%s_rho.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error("Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes * ntypes);

  /* write info block */
  for (i = 0; i < ntypes; i++) {
    for (j = 0; j < ntypes; j++) {
      col1 = (ntypes * (ntypes + 1)) / 2 + j;
      col2 = i * ntypes + j;
#ifdef APOT
      r2begin[col2] = SQR((plotmin == 0 ? 0.1 : plotmin));
#else
      /* Extrapolation possible  */
      r2begin[col2] = SQR(MAX(pt->begin[col1] - extend * pt->step[col1], 0));
#endif /* APOT */
      r2end[col2] = SQR(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2],
	r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  for (i = 0; i < ntypes; i++) {
    for (j = 0; j < ntypes; j++) {
      col1 = (ntypes * (ntypes + 1)) / 2 + j;
      col2 = i * ntypes + j;
      r2 = r2begin[col2];
      for (k = 0; k < imdpotsteps; k++) {
#ifdef APOT
	apot_table.fvalue[col1] (sqrt(r2), apot_table.values[col1], &temp);
	temp =
	  smooth_pot[col1] ? temp * cutoff(sqrt(r2), apot_table.end[col1],
	  apot_table.values[col1][apot_table.n_par[col1] - 1]) : temp;
	fprintf(outfile, "%.16e\n", temp);
#else
	fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, sqrt(r2)));
#endif /* APOT */
	r2 += r2step[col2];
      }
      fprintf(outfile, "%.16e\n", 0.0);
      fprintf(outfile, "\n");
    }
  }
  fclose(outfile);
  printf("IMD transfer function written to \t%s\n", filename);

  /* write embedding function (over r) */
  sprintf(filename, "%s_F.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error("Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes);

  /* write info block */
  for (i = 0; i < ntypes; i++) {
    col1 = (ntypes * (ntypes + 3)) / 2 + i;
#ifdef APOT
    r2begin[i] = 0;
    r2end[i] = pt->end[col1];
#else
    /* pad with zeroes */
    r2begin[i] = pt->begin[col1] - extend * pt->step[col1];
    /* extrapolation */
    r2end[i] = pt->end[col1] + extend * pt->step[col1];
#endif
    r2step[i] = (r2end[i] - r2begin[i]) / imdpotsteps;
    fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[i], r2end[i], r2step[i]);
  }
  fprintf(outfile, "\n");

  /* write data */
  for (i = 0; i < ntypes; i++) {
    r2 = r2begin[i];
    col1 = (ntypes * (ntypes + 3)) / 2 + i;
    root =
      (pt->begin[col1] >
      0) ? pt->table[pt->first[col1]] / sqrt(pt->begin[col1]) : 0.;
    root +=
      (pt->end[col1] <
      0) ? pt->table[pt->last[col1]] / sqrt(-pt->end[col1]) : 0;
    for (k = 0; k <= imdpotsteps; k++) {
#ifdef APOT
      apot_table.fvalue[col1] (r2, apot_table.values[col1], &temp);
#else
#ifdef WZERO
      if (r2 < pt->begin[col1] && pt->begin[col1] > 0)
	if (r2 <= 0)
	  temp = 100 * (root / fabs(root)) * r2;	/* steep decline */
	else
	  temp = root * sqrt(r2);	/* sqrt-like shape */
      else if (r2 > pt->end[col1] && pt->end[col1] < 0)
	if (r2 >= 0)
	  temp = -100. * (root / fabs(root)) * r2;	/* steep decline */
	else
	  temp = root * sqrt(-r2);	/* sqrt-like shape */
      else {
#ifdef PARABEL
	temp = parab(pt, pt->table, col1, r2);
#else
	temp = splint_ne(pt, pt->table, col1, r2);
#endif
      }
#else /* WZERO */
      temp = splint_ne(pt, pt->table, col1, r2);
#endif /* WZERO */
      temp2 = r2 - pt->end[col1];
      temp += (temp2 > 0.) ? 5e2 * (temp2 * temp2 * temp2) : 0.;
#ifdef NEWSCALE
      temp -= lambda[i] * r2;
#endif /* NEWSCALE */
#endif /* APOT */
      fprintf(outfile, "%.16e\n", temp);
      r2 += r2step[i];
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);
  printf("IMD embedding function written to \t%s\n", filename);
#endif

#ifdef ADP
  /* write dipole function (over r^2) */
  sprintf(filename, "%s_upot.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error("Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes * ntypes);

  /* write info block */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col1 += paircol + 2 * ntypes;
      col2 = i * ntypes + j;
      /* Extrapolation possible  */
#ifdef APOT
      r2begin[col2] = SQR((plotmin == 0 ? 0.1 : plotmin));
#else
      r2begin[col2] = SQR(MAX(pt->begin[col1] - extend * pt->step[col1], 0));
#endif
      r2end[col2] = SQR(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2],
	r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col1 += paircol + 2 * ntypes;
      col2 = i * ntypes + j;
      r2 = r2begin[col2];
      for (k = 0; k < imdpotsteps; k++) {
#ifdef APOT
	apot_table.fvalue[col1] (sqrt(r2), apot_table.values[col1], &temp);
	temp =
	  smooth_pot[col1] ? temp * cutoff(sqrt(r2), apot_table.end[col1],
	  apot_table.values[col1][apot_table.n_par[col1] - 1]) : temp;
	fprintf(outfile, "%.16e\n", temp);
#else
	fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, sqrt(r2)));
#endif /* APOT */
	r2 += r2step[col2];
      }
      fprintf(outfile, "%.16e\n", 0.0);
      fprintf(outfile, "\n");
    }
  }
  fclose(outfile);
  printf("IMD dipole potential written to \t%s\n", filename);

  /* write quadrupole function (over r^2) */
  sprintf(filename, "%s_wpot.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error("Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes * ntypes);

  /* write info block */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col1 += 2 * paircol + 2 * ntypes;
      col2 = i * ntypes + j;
      /* Extrapolation possible  */
#ifdef APOT
      r2begin[col2] = SQR((plotmin == 0 ? 0.1 : plotmin));
#else
      r2begin[col2] = SQR(MAX(pt->begin[col1] - extend * pt->step[col1], 0));
#endif
      r2end[col2] = SQR(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2],
	r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col1 += 2 * paircol + 2 * ntypes;
      col2 = i * ntypes + j;
      r2 = r2begin[col2];
      for (k = 0; k < imdpotsteps; k++) {
#ifdef APOT
	apot_table.fvalue[col1] (sqrt(r2), apot_table.values[col1], &temp);
	temp =
	  smooth_pot[col1] ? temp * cutoff(sqrt(r2), apot_table.end[col1],
	  apot_table.values[col1][apot_table.n_par[col1] - 1]) : temp;
	fprintf(outfile, "%.16e\n", temp);
#else
	fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, sqrt(r2)));
#endif /* APOT */
	r2 += r2step[col2];
      }
      fprintf(outfile, "%.16e\n", 0.0);
      fprintf(outfile, "\n");
    }
  }
  fclose(outfile);
  printf("IMD quadrupole potential written to \t%s\n", filename);
#endif

  free(r2begin);
  free(r2end);
  free(r2step);
}

/****************************************************************
 *
 *  write plot version of potential table
 *
 ****************************************************************/

void write_plotpot_pair(pot_table_t *pt, char *filename)
{
  FILE *outfile;
  int   i, j;
#ifndef APOT
  int   k = 0, l;
#else
  real  h;
#endif
  real  r, r_step, temp;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error("Could not open file %s\n", filename);

  /* write data */
#ifndef APOT
  for (i = 0; i < ntypes; i++)
    for (j = i; j < ntypes; j++) {
      r = pt->begin[k];
      r_step = (pt->end[k] - r) / (NPLOT - 1);
      for (l = 0; l < NPLOT - 1; l++) {
#ifdef NEWSCALE
	fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, k, r)
	  + (r <= pt->end[paircol + i] ? splint_ne(pt, pt->table, paircol + i,
	      r) * lambda[j] : 0.)
	  + (r <= pt->end[paircol + j] ? splint_ne(pt, pt->table, paircol + j,
	      r) * lambda[i] : 0.));
#else
	fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, k, r));
#endif /* NEWSCALE */
	r += r_step;
      }
      fprintf(outfile, "%e %e\n\n\n", r, 0.0);
      k++;
    }
#if defined EAM || defined ADP
  for (i = paircol; i < paircol + ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
  for (i = paircol + ntypes; i < paircol + 2 * ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (l = 0; l < NPLOT; l++) {
#ifdef PARABEL
      temp = parab(pt, pt->table, i, r);
#else
      temp = splint_ne(pt, pt->table, i, r);
#endif
#ifdef NEWSCALE
      temp -= lambda[i - (paircol + ntypes)] * r;
#endif /* NEWSCALE */
      fprintf(outfile, "%e %e\n", r, temp);
      r += r_step;
    }
    fprintf(outfile, "\n\n\n");
  }
#endif /* EAM || ADP */
#ifdef ADP
  for (i = paircol + 2 * ntypes; i < 2 * (paircol + ntypes); i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - r) / (NPLOT - 1);
    for (l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
    k++;
  }
  for (i = 2 * (paircol + ntypes); i < 3 * paircol + 2 * ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - r) / (NPLOT - 1);
    for (l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
    k++;
  }
#endif
#else /* APOT */
  for (i = 0; i < apot_table.number; i++) {
    r = (plotmin == 0 ? 0.1 : plotmin);
    r_step = (apot_table.end[i] - r) / (NPLOT - 1);
    h = apot_table.values[i][apot_table.n_par[i] - 1];
    for (j = 0; j < NPLOT; j++) {
      apot_table.fvalue[i] (r, apot_table.values[i], &temp);
      temp = smooth_pot[i] ? temp * cutoff(r, apot_table.end[i], h) : temp;
      if (isnan(temp))
	temp = 10e30;
      fprintf(outfile, "%e %e\n", r, temp);
      r += r_step;
    }
    if (i != (apot_table.number - 1))
      fprintf(outfile, "\n\n");
  }
#endif /* APOT */

  fclose(outfile);
  printf("Potential plotting data written to \t%s\n", filename);
}

/****************************************************************
 *
 *  write alternate plot version of potential table
 *  (same intervals for all pair and transfer functions)
 *
 ****************************************************************/

void write_altplot_pair(pot_table_t *pt, char *filename)
{
  int   i, j, k, l;
  real  r, rmin = 100., rmax = 0., r_step;
#ifdef EAM
  real  temp;
#endif /* EAM */
  FILE *outfile;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error("Could not open file %s\n", filename);

  /* write data */
  k = 0;
  for (i = 0; i < ntypes; i++) {
    for (j = i; j < ntypes; j++) {
      rmin = MIN(rmin, pt->begin[k]);
      rmax = MAX(rmax, pt->end[k]);
      k++;
    }
    rmin = MIN(rmin, pt->begin[paircol + i]);
    rmax = MAX(rmax, pt->end[paircol + i]);
  }
  k = 0;
  r_step = (rmax - rmin) / (NPLOT - 1);
  for (i = 0; i < ntypes; i++)
    for (j = i; j < ntypes; j++) {
      r = rmin;
      for (l = 0; l < NPLOT - 1; l++) {
#ifdef NEWSCALE
	fprintf(outfile, "%e %e\n", r, (r <= pt->end[k] ? splint_ne(pt,
	      pt->table, k, r) : 0.)
	  + (r <= pt->end[paircol + i] ? splint_ne(pt, pt->table, paircol + i,
	      r) * lambda[j] : 0.)
	  + (r <= pt->end[paircol + j] ? splint_ne(pt, pt->table, paircol + j,
	      r) * lambda[i] : 0.));
#else
	fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, k, r));
#endif /* NEWSCALE */
	r += r_step;
      }
      fprintf(outfile, "%e %e\n\n\n", r, 0.0);
      k++;
    }
#ifdef EAM
  j = k;
  for (i = j; i < j + ntypes; i++) {
    r = rmin;
    for (l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, r <= pt->end[i] ? splint_ne(pt, pt->table,
	  i, r) : 0);
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
  for (i = j + ntypes; i < j + 2 * ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (l = 0; l < NPLOT; l++) {
#ifdef PARABEL
      temp = parab(pt, pt->table, i, r);
#else
      temp = splint_ne(pt, pt->table, i, r);
#endif
#ifdef NEWSCALE
      temp -= lambda[i - (j + ntypes)] * r;
#endif /* NEWSCALE */
      fprintf(outfile, "%e %e\n", r, temp);
      r += r_step;
    }
    fprintf(outfile, "\n\n\n");
  }
#endif
  fclose(outfile);
  printf("Potential plotting data written to %s\n", filename);
}

#ifdef PDIST

/****************************************************************
 *
 * write_pairdist(pot_table_t *pt, char *filename)
 *    - write distribution function of function access
 *
 ****************************************************************/

void write_pairdist(pot_table_t *pt, char *filename)
{
  int  *freq;			/* frequency... */
  int   h, i, j, k, l, typ1, typ2, col;
  real  rr;
  atom_t *atom;
  neigh_t *neigh;
  FILE *outfile;
  char  msg[255];

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error("Could not open file %s\n", filename);

  /* initialize distribution vector */
  freq = (int *)malloc(ndimtot * sizeof(int));
  for (i = 0; i < ndimtot; i++)
    freq[i] = 0;

  for (h = firstconf; h < firstconf + myconf; h++) {
    for (i = 0; i < inconf[h]; i++) {
      atom = atoms + i + cnfstart[h];
      typ1 = atom->typ;

      /* pair potentials */
      for (j = 0; j < atom->n_neigh; j++) {
	neigh = atom->neigh + j;
	typ2 = neigh->typ;
	col = (typ1 <= typ2) ? typ1 * ntypes + typ2 - ((typ1 * (typ1 + 1)) / 2)
	  : typ2 * ntypes + typ1 - ((typ2 * (typ2 + 1)) / 2);
	/* this has already been calculated */
	if (neigh->r < pt->end[col])
	  freq[neigh->slot[0]]++;
#ifdef EAM
	/* transfer function */
	col = paircol + typ2;
	if (neigh->r < pt->end[col])
	  freq[neigh->slot[1]]++;
#endif /* EAM */
      }
#ifdef EAM
      /* embedding function - get index first */
      col = paircol + ntypes + typ1;
      if (format == 3) {
	rr = atom->rho - pt->begin[col];
#ifdef NORESCALE
	if (rr < 0)
	  rr = 0;		/* extrapolation */
	j = MIN((int)(rr * pt->invstep[col]) + pt->first[col], pt->last[col]);
#else
	if (rr < 0)
	  error("short distance");
	j = (int)(rr * pt->invstep[col]) + pt->first[col];
#endif
      } else {			/* format ==4 */
	rr = atom->rho;
	k = pt->first[col];
	l = pt->last[col];
	while (l - k > 1) {
	  j = (k + l) >> 1;
	  if (pt->xcoord[j] > rr)
	    l = j;
	  else
	    k = j;
	}
	j = k;
      }
      freq[j]++;
#endif /* EAM */
    }
  }
  /* finished calculating data - write it to output file */
  j = 0;
  for (col = 0; col < pt->ncols; col++) {
    for (i = pt->first[col]; i < pt->last[col]; i++) {
      rr = 0.5 * (pt->xcoord[i] + pt->xcoord[i + 1]);
      fprintf(outfile, "%f %d\n", rr, freq[i]);
    }
    fprintf(outfile, "\n\n");
  }
  fclose(outfile);
  printf("Distribution data written to\t\t%s\n", filename);
}

#endif /* PDIST */

#ifdef COULOMB

/****************************************************************
 *
 * write endpot for IMD with electrostatics
 *
 ****************************************************************/

void write_coul2imd()
{
  apot_table_t *apt = &apot_table;
  int   i, j;
  FILE *outfile;
  char *filename = "coul2imd";

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error("Could not open file %s\n", filename);

  fprintf(outfile, "charge\t\t");
  for (i = 0; i < ntypes; i++)
    fprintf(outfile, "%f\t", apt->charge[i]);
  fprintf(outfile, "\n");

  if ((strcmp(apt->names[0], "ms") == 0) && (strcmp(apt->names[1], "ms") == 0)
    && (strcmp(apt->names[2], "ms") == 0)) {
    fprintf(outfile, "ms_D\t\t");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][0]);
    fprintf(outfile, "\n");
    fprintf(outfile, "ms_gamma\t");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][1]);
    fprintf(outfile, "\n");
    fprintf(outfile, "ms_r0\t\t");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][2]);
    fprintf(outfile, "\n");
  } else if ((strcmp(apt->names[0], "buck") == 0)
    && (strcmp(apt->names[1], "buck") == 0)
    && (strcmp(apt->names[2], "buck") == 0)) {
    fprintf(outfile, "buck_a\t");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][0]);
    fprintf(outfile, "\n");
    fprintf(outfile, "buck_sigma\t");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][1]);
    fprintf(outfile, "\n");
    fprintf(outfile, "buck_c\t");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][2]);
    fprintf(outfile, "\n");
  } else {
    printf
      ("Only Morse-Stretch or Buckingham potential supported for coul2imd yet.\n");
    printf("No short-range potential printed for coul2imd.");
  }

  fprintf(outfile, "\n");
  fprintf(outfile, "ew_rcut\t\t%f\n", dp_cut);
  fprintf(outfile, "ew_kappa\t%f\n", dp_kappa);
  fprintf(outfile, "r_cut\t\t");
  for (i = 0; i < apt->number; i++)
    fprintf(outfile, "%f\t", apot_table.end[0]);
  fprintf(outfile, "\n");
  fprintf(outfile, "\n");

#ifdef DIPOLE
  fprintf(outfile, "dp_alpha\t");
  for (i = 0; i < ntypes; i++)
    fprintf(outfile, "%f\t", apt->dp_alpha[i]);
  fprintf(outfile, "\n");
  fprintf(outfile, "dp_b\t\t");
  for (i = 0; i < apt->number; i++)
    fprintf(outfile, "%f\t", apt->dp_b[i]);
  fprintf(outfile, "\n");
  fprintf(outfile, "dp_c\t\t");
  for (i = 0; i < apt->number; i++)
    fprintf(outfile, "%f\t", apt->dp_c[i]);
  fprintf(outfile, "\n");
#endif /* DIPOLE */

  fprintf(outfile, "\n");
  fclose(outfile);
}

#endif /* COULOMB */
