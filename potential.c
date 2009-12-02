/****************************************************************
*
*  potential.c: Routines for reading, writing and interpolating a
*      potential table
*
*****************************************************************/
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
/****************************************************************
* $Revision: 1.76 $
* $Date: 2009/12/02 11:02:26 $
*****************************************************************/

#define NPLOT 1000
#define REPULSE
#ifndef POTSCALE
#include "potfit.h"
#else
#include "potscale.h"
#endif

/******************************************************************************
*
* read potential table
*
******************************************************************************/

#ifdef APOT
void read_pot_table(pot_table_t *pt, apot_table_t *apt, char *filename,
		    int ncols)
#else
void read_pot_table(pot_table_t *pt, char *filename, int ncols)
#endif
{
  FILE *infile;
  char  buffer[1024], msg[255], *res, *str;
  int   have_format = 0, end_header = 0;
  int   size, i, j, k, l, *nvals;
  real *val;

  /* open file */
  infile = fopen(filename, "r");
  if (NULL == infile) {
    sprintf(msg, "Could not open file %s\n", filename);
    error(msg);
  }

  /* read the header */
  do {
    /* read one line */
    res = fgets(buffer, 1024, infile);
    if (NULL == res) {
      sprintf(msg, "Unexpected end of file in %s", filename);
      error(msg);
    }
    /* check if it is a header line */
    if (buffer[0] != '#') {
      sprintf(msg, "Header corrupt in file %s", filename);
      error(msg);
    }
    /* stop after last header line */
    if (buffer[1] == 'E') {
      end_header = 1;
    }
    /* invariant potentials */
    else if (buffer[1] == 'I') {
      if (have_format) {
	/* gradient complete */
	for (i = 0; i < size; i++) {
	  str = strtok(((i == 0) ? buffer + 2 : NULL), " \t\r\n");
	  if (str == NULL) {
	    sprintf(msg, "Not enough items in #I header line.");
	    error(msg);
	  } else
	    ((int *)invar_pot)[i] = atoi(str);
	}
	have_invar = 1;
      } else {
	sprintf(msg, "#I needs to be specified after #F in file %s",
		filename);
	error(msg);
      }
    }
#ifndef APOT
    else if (buffer[1] == 'G') {
      if (have_format) {
	/* gradient complete */
	for (i = 0; i < size; i++) {
	  str = strtok(((i == 0) ? buffer + 2 : NULL), " \t\r\n");
	  if (str == NULL) {
	    sprintf(msg, "Not enough items in #G header line.");
	    error(msg);
	  } else
	    ((int *)gradient)[i] = atoi(str);
	}
	have_grad = 1;
      } else {
	sprintf(msg, "#G needs to be specified after #F in file %s",
		filename);
	error(msg);
      }
    }
#endif

    /* see if it is the format line */
    else if (buffer[1] == 'F') {
      /* format complete? */
      if (2 != sscanf((const char *)(buffer + 2), "%d %d", &format, &size)) {
	sprintf(msg, "Corrupt format header line in file %s", filename);
	error(msg);
      }

      /* right number of columns? */
#ifdef EAM
      if (size == ncols + 2 * ntypes) {
	printf("Using EAM potential from file %s\n", filename);
      }
#else
      if (size == ncols) {
	printf("Using pair potential from file %s\n", filename);
      }
#endif
      else {
#ifdef EAM
	sprintf(msg,
		"Wrong number of data columns in file %s,\n should be %d for EAM, but are %d",
		filename, ncols + 2 * ntypes, size);
#else
	sprintf(msg,
		"Wrong number of data columns in file %s,\n should be %d (pair potentials), but are %d",
		filename, ncols, size);
#endif
	error(msg);
      }
      /* recognized format? */
      if ((format != 0) && (format != 3) && (format != 4) && (format != 5)) {
	sprintf(msg, "Unrecognized format specified for file %s", filename);
	error(msg);
      }
      gradient = (int *)malloc(size * sizeof(int));
      invar_pot = (int *)malloc(size * sizeof(int));
      smooth_pot = (int *)malloc(size * sizeof(int));
      reg_for_free(gradient, "gradient");
      reg_for_free(invar_pot, "invar_pot");
      reg_for_free(smooth_pot, "smooth_pot");
      for (i = 0; i < size; i++) {
	gradient[i] = 0;
	invar_pot[i] = 0;
	smooth_pot[i] = 0;
      }
      have_format = 1;
    }
  } while (!end_header);

  /* did we have a format in the header? */
  if (!have_format) {
    sprintf(msg, "Format not specified in header of file %s", filename);
    error(msg);
  } else if (format != 0)
    printf("Potential file format %d.\n", format);
  else
    printf("Potential file format %d (analytical potentials) detected.\n",
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
  reg_for_free(pt->begin, "pt->begin");
  reg_for_free(pt->end, "pt->end");
  reg_for_free(pt->step, "pt->step");
  reg_for_free(pt->invstep, "pt->invstep");
  reg_for_free(pt->first, "pt->first");
  reg_for_free(pt->last, "pt->last");
  if ((pt->begin == NULL) || (pt->end == NULL) || (pt->step == NULL) ||
      (pt->invstep == NULL) || (pt->first == NULL) || (pt->last == NULL) ||
      (nvals == NULL)) {
    sprintf(msg, "Cannot allocate info block for potential table %s",
	    filename);
    error(msg);
  }
#ifdef APOT
#ifdef EAM
  disable_cp = 1;
#endif
  /* allocate memory for analytic potential table */
  apt->number = size;
  apt->total_par = 0;
  apt->n_par = (int *)malloc(size * sizeof(int));
  apt->begin = (real *)malloc(size * sizeof(real));
  apt->end = (real *)malloc(size * sizeof(real));
  apt->param_name = (char ***)malloc(size * sizeof(char **));
  apt->fvalue = (fvalue_pointer *) malloc(size * sizeof(fvalue_pointer));
  apt->co_pol = (real **)malloc(size * sizeof(real *));
  reg_for_free(apt->n_par, "apt->n_par");
  reg_for_free(apt->begin, "apt->begin");
  reg_for_free(apt->end, "apt->end");
  reg_for_free(apt->param_name, "apt->param_name");
  reg_for_free(apt->fvalue, "apt->fvalue");
  reg_for_free(apt->co_pol, "apt->co_pol");
  if (!disable_cp) {
    apt->values = (real **)malloc((size + 1) * sizeof(real *));
    apt->values[size] = (real *)malloc(ntypes * sizeof(real));
    apt->invar_par = (int **)malloc(size * sizeof(int *));
    apt->chempot = apt->values[size];
    reg_for_free(apt->chempot, "apt->chempot");
    apt->pmin = (real **)malloc((size + 1) * sizeof(real *));
    apt->pmin[size] = (real *)malloc(ntypes * sizeof(real));
    reg_for_free(apt->pmin[size], "apt->pmin[size]");
    apt->pmax = (real **)malloc((size + 1) * sizeof(real *));
    apt->pmax[size] = (real *)malloc(ntypes * sizeof(real));
    reg_for_free(apt->pmax[size], "apt->pmax[size}");
  } else {
    apt->values = (real **)malloc(size * sizeof(real *));
    apt->invar_par = (int **)malloc(size * sizeof(int *));
    apt->pmin = (real **)malloc(size * sizeof(real *));
    apt->pmax = (real **)malloc(size * sizeof(real *));
  }
  reg_for_free(apt->values, "apt->values");
  reg_for_free(apt->invar_par, "apt->invar_par");
  reg_for_free(apt->pmin, "apt->pmin");
  reg_for_free(apt->pmax, "apt->pmax");
  apt->names = (char **)malloc(size * sizeof(char *));
  pot_list_length = (int *)malloc(size * sizeof(int));
  pot_list = (int ***)malloc(size * sizeof(int **));
  reg_for_free(apt->names, "apt->names");
  reg_for_free(pot_list_length, "pot_list_length");
  reg_for_free(pot_list, "pot_list");
  for (i = 0; i < size; i++) {
    apt->names[i] = (char *)malloc(20 * sizeof(char));
    reg_for_free(apt->names[i], "apt->names[i]");
    pot_list[i] = (int **)malloc(sizeof(int *));
    reg_for_free(pot_list[i], "pot_list[i]");
    pot_list_length[i] = 0;
  }
  if ((apt->n_par == NULL) || (apt->begin == NULL) || (apt->end == NULL)
      || (apt->fvalue == NULL) || (apt->names == NULL) || (apt->pmin == NULL)
      || (apt->pmax == NULL) || (apt->param_name == NULL)
      || (apt->values == NULL) || (pot_list_length == NULL)
      || (pot_list == NULL)) {
    sprintf(msg, "Cannot allocate info block for analytic potential table %s",
	    filename);
    error(msg);
  }
#endif
  switch (format) {
#ifdef APOT
      case 0:
	read_apot_table(pt, apt, filename, infile);
	break;
#endif
      case 3:
	read_pot_table3(pt, size, ncols, nvals, filename, infile);
	break;
      case 4:
	read_pot_table4(pt, size, ncols, nvals, filename, infile);
	break;
      case 5:
	read_pot_table5(pt, size, ncols, nvals, filename, infile);
  }
  fclose(infile);


#ifndef POTSCALE		/* not needed in potscale program */
  /* compute rcut and rmin */
  rcut = (real *)malloc(ntypes * ntypes * sizeof(real));
  reg_for_free(rcut, "rcut");
  if (NULL == rcut)
    error("Cannot allocate rcut");
  rmin = (real *)malloc(ntypes * ntypes * sizeof(real));
  reg_for_free(rmin, "rmin");
  if (NULL == rmin)
    error("Cannot allocate rmin");
#ifdef APOT
#ifdef EAM
  pot_index =
    (int *)malloc(((ntypes * (ntypes + 1) / 2) + ntypes) * sizeof(int));
  reg_for_free(pot_index, "pot_index");
  if (NULL == pot_index)
    error("Cannot allocate pot_index");
  for (i = 0; i < (ntypes * (ntypes + 1) / 2 + ntypes); i++)
    pot_index[i] = ntypes * ntypes;
#else /* EAM */
  pot_index = (int *)malloc(ntypes * (ntypes + 1) / 2 * sizeof(int));
  reg_for_free(pot_index, "pot_index");
  if (NULL == pot_index)
    error("Cannot allocate pot_index");
  for (i = 0; i < ntypes * (ntypes + 1) / 2; i++)
    pot_index[i] = ntypes * ntypes;
#endif /* EAM */
#endif /* APOT */
  for (i = 0; i < ntypes; i++)
    for (j = 0; j < ntypes; j++) {
      k = (i <= j) ? i * ntypes + j - ((i * (i + 1)) / 2)
	: j * ntypes + i - ((j * (j + 1)) / 2);
      rmin[i * ntypes + j] = pt->begin[k];
#ifdef APOT
      pot_index[k] = MIN(pot_index[k], i * ntypes + j);
#endif /* APOT */
      rcut[i * ntypes + j] = pt->end[k];
    }
#if defined EAM && defined APOT
  j = 0;
  for (i = 0; i < ntypes; i++) {
    pot_index[++k] = pot_index[j];
    j += ntypes - i;
  }
#endif
#ifdef EAM
  for (i = 0; i < ntypes; i++) {
    for (j = 0; j < ntypes; j++) {
      rcut[i * ntypes + j] = MAX(rcut[i * ntypes + j],
				 pt->end[(ntypes * (ntypes + 1)) / 2 + i]);
      rcut[i * ntypes + j] = MAX(rcut[i * ntypes + j],
				 pt->end[(ntypes * (ntypes + 1)) / 2 + j]);
      rmin[i * ntypes + j] = MIN(rmin[i * ntypes + j],
				 pt->begin[(ntypes * (ntypes + 1)) / 2 + i]);
      rmin[i * ntypes + j] = MIN(rmin[i * ntypes + j],
				 pt->begin[(ntypes * (ntypes + 1)) / 2 + j]);
    }
  }
#endif
#endif /* POTSCALE */
  paircol = (ntypes * (ntypes + 1)) / 2;
#ifndef APOT
  /* read maximal changes file */
  maxchange = (real *)malloc(pt->len * sizeof(real));
  reg_for_free(maxchange, "maxchange");
  if (usemaxch) {
    /* open file */
    infile = fopen(maxchfile, "r");
    if (NULL == infile) {
      sprintf(msg, "Could not open file %s\n", maxchfile);
      error(msg);
    }
    val = maxchange;
    for (i = 0; i < pt->len; i++) {
      if (1 > fscanf(infile, " %lf\n", val)) {
	sprintf(msg, "Premature end of maxch file %s", maxchfile);
	error(msg);
      } else
	val++;
    }
    fclose(infile);
  }
#endif
  for (i = 0; i < ntypes; i++) {
    for (j = 0; j < ntypes; j++) {
      rcutmax = MAX(rcutmax, rcut[i + ntypes * j]);
    }
  }

  free(nvals);
  return;
}

/*****************************************************************************
*
*  read potential in analytic format:
*  	for more information an how to specify an analytic potential
*  	please check the documentation
*
******************************************************************************/

#ifdef APOT
void read_apot_table(pot_table_t *pt, apot_table_t *apt, char *filename,
		     FILE *infile)
{
  int   c, i, j, k, l, comment;
  char  msg[255];
  char  buffer[255];
  char  name[255];
  char *token;
  real *val, *list, temp;
  fpos_t filepos;

#ifndef EAM
  /* read cp */
  fgetpos(infile, &filepos);
  if (!disable_cp) {
    for (i = 0; i < ntypes; i++) {
      if (4 > fscanf(infile, "%s %lf %lf %lf", buffer, &apt->chempot[i],
		     &apt->pmin[apt->number][i],
		     &apt->pmax[apt->number][i])) {
	sprintf(msg, "Could not read chemical potential for atomtype #%d\n",
		i);
	error(msg);
      }

      /* split cp and _# */
      token = strchr(buffer, '_');
      if (token != NULL) {
	strncpy(name, buffer, strlen(buffer) - strlen(token));
	name[strlen(buffer) - strlen(token)] = '\0';
      }
      if (strcmp("cp", name) != 0) {
	fprintf(stderr, "Found \"%s\" instead of \"cp\"\n", name);
	sprintf(msg, "No chemical potentials found in %s.\n", filename);
	error(msg);
      }
    }
    fgetpos(infile, &filepos);
    if (2 > fscanf(infile, "%s %d", buffer, &compnodes)) {
      if (strcmp("type", buffer) == 0)
	compnodes = -1;
      else {
	sprintf(msg,
		"Could not read number of composition nodes from potential file.\n");
	error(msg);
      }
    }
    if (strcmp(buffer, "cn") != 0 && ntypes > 1 && compnodes != -1) {
      sprintf(msg,
	      "No composition nodes found in %s.\n(cn 0 is also allowed)\n",
	      filename);
      error(msg);
    }
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
	if (4 >
	    fscanf(infile, "%lf %lf %lf %lf", &compnodelist[j],
		   &apt->chempot[ntypes + j],
		   &apt->pmin[apt->number][ntypes + j],
		   &apt->pmax[apt->number][ntypes + j])) {
	  sprintf(msg, "Could not read composition node %d\n", j + 1);
	  error(msg);
	}
	if (apt->pmin[apt->number][ntypes + j] > apt->chempot[ntypes + j]
	    || apt->pmax[apt->number][ntypes + j] < apt->chempot[ntypes + j]) {
	  sprintf(msg, "composition node %d is out of bounds.\n", j + 1);
	  error(msg);
	}
      }

      /* check compnodes for valid values */
      if (ntypes == 2) {
	for (j = 0; j < compnodes; j++)
	  if (compnodelist[j] > 1 || compnodelist[j] < 0) {
	    sprintf(msg,
		    "Composition node %d is %f but should be inside [0,1].\n",
		    j + 1, compnodelist[j]);
	    error(msg);
	  }
      }
    }
    if (compnodes != -1)
      printf
	("Enabled chemical potentials with %d extra composition node(s).\n",
	 compnodes);
    if (compnodes == -1)
      compnodes = 0;
  }
  fsetpos(infile, &filepos);
#endif

  do {
    fgetpos(infile, &filepos);
    fscanf(infile, "%s", buffer);
  } while (strcmp(buffer, "type") != 0);
  fsetpos(infile, &filepos);

  for (i = 0; i < apt->number; i++) {
    /* read type */
    if (2 > fscanf(infile, "%s %s", buffer, name)) {
      sprintf(msg, "Premature end of potential file %s", filename);
      error(msg);
    }
    if (strcmp(buffer, "type") != 0) {
      sprintf(msg,
	      "Unknown keyword in file %s, expected \"type\" but found \"%s\".",
	      filename, buffer);
      error(msg);
    }

    /* split name and _sc */
    token = strrchr(name, '_');
    if (token != NULL && strcmp(token + 1, "sc") == 0) {
      strncpy(buffer, name, strlen(name) - 3);
      buffer[strlen(name) - 3] = '\0';
      strcpy(name, buffer);
      smooth_pot[i] = 1;
      do_smooth = 1;
    }

    if (apot_parameters(name) == -1) {
      sprintf(msg,
	      "Unknown function type in file %s, please define \"%s\" in functions.c.",
	      filename, name);
      error(msg);
    }

    strcpy(apt->names[i], name);
    apt->n_par[i] = apot_parameters(name);
    apt->total_par += apt->n_par[i];

    /* read cutoff */
    if (i < (ntypes * (ntypes + 1) / 2 + ntypes)) {
      if (2 > fscanf(infile, "%s %lf", buffer, &apt->end[i])) {
	printf("%s\n", buffer);
	sprintf(msg,
		"Could not read cutoff for potential #%d in file %s\nAborting",
		i, filename);
	error(msg);
      }
      if (strcmp(buffer, "cutoff") != 0) {
	sprintf(msg,
		"No cutoff for potential #%d found after \"type\" keyword in file %s\nAborting",
		i, filename);
	error(msg);
      }
    } else {
      fgetpos(infile, &filepos);
      fscanf(infile, "%s", buffer);
      if (strncmp(buffer, "cutoff", 6) != 0)
	fsetpos(infile, &filepos);
#ifdef DEBUG
      else
	fprintf(stderr, "Ignoring cutoff for embedding function %d\n", i);
#endif
      apt->end[i] = 2;
    }

    /* set small begin to prevent division by zero-errors */
    apt->begin[i] = 0.0001;
    apt->values[i] = (real *)malloc(apt->n_par[i] * sizeof(real));
    reg_for_free(apt->values[i], "apt->values[i]");
    apt->invar_par[i] = (int *)malloc((apt->n_par[i] + 1) * sizeof(int));
    reg_for_free(apt->invar_par[i], "apt->invar_par[i]");
    apt->pmin[i] = (real *)malloc(apt->n_par[i] * sizeof(real));
    reg_for_free(apt->pmin[i], "apt->pmin[i]");
    apt->pmax[i] = (real *)malloc(apt->n_par[i] * sizeof(real));
    reg_for_free(apt->pmax[i], "apt->pmax[i]");
    apt->param_name[i] = (char **)malloc(apt->n_par[i] * sizeof(char *));
    reg_for_free(apt->param_name[i], "apt->param_name[i]");
    if (smooth_pot[i] == 1)
      apt->co_pol[i] = (real *)malloc(5 * sizeof(real));
    reg_for_free(apt->co_pol[i], "apt->co_pol[i]");

    if (NULL == apt->values[i] || NULL == apt->pmin[i]
	|| NULL == apt->pmax[i] || NULL == apt->param_name[i]) {
      sprintf(msg,
	      "Cannot allocate memory for potential paramters.\nAborting");
      error(msg);
    }

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
      /* Warning: The array for parameter names holds only 30 chars
       * so there is the chance of a buffer overflow.
       * If you need longer parameter names increase the size of
       * the array.
       */
      apt->param_name[i][j] = (char *)malloc(30 * sizeof(char));
      reg_for_free(apt->param_name[i][j], "apt->param_name[i][j]");
      if (NULL == apt->param_name[i][j])
	error("Error in allocating memory for parameter name");
      if (4 > fscanf(infile, "%s %lf %lf %lf", apt->param_name[i][j],
		     &apt->values[i][j], &apt->pmin[i][j],
		     &apt->pmax[i][j])) {
	if (strcmp(apt->param_name[i][j], "type") == 0) {
	  sprintf(msg,
		  "Not enough parameters for potential #%d in file %s specified.\nYou specified %d parameters, but needed are %d.\nAborting",
		  i + 1, filename, j, apt->n_par[i]);
	  error(msg);
	}
	sprintf(msg,
		"Could not read parameter #%d of potential #%d in file %s",
		j + 1, i + 1, filename);
	error(msg);
	printf("%s\n", apt->param_name[i][0]);
      }
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
	if (!opt) {
	  if (apt->values[i][j] < apt->pmin[i][j]) {
	    fprintf(stderr, "\n --> Warning <--\n");
	    fprintf(stderr,
		    "Starting value for paramter #%d in potential #%d is smaller than the specified minimum.\nResetting it from %f to pmin (%f)\n",
		    j + 1, i + 1, apt->values[i][j], apt->pmin[i][j]);
	    apt->values[i][j] = apt->pmin[i][j];
	  }
	  if (apt->values[i][j] > apt->pmax[i][j]) {
	    fprintf(stderr, "\n --> Warning <--\n");
	    fprintf(stderr,
		    "Starting value for paramter #%d in potential #%d is bigger than the specified maximumm.\nResetting it from %f to pmax (%f)\n",
		    j + 1, i + 1, apt->values[i][j], apt->pmax[i][j]);
	    apt->values[i][j] = apt->pmax[i][j];
	  }
	} else {
	  apt->values[i][j] = (apt->pmin[i][j] + apt->pmax[i][j]) / 2.;
	  fprintf(stderr, "\n --> Warning <--\n");
	  if (apt->values[i][j] == 0) {
	    apt->values[i][j] = apt->pmax[i][j] * 0.25;
	    fprintf(stderr,
		    "Starting value for paramter #%d in potential #%d is outside of specified adjustment range.\nAutosetting it to %f (pmax*.25)\n",
		    j + 1, i + 1, apt->values[i][j]);
	  } else {
	    fprintf(stderr,
		    "Starting value for paramter #%d in potential #%d is outside of specified adjustment range.\nAutosetting it to %f ((pmin+pmax)/2)\n",
		    j + 1, i + 1, apt->values[i][j]);
	  }
	}
      }
    }
  }
  printf("\n");

  /* assign the potential functions to the function pointers */
  if (apot_assign_functions(apt) == -1) {
    sprintf(msg, "Could not assign the function pointers.\nAborting");
    error(msg);
  }

  if (!disable_cp) {
    cp_start = apt->total_par + ntypes * (ntypes + 1);
    apt->total_par += (ntypes + compnodes);
  }

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
  if (!disable_cp) {
    pt->len += (ntypes + compnodes);
  }

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
      || (apt->idxparam == NULL)) {
    sprintf(msg, "Cannot allocate memory for potential table.\nAborting");
    error(msg);
  }
  for (i = 0; i < pt->len; i++) {
    pt->table[i] = 0;
    calc_list[i] = 0;
    pt->idx[i] = 0;
  }
  k = 0;
  l = 0;
  val = pt->table;
  list = calc_list;
  for (i = 0; i < apt->number; i++) {
    val += 2;
    list += 2;
    l += 2;
    for (j = 0; j < apt->n_par[i]; j++) {
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
    if (!invar_pot[i]) {
      pt->idxlen += apt->n_par[i] - apt->invar_par[i][apt->n_par[i]];
    }
  }
#ifndef EAM
  if (!disable_cp) {
    init_chemical_potential(ntypes);
    i = apt->number;
    for (j = 0; j < (ntypes + compnodes); j++) {
      *(val + j) = apt->values[i][j];
      pt->idx[k] = l++;
      apt->idxpot[k] = i;
      apt->idxparam[k++] = j;
    }
    pt->idxlen += (ntypes + compnodes);
  }
#endif

  init_calc_table(pt, &calc_pot);
  return;
}
#endif

/*****************************************************************************
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
******************************************************************************/

void read_pot_table3(pot_table_t *pt, int size, int ncols, int *nvals,
		     char *filename, FILE *infile)
{
  int   i, j, k, l;
  char  msg[255];
  real *val;

  /* read the info block of the function table */
  for (i = 0; i < size; i++) {
    if (3 >
	fscanf(infile, "%lf %lf %d", &pt->begin[i], &pt->end[i], &nvals[i])) {
      sprintf(msg, "Premature end of potential file %s", filename);
      error(msg);
    }
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
  if ((NULL == pt->table) || (NULL == pt->idx) || (NULL == pt->d2tab)) {
    error("Cannot allocate memory for potential table");
  }

  /* input loop */
  val = pt->table;
  k = 0;
  l = 0;
  for (i = 0; i < ncols; i++) {	/* read in pair pot */
    if (have_grad) {		/* read gradient */
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      }
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
    for (j = 0; j < nvals[i]; j++) {	/* read values */
      if (1 > fscanf(infile, "%lf\n", val)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      } else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;
    }
  }
#ifdef EAM
  for (i = ncols; i < ncols + ntypes; i++) {	/* read in rho */
    if (have_grad) {		/* read gradient */
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      }
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
    for (j = 0; j < nvals[i]; j++) {	/* read values */
      if (1 > fscanf(infile, "%lf\n", val)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      } else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if ((!invar_pot[i]) && (j < nvals[i] - 1))
	pt->idx[k++] = l++;
      else
	l++;
    }
  }
  for (i = ncols + ntypes; i < ncols + 2 * ntypes; i++) {	/* read in F */
    if (have_grad) {		/* read gradient */
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      }
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
    for (j = 0; j < nvals[i]; j++) {	/* read values */
      if (1 > fscanf(infile, "%lf\n", val)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      } else
	val++;
      pt->xcoord[l] = pt->begin[i] + j * pt->step[i];
      if (!invar_pot[i])
	pt->idx[k++] = l++;
      else
	l++;
    }
  }

#endif
  pt->idxlen = k;
  init_calc_table(pt, &calc_pot);

}

/*****************************************************************************
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
******************************************************************************/

void read_pot_table4(pot_table_t *pt, int size, int ncols, int *nvals,
		     char *filename, FILE *infile)
{
  int   i, k, l, j;
  char  msg[255];
  real *val, *ord;

  /* read the info block of the function table */
  for (i = 0; i < size; i++) {
    if (1 > fscanf(infile, "%d", &nvals[i])) {
      sprintf(msg, "Premature end of potential file %s", filename);
      error(msg);
    }
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
  if ((NULL == pt->table) || (NULL == pt->idx) || (NULL == pt->d2tab)) {
    error("Cannot allocate memory for potential table");
  }

  /* input loop */
  val = pt->table;
  ord = pt->xcoord;
  k = 0;
  l = 0;
  for (i = 0; i < ncols; i++) {	/* read in pair pot */
    if (have_grad) {		/* read gradient */
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      }
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
    for (j = 0; j < nvals[i]; j++) {	/* read values */
      if (2 > fscanf(infile, "%lf %lf\n", ord, val)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      } else {
	val++;
	ord++;
      }
      if ((j > 0) && (*(ord - 1) <= *(ord - 2))) {
	sprintf(msg, "Abscissa not monotonous in potential %d.", i);
	error(msg);
      }
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
#ifdef EAM
  for (i = ncols; i < ncols + ntypes; i++) {	/* read in rho */
    if (have_grad) {		/* read gradient */
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      }
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
    for (j = 0; j < nvals[i]; j++) {	/* read values */
      if (2 > fscanf(infile, "%lf %lf\n", ord, val)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      } else {
	ord++;
	val++;
      }
      if ((j > 0) && (*(ord - 1) <= *(ord - 2))) {
	sprintf(msg, "Abscissa not monotonous in potential %d.", i);
	error(msg);
      }
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
  for (i = ncols + ntypes; i < ncols + 2 * ntypes; i++) {	/* read in F */
    if (have_grad) {		/* read gradient */
      if (2 > fscanf(infile, "%lf %lf\n", val, val + 1)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      }
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
    for (j = 0; j < nvals[i]; j++) {	/* read values */
      if (1 > fscanf(infile, "%lf %lf\n", ord, val)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      } else {
	ord++;
	val++;
      }
      if ((j > 0) && (*(ord - 1) <= *(ord - 2))) {
	sprintf(msg, "Abscissa not monotonous in potential %d.", i);
	error(msg);
      }
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

#endif
  pt->idxlen = k;
  init_calc_table(pt, &calc_pot);
}

/*****************************************************************************
*
*  read potential in fifth format:
*
*  Functions specified by gaussians.
*
*  THIS IS BY NO MEANS FINISHED! IT WILL NOT WORK!
*
*  Header:  one line for each function with
*           rbegin rstart npoints
*
*  Table: Center, width, amplitude of Gaussians,
*         functions separated by blank lines
*
*
*  d2tab is used to store the width of the Gaussians.
*
******************************************************************************/

void read_pot_table5(pot_table_t *pt, int size, int ncols, int *nvals,
		     char *filename, FILE *infile)
{
  int   i, j, k, l;
  char  msg[255];
  real *val, *ord, *width;

  /* THIS FUNCTIONALITY IS NOT YET IMPLEMENTED */

  error("Potential format 5 is not yet implemented");

  /* read the info block of the function table */
  for (i = 0; i < size; i++) {
    if (3 >
	fscanf(infile, "%lf %lf %d", &pt->begin[i], &pt->end[i], &nvals[i])) {
      sprintf(msg, "Premature end of potential file %s", filename);
      error(msg);
    }
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
  pt->xcoord = (real *)malloc(pt->len * sizeof(real));
  pt->d2tab = (real *)malloc(pt->len * sizeof(real));
  pt->idx = (int *)malloc(pt->len * sizeof(int));
  if ((NULL == pt->table) || (NULL == pt->idx) || (NULL == pt->d2tab)) {
    error("Cannot allocate memory for potential table");
  }

  /* input loop */
  val = pt->table;
  ord = pt->xcoord;
  width = pt->d2tab;
  k = 0;
  l = 0;
  for (i = 0; i < ncols; i++) {	/* read in pair pot */
/*     if (have_grad) { 		/\* read gradient *\/ */
/*       if (2>fscanf(infile, "%lf %lf\n", val, val+1)){ */
/*         sprintf(msg, "Premature end of potential file %s", filename); */
/* 	error(msg); */
/*       }  */
/*     } else { */
    *val = 1e30;
    *(val + 1) = 0.;
/*     }  */
    val += 2;
    ord += 2;
    width += 2;
/*     if (gradient[i]>>1) pt->idx[k++]=l++; else l++; */
/*     if (gradient[i]% 2) pt->idx[k++]=l++; else l++; */
    l += 2;
    for (j = 0; j < nvals[i]; j++) {	/* read values */
      if (3 > fscanf(infile, "%lf %lf %lf\n", ord, width, val)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      } else {
	val++;
	ord++;
	width++;
      }
//      pt->xcoord[l]=pt->begin[i] + j * pt->step[i] ;
//      if (j<nvals[i]-1) pt->idx[k++] = l++;
      if (!invar_pot[i])
	pt->idx[k++] = l++;
      else
	l++;
//    else l++;
    }
  }
#ifdef EAM
  for (i = ncols; i < ncols + ntypes; i++) {	/* read in rho */
/*     if (have_grad) { 		/\* read gradient *\/ */
/*       if (2>fscanf(infile, "%lf %lf\n", val, val+1)){ */
/*         sprintf(msg, "Premature end of potential file %s", filename); */
/* 	error(msg); */
/*       }  */
/*     } else { */
    *val = 1e30;
    *(val + 1) = 0.;
/*     }  */
    val += 2;
    ord += 2;
    width += 2;
/*     if (gradient[i]>>1) pt->idx[k++]=l++; else l++; */
/*     if (gradient[i]% 2) pt->idx[k++]=l++; else l++; */
    l += 2;
    for (j = 0; j < nvals[i]; j++) {	/* read values */
      if (3 > fscanf(infile, "%lf %lf %lf\n", ord, width, val)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      } else {
	ord++;
	val++;
	width++;
      }
/*       if ((j>0) && (*(ord-1) <= *(ord-2))) { */
/* 	sprintf(msg, "Abscissa not monotonous in potential %d.",i); */
/* 	error(msg); */
/*       } */
/*       if (j<nvals[i]-1) pt->idx[k++] = l++; */
      if (!invar_pot[i])
	pt->idx[k++] = l++;
      else
	l++;
/*       else l++; */
    }
/*     pt->begin[i]=pt->xcoord[pt->first[i]]; */
/*     pt->end[i]  =pt->xcoord[pt->last[i]]; */
/*     /\* pt->step is average step length.. *\/ */
/*     pt->step[i] = (pt->end[i]-pt->begin[i])/((real) nvals[i]-1); */
/*     pt->invstep[i] = 1./pt->step[i];  */

  }
  for (i = ncols + ntypes; i < ncols + 2 * ntypes; i++) {	/* read in F */
/*     if (have_grad) { 		/\* read gradient *\/ */
/*       if (2>fscanf(infile, "%lf %lf\n", val, val+1)){ */
/*         sprintf(msg, "Premature end of potential file %s", filename); */
/* 	error(msg); */
/*       }  */
/*     } else { */
    *val = 1e30;
    *(val + 1) = 1.e30;
/*     }  */
    val += 2;
    ord += 2;
    width += 2;
/*     if (gradient[i]>>1) pt->idx[k++]=l++; else l++; */
/*     if (gradient[i]% 2) pt->idx[k++]=l++; else l++; */
    l += 2;
    for (j = 0; j < nvals[i]; j++) {	/* read values */
      if (3 > fscanf(infile, "%lf %lf %lf\n", ord, width, val)) {
	sprintf(msg, "Premature end of potential file %s", filename);
	error(msg);
      } else {
	ord++;
	val++;
	width++;
      }
/*       if ((j>0) && (*(ord-1) <= *(ord-2))) { */
/* 	sprintf(msg, "Abscissa not monotonous in potential %d.",i); */
/* 	error(msg); */
/*       } */
      if (!invar_pot[i])
	pt->idx[k++] = l++;
      else
	l++;
    }
/*     pt->begin[i]=pt->xcoord[pt->first[i]]; */
/*     pt->end[i]  =pt->xcoord[pt->last[i]]; */
/*     /\* pt->step is average step length.. *\/ */
/*     pt->step[i] = (pt->end[i]-pt->begin[i])/((real) nvals[i]-1); */
/*     pt->invstep[i] = 1./pt->step[i];  */
  }

#endif
  pt->idxlen = k;
  init_calc_table(pt, &calc_pot);
  update_calc_table(pt->table, calc_pot.table, 1);
}

/*****************************************************************************
*
*  init_calc_table: Initialize table used for calculation.
*
*  *  Header:  one line for each function with
*           rbegin rstart npoints
*
*  Table: Center, width, amplitude of Gaussians,
*         functions separated by blank lines
*
******************************************************************************/

void init_calc_table(pot_table_t *optt, pot_table_t *calct)
{
  int   i, j, size, x = 0, index;
  int  *sp;
  real *val, *ord, f;

  switch (format) {
#ifdef APOT
      case 0:
	{
	  size = apot_table.number;
	  calct->len =
	    size * APOT_STEPS + 2 * optt->ncols + ntypes + compnodes;
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
	  if (calct->first == NULL || calct->last == NULL
	      || calct->step == NULL || calct->invstep == NULL
	      || calct->xcoord == NULL || calct->table == NULL
	      || calct->d2tab == NULL || calct->idx == NULL)
	    error("Cannot allocate info block for calc potential table\n");
	  for (i = 0; i < size; i++) {
	    val = apot_table.values[i];
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
	      calct->table[index] = f;
	      calct->idx[i * APOT_STEPS + j] = index;
	    }
	  }
	}
	break;
#endif
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
	break;
      case 5:			/* Here's work to do */
	size = optt->ncols;
	sp = (int *)malloc(size * sizeof(int));
	/* allocate info block of function table */
	calct->ncols = size;
	calct->len = 0;
	calct->begin = (real *)malloc(size * sizeof(real));
	calct->end = (real *)malloc(size * sizeof(real));
	calct->step = (real *)malloc(size * sizeof(real));
	calct->invstep = (real *)malloc(size * sizeof(real));
	calct->first = (int *)malloc(size * sizeof(int));
	calct->last = (int *)malloc(size * sizeof(int));
	if ((calct->begin == NULL) || (calct->end == NULL) ||
	    (calct->step == NULL) || (calct->invstep == NULL) ||
	    (calct->first == NULL) || (calct->last == NULL)) {
	  error("Cannot allocate info block for calc potential table\n");
	}
	for (i = 0; i < size; i++) {
	  calct->begin[i] = optt->begin[i];
	  calct->end[i] = optt->end[i];
	  sp[i] = 10 * (optt->last[i] - optt->first[i] + 1);	/* oversample * 10 */
	  calct->step[i] = (calct->end[i] - calct->begin[i]) / (sp[i] - 1);
	  calct->invstep[i] = 1. / calct->step[i];
	  if (i == 0)
	    calct->first[i] = 2;
	  else
	    calct->first[i] = calct->last[i - 1] + 3;
	  calct->last[i] = calct->first[i] + sp[i] - 1;
	  calct->len = calct->last[i] + 1;
	}
	calct->table = (real *)malloc(calct->len * sizeof(real));
	calct->xcoord = (real *)malloc(calct->len * sizeof(real));
	calct->d2tab = (real *)malloc(calct->len * sizeof(real));
	calct->idx = (int *)malloc(calct->len * sizeof(int));
	if ((NULL == calct->table) || (NULL == calct->idx) ||
	    (NULL == calct->d2tab)) {
	  error("Cannot allocate memory for potential table");
	}
  }
}

void update_calc_table(real *xi_opt, real *xi_calc, int do_all)
{
  int   i, j, k, l, size, change;
  int  *sp;
  real  r, x0, temp;
  real *val, *ord, *cut;
  static real *params, *list;

  switch (format) {
#ifdef APOT
      case 0:
	{
	  cut = apot_table.end;
	  if (do_smooth && params == NULL) {
	    params = (real *)malloc(4 * sizeof(real));
	    reg_for_free(params, "params");
	  }
	  val = xi_opt;
	  list = calc_list + 2;
	  for (i = 0; i < calc_pot.ncols; i++) {
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
	    if (change || do_all) {
	      if ((!invar_pot[i] && !smooth_pot[i])
		  || i >= (ntypes * (ntypes + 1) / 2 + ntypes)) {
		for (j = 0; j < APOT_STEPS; j++) {
		  k = i * APOT_STEPS + (i + 1) * 2 + j;
		  apot_table.fvalue[i] (calc_pot.xcoord[k], val, xi_calc + k);
		}
	      } else if (!invar_pot[i] && smooth_pot[i]) {
		k = i * APOT_STEPS + (i + 1) * 2;
		x0 = smooth(apot_table.fvalue[i], cut[i], val, params);
		if (myid == 0) {
		  apot_table.co_pol[i][0] = x0;
		  for (j = 0; j < 4; j++)
		    apot_table.co_pol[i][j + 1] = params[j];
		}
		for (j = 0; j < APOT_STEPS; j++) {
		  temp = j * calc_pot.step[i] + apot_table.begin[i];
		  if (temp <= x0)
		    apot_table.fvalue[i] (temp, val, xi_calc + k + j);
		  else
		    *(xi_calc + k + j) =
		      params[0] * temp * temp * temp +
		      temp * temp * params[1] + temp * params[2] + params[3];
		  calc_pot.xcoord[k + j] =
		    calc_pot.begin[i] + j * calc_pot.step[i];
		}
	      }
	    }
	    val += apot_table.n_par[i];
	    list += apot_table.n_par[i] + 2;
	  }
	}
	return;
#endif
      case 3:			/* fall through */
      case 4:{
	  /* Nothing to do */
	  return;
	}
      case 5:{
	  error("This potential is not yet implemented");
	  size = calc_pot.ncols;
	  /* calculation loop */
	  val = xi_calc;
	  ord = calc_pot.xcoord;
	  k = 0;
	  l = 0;
	  for (i = 0; i < size; i++) {	/* Functions */
	    /* Set gradients */
	    *val = 0;
	    *(val + 1) = 0;
	    for (j = opt_pot.first[i]; j <= opt_pot.last[i]; j++)
	      *val +=
		2 * opt_pot.table[j] * (calc_pot.begin[i] -
					opt_pot.xcoord[j]) /
		SQRREAL(opt_pot.d2tab[j]) *
		exp(SQR
		    ((calc_pot.begin[i] -
		      opt_pot.xcoord[j]) / opt_pot.d2tab[j]));
#ifdef EAM
	    if (i >= size - ntypes)	/* Embedding function */
	      for (j = opt_pot.first[i]; j <= opt_pot.last[i]; j++)
		*val +=
		  2 * opt_pot.table[j] * (calc_pot.end[i] -
					  opt_pot.xcoord[j]) /
		  SQRREAL(opt_pot.d2tab[j]) *
		  exp(SQRREAL
		      ((calc_pot.end[i] -
			opt_pot.xcoord[j]) / opt_pot.d2tab[j]));

#endif /* EAM */
	    val += 2;
	    ord += 2;
	    /* set values */
	    r = calc_pot.begin[i];
	    for (k = calc_pot.first[i]; k < calc_pot.last[i]; k++) {
	      *val = 0;
	      *ord = r;
	      for (j = opt_pot.first[i]; j <= opt_pot.last[i]; j++)
		*val += opt_pot.table[j] *
		  exp(SQRREAL((r - opt_pot.xcoord[j]) / opt_pot.d2tab[j]));
	      val++;
	      ord++;
	      r += calc_pot.step[i];
	    }
	    *ord = r;
	    *val = 0;		/* Cut off at r=r_cut */
#ifdef EAM
	    if (i >= size - ntypes)	/* Embedding function not cut off */
	      for (j = opt_pot.first[i]; j <= opt_pot.last[i]; j++)
		*val += opt_pot.table[j] *
		  exp(SQRREAL((r - opt_pot.xcoord[j]) / opt_pot.d2tab[j]));
#endif /* EAM */
	  }
	}
	return;
  }
}


#ifndef POTSCALE
#ifdef OLDCODE
/*****************************************************************************
*
*  Evaluate derivative of potential with quadratic interpolation.
*  col is typ1 * ntypes + typ2.
*
******************************************************************************/

real grad2(pot_table_t *pt, real *xi, int col, real r)
{
  real  rr, istep, chi, p0, p1, p2, dv, d2v;
  int   k;

  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0)
    error("short distance!");

  /* indices into potential table */
  istep = pt->invstep[col];
  k = (int)(rr * istep);
  chi = (rr - k * pt->step[col]) * istep;
  k += pt->first[col];

  /* intermediate values */
  p0 = (k <= pt->last[col]) ? xi[k++] : 0.0;
  p1 = (k <= pt->last[col]) ? xi[k++] : 0.0;
  p2 = (k <= pt->last[col]) ? xi[k++] : 0.0;
  dv = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* return the derivative */
  return istep * (dv + (chi - 0.5) * d2v);
}

/*****************************************************************************
*
*  Evaluate derivative of potential with cubic interpolation.
*  col is typ1 * ntypes + typ2.
*
******************************************************************************/

real grad3(pot_table_t *pt, real *xi, int col, real r)
{
  real  rr, istep, chi, p0, p1, p2, p3;
  real  dfac0, dfac1, dfac2, dfac3;
  int   k;

  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0)
    error("short distance!");

  /* indices into potential table */
  istep = pt->invstep[col];
  k = (int)(rr * istep);
  if (k == 0)
    return grad2(pt, xi, col, r);	/* parabolic fit if on left border */
  chi = (rr - k * pt->step[col]) * istep;
  k += pt->first[col];
  k--;

  /* intermediate values */
  if (k <= pt->last[col])
    p0 = xi[k++];
  else
    return 0.0;
  if (k <= pt->last[col])
    p1 = xi[k++];
  else
    return 0.0;
  if (k <= pt->last[col])
    p2 = xi[k++];
  else
    return 0.0;
  if (k <= pt->last[col])
    p3 = xi[k];
  else {			/* p2 = 0.0;  dp2 = 0.0; */
    dfac0 = -0.25 * (3.0 * chi - 1.0) * (chi - 1.0);
    dfac1 = (3.0 * chi + 1.0) * (chi - 1.0);
    /* dfac2 = -0.25 * (9.0 * chi + 5.0) * (chi - 1.0); */
    /* dfac3 =  0.5  * (3.0 * (chi*chi - 1)); */
    return istep * (dfac0 * p0 + dfac1 * p1);
  }

  /* factors for the interpolation of the 1. derivative */
  dfac0 = -(1.0 / 6.0) * ((3.0 * chi - 6.0) * chi + 2.0);
  dfac1 = 0.5 * ((3.0 * chi - 4.0) * chi - 1.0);
  dfac2 = -0.5 * ((3.0 * chi - 2.0) * chi - 2.0);
  dfac3 = 1.0 / 6.0 * (3.0 * chi * chi - 1.0);

  /* return the derivative */
  return istep * (dfac0 * p0 + dfac1 * p1 + dfac2 * p2 + dfac3 * p3);
}


/*****************************************************************************
*
*  Evaluate potential with quadratic interpolation.
*  col is typ1 * ntypes + typ2.
*
******************************************************************************/

real pot2(pot_table_t *pt, int col, real r)
{
  real  rr, istep, chi, p0, p1, p2, dv, d2v;
  int   k;

  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0)
    error("short distance!");

  /* indices into potential table */
  istep = pt->invstep[col];
  k = (int)(rr * istep);
  chi = (rr - k * pt->step[col]) * istep;
  k += pt->first[col];

  /* intermediate values */
  p0 = (k <= pt->last[col]) ? pt->table[k++] : 0.0;
  p1 = (k <= pt->last[col]) ? pt->table[k++] : 0.0;
  p2 = (k <= pt->last[col]) ? pt->table[k++] : 0.0;
  dv = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* return the potential value */
  return p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;
}


/*****************************************************************************
*
*  Evaluate potential with cubic interpolation.
*  col is typ1 * ntypes + typ2.
*
******************************************************************************/

real pot3(pot_table_t *pt, int col, real r)
{
  real  rr, istep, chi, p0, p1, p2, p3;
  real  fac0, fac1, fac2, fac3;
  int   k;

  /* check for distances shorter than minimal distance in table */
  rr = r - pt->begin[col];
  if (rr < 0)
    error("short distance!");

  /* indices into potential table */
  istep = pt->invstep[col];
  k = (int)(rr * istep);
  if (k == 0)
    return pot2(pt, col, r);	/* parabolic fit if on left border */
  chi = (rr - k * pt->step[col]) * istep;
  k += pt->first[col];
  k--;

  /* intermediate values */
  if (k <= pt->last[col])
    p0 = pt->table[k++];
  else
    return 0.0;
  if (k <= pt->last[col])
    p1 = pt->table[k++];
  else
    return 0.0;
  if (k <= pt->last[col])
    p2 = pt->table[k++];
  else
    return 0.0;
  if (k <= pt->last[col])
    p3 = pt->table[k];
  else {			/* p2 = 0.0; dp2 = 0.0 */
    fac0 = -0.25 * chi * SQR(chi - 1.0);
    fac1 = (chi * chi - 1) * (chi - 1);
    /* fac2 = -0.25 * chi * (chi + 1) * (3.0*chi - 5.0); */
    /* fac3 = -0.5  * (chi*chi - 1) * chi;           */
    return fac0 * p0 + fac1 * p1;
  }				/* go smoothly: interpolate with f=f'=0 at chi=1. */

  /* factors for the interpolation */
  fac0 = -(1.0 / 6.0) * chi * (chi - 1.0) * (chi - 2.0);
  fac1 = 0.5 * (chi * chi - 1.0) * (chi - 2.0);
  fac2 = -0.5 * chi * (chi + 1.0) * (chi - 2.0);
  fac3 = (1.0 / 6.0) * chi * (chi * chi - 1.0);

  /* return the potential value */
  return fac0 * p0 + fac1 * p1 + fac2 * p2 + fac3 * p3;
}
#endif /* OLDCODE */

#ifdef PARABEL
/*****************************************************************************
*
*  Evaluate value from parabole through three points.
*  Extrapolates for all k.
*
******************************************************************************/

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

/*****************************************************************************
*
*  Evaluate value from parabole through three points.
*  Extrapolates for all k. Nonequidistant points.
*
******************************************************************************/

real parab_ne(pot_table_t *pt, real *xi, int col, real r)
{
  real  x0, x1, x2, chi0, chi1, chi2, p0, p1, p2;
  int   k;

  /* renorm to beginning of table */
//  rr = r - pt->begin[col];
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
//  dv  = p1 - p0;
//  d2v = p2 - 2 * p1 + p0;

  /* return the potential value */
  return chi1 * chi2 * p0 - chi0 * chi2 * p1 + chi0 * chi1 * p2;

}

/*****************************************************************************
*
*  Evaluate deritvative from parabole through three points.
*  Extrapolates for all k.
*
******************************************************************************/

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

/*****************************************************************************
*
*  Evaluate deritvative from parabole through three points.
*  Extrapolates for all k.
*
******************************************************************************/

real parab_grad_ne(pot_table_t *pt, real *xi, int col, real r)
{
  real  h0, h1, h2, x0, x1, x2, chi0, chi1, chi2, p0, p1, p2, dv, d2v;
  int   k;

  /* renorm to beginning of table */
//  rr = r - pt->begin[col];
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

  /* intermediate values */
//  dv  = p1 - p0;
//  d2v = p2 - 2 * p1 + p0;

  /* return the potential value */
  return (chi2 / h1 + chi1 / h2) * p0
    - (chi0 / h2 + chi2 / h0) * p1 + (chi0 / h1 + chi1 / h0) * p2;

}

/*****************************************************************************
*
*  Evaluate value and deritvative from parabole through three points.
*  Extrapolates for all k.
*
******************************************************************************/

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

/*****************************************************************************
*
*  Evaluate value and deritvative from parabole through three points.
*  Extrapolates for all k.
*
******************************************************************************/

real parab_comb_ne(pot_table_t *pt, real *xi, int col, real r, real *grad)
{
  real  h0, h1, h2, x0, x1, x2, chi0, chi1, chi2, p0, p1, p2, dv, d2v;
  int   k;

  /* renorm to beginning of table */
//  rr = r - pt->begin[col];
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

  /* intermediate values */
//  dv  = p1 - p0;
//  d2v = p2 - 2 * p1 + p0;

  /* return the potential value */
  *grad = (chi2 / h1 + chi1 / h2) * p0
    - (chi0 / h2 + chi2 / h0) * p1 + (chi0 / h1 + chi1 / h0) * p2;

  return chi1 * chi2 * p0 - chi0 * chi2 * p1 + chi0 * chi1 * p2;
}

#endif

#ifdef APOT

void write_apot_table(apot_table_t *apt, char *filename)
{
  int   i, j;
  FILE *outfile;
  char  msg[255];

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile) {
    sprintf(msg, "Could not open file %s\n", filename);
    error(msg);
  }

  /* write header */
  fprintf(outfile, "#F 0 %d", apt->number);
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
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, " %d", invar_pot[i]);
  }
  fprintf(outfile, "\n#E\n\n");

  if (!disable_cp) {
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, "cp_%s %f %f %f\n", elements[i], apt->chempot[i],
	      apt->pmin[apt->number][i], apt->pmax[apt->number][i]);
    if (compnodes > 0)
      fprintf(outfile, "cn %d\n", compnodes);
    for (j = 0; j < compnodes; j++)
      fprintf(outfile, "%f %f %f %f\n", compnodelist[j],
	      apt->chempot[ntypes + j], apt->pmin[apt->number][ntypes + j],
	      apt->pmax[apt->number][ntypes + j]);
    fprintf(outfile, "\n");
  }

  /* write data */
  for (i = 0; i < apt->number; i++) {
    if (smooth_pot[i]) {
      fprintf(outfile, "type %s_sc\n", apt->names[i]);
    } else {
      fprintf(outfile, "type %s\n", apt->names[i]);
    }
    fprintf(outfile, "cutoff %f\n", apot_table.end[i]);
    fprintf(outfile, "# rmin %f\n", apt->begin[i]);
    if (smooth_pot[i]) {
      fprintf(outfile, "# r_0=%f", apt->co_pol[i][0]);
      for (j = 1; j < 5; j++)
	fprintf(outfile, " a_%d=%f", j - 1, apt->co_pol[i][j]);
      fprintf(outfile, "\n");
    }
    for (j = 0; j < apt->n_par[i]; j++) {
      fprintf(outfile, "%s %.10f %f %f\n", apt->param_name[i][j],
	      apt->values[i][j], apt->pmin[i][j], apt->pmax[i][j]);
    }
    if (i != (apt->number - 1))
      fprintf(outfile, "\n");
  }
  fclose(outfile);
}
#endif

/*****************************************************************************
*
*  write potential table (format 3)
*
******************************************************************************/

void write_pot_table3(pot_table_t *pt, char *filename)
{
  FILE *outfile, *outfile2;
  char  msg[255];
  int   i, j, flag = 0;
  real  r;

  if (plotpointfile != "\0")
    flag = 1;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile) {
    sprintf(msg, "Could not open file %s\n", filename);
    error(msg);
  }

  /* if needed: open file for plotpoints */
  if (flag) {
    outfile2 = fopen(plotpointfile, "w");
    if (NULL == outfile) {
      sprintf(msg, "Could not open file %s\n", filename);
      error(msg);
    }
  }

  /* write header */
  fprintf(outfile, "#F 3 %d", pt->ncols);
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
    fprintf(outfile, "%.16e %.16e %d\n",
	    pt->begin[i], pt->end[i], pt->last[i] - pt->first[i] + 1);
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

/*****************************************************************************
*
*  write potential table (format 4)
*
******************************************************************************/

void write_pot_table4(pot_table_t *pt, char *filename)
{
  FILE *outfile, *outfile2;
  char  msg[255];
  int   i, j, flag = 0;
  real  r;

  if (plotpointfile != "\0")
    flag = 1;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile) {
    sprintf(msg, "Could not open file %s\n", filename);
    error(msg);
  }

  /* if needed: open file for plotpoints */
  if (flag) {
    outfile2 = fopen(plotpointfile, "w");
    if (NULL == outfile) {
      sprintf(msg, "Could not open file %s\n", filename);
      error(msg);
    }
  }

  /* write header */
  fprintf(outfile, "#F 4 %d", pt->ncols);
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

#endif /* POTSCALE */

/*****************************************************************************
*
*  write potential table for IMD (format 2)
*
******************************************************************************/

void write_pot_table_imd(pot_table_t *pt, char *prefix)
{
  FILE *outfile;
  char  msg[255];
  char  filename[255];
  real *r2begin, *r2end, *r2step, temp2, r2, temp, root, x0;
  int   i, j, k, m, m2, col1, col2;
#ifdef APOT
  real  params[4];
#endif

  sprintf(filename, "%s_phi.imd.pt", prefix);
  /* allocate memory */
  r2begin = (real *)malloc(ntypes * ntypes * sizeof(real));
  r2end = (real *)malloc(ntypes * ntypes * sizeof(real));
  r2step = (real *)malloc(ntypes * ntypes * sizeof(real));
  if ((r2begin == NULL) || (r2end == NULL) || (r2step == NULL))
    error("Cannot allocate memory in  write_pot_table_imd");

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile) {
    sprintf(msg, "Could not open file %s\n", filename);
    error(msg);
  }

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
      fprintf(outfile, "%.16e %.16e %.16e\n",
	      r2begin[col2], r2end[col2], r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

#ifndef APOT
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
#ifdef NEWSCALE
	/* Pair potentials corrected so that U'(1)   =0 with NORESCALE */
	/*                               and U'(n_av)=0 without */
	fprintf(outfile, "%.16e\n",
		splint_ne(pt, pt->table, col1,
			  sqrt(r2)) + (sqrt(r2) <=
				       pt->end[paircol + j] ? lambda[i] *
				       splint_ne(pt, pt->table, paircol + j,
						 sqrt(r2)) : 0.) +
		(sqrt(r2) <= pt->end[paircol + i] ? lambda[j] *
		 splint_ne(pt, pt->table, paircol + i, sqrt(r2)) : 0.));
#else /* NEWSCALE */
	fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, sqrt(r2)));
#endif /* NEWSCALE */
	r2 += r2step[col2];
      }
      fprintf(outfile, "%.16e\n", 0.0);
      fprintf(outfile, "\n");
    }
  }
#else
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col2 = i * ntypes + j;
      r2 = r2begin[col2];
      if (do_smooth && smooth_pot[col1])
	x0 =
	  smooth(apot_table.fvalue[col1], apot_table.end[col1],
		 apot_table.values[col1], params);
      else
	x0 = rcut[col2];
      for (k = 0; k < imdpotsteps; k++) {
	if (sqrt(r2) < x0)
	  apot_table.fvalue[col1] (sqrt(r2), apot_table.values[col1], &temp);
	else
	  temp =
	    params[0] * r2 * sqrt(r2) + params[1] * r2 +
	    params[2] * sqrt(r2) + params[3];
	fprintf(outfile, "%.16e\n", temp);
	r2 += r2step[col2];
      }
      fprintf(outfile, "%.16e", 0.0);
      if (!((i == (ntypes - 1)) && (j == (ntypes - 1))))
	fprintf(outfile, "\n\n\n");
    }
  }
#endif
  fclose(outfile);
  printf("IMD pair potential data written to %s\n", filename);
#ifdef EAM
  /* write rho_r2 */
  sprintf(filename, "%s_rho.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile) {
    sprintf(msg, "Could not open file %s\n", filename);
    error(msg);
  }

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes * ntypes);

  /* write info block */
  for (i = 0; i < ntypes; i++) {
    for (j = 0; j < ntypes; j++) {
      col1 = (ntypes * (ntypes + 1)) / 2 + j;
      col2 = i * ntypes + j;
      /* Extrapolation possible  */
      r2begin[col2] = SQR(MAX(pt->begin[col1] - extend * pt->step[col1], 0));
      r2end[col2] = SQR(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n",
	      r2begin[col2], r2end[col2], r2step[col2]);
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
	fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, sqrt(r2)));
	r2 += r2step[col2];
      }
      fprintf(outfile, "%.16e\n", 0.0);
      fprintf(outfile, "\n");
    }
  }
  fclose(outfile);
  printf("IMD electron transfer data written to %s\n", filename);
  /* write F_rho */
  sprintf(filename, "%s_F.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile) {
    sprintf(msg, "Could not open file %s\n", filename);
    error(msg);
  }

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes);

  /* write info block */
  for (i = 0; i < ntypes; i++) {
    col1 = (ntypes * (ntypes + 3)) / 2 + i;
    /* pad with zeroes */
    r2begin[i] = pt->begin[col1] - extend * pt->step[col1];
    /* extrapolation */
    r2end[i] = pt->end[col1] + extend * pt->step[col1];
    r2step[i] = (r2end[i] - r2begin[i]) / imdpotsteps;
    fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[i], r2end[i], r2step[i]);
  }
  fprintf(outfile, "\n");

  /* write data */
  for (i = 0; i < ntypes; i++) {
    r2 = r2begin[i];
    col1 = (ntypes * (ntypes + 3)) / 2 + i;
    root = (pt->begin[col1] > 0) ?
      pt->table[pt->first[col1]] / sqrt(pt->begin[col1]) : 0.;
    root += (pt->end[col1] < 0) ?
      pt->table[pt->last[col1]] / sqrt(-pt->end[col1]) : 0;
    for (k = 0; k <= imdpotsteps; k++) {
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
#ifdef REPULSE
      temp2 = r2 - pt->end[col1];
      temp += (temp2 > 0.) ? 5e2 * (temp2 * temp2 * temp2) : 0.;
#endif
#ifdef NEWSCALE
      temp -= lambda[i] * r2;
#endif /* NEWSCALE */
      fprintf(outfile, "%.16e\n", temp);
      r2 += r2step[i];
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);
  printf("IMD embedding data written to %s\n", filename);
#endif
  free(r2begin);
  free(r2end);
  free(r2step);
}

/*****************************************************************************
*
*  write plot version of potential table
*
******************************************************************************/

void write_plotpot_pair(pot_table_t *pt, char *filename)
{
  FILE *outfile;
  char  msg[255];
  int   i, j, k, l, col;
  real  r, r_step, temp, x0;
  real  params[4];

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile) {
    sprintf(msg, "Could not open file %s\n", filename);
    error(msg);
  }

  /* write data */
  k = 0;

#ifndef APOT
  for (i = 0; i < ntypes; i++)
    for (j = i; j < ntypes; j++) {
      r = pt->begin[k];
      r_step = (pt->end[k] - r) / (NPLOT - 1);
      for (l = 0; l < NPLOT - 1; l++) {
#ifdef NEWSCALE
	fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, k, r)
		+ (r <= pt->end[paircol + i] ?
		   splint_ne(pt, pt->table, paircol + i, r) * lambda[j] : 0.)
		+ (r <=
		   pt->end[paircol + j] ? splint_ne(pt, pt->table,
						    paircol + j,
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
#endif
#else
  for (i = 0; i < apot_table.number; i++) {
    if (i < paircol + ntypes)
      r = (plotmin == 0 ? 0.1 : plotmin);
    else
      r = apot_table.begin[i];
    j = i * (i + 1) / 2;
    if (!invar_pot[i] && smooth_pot[i]) {
      x0 =
	smooth(apot_table.fvalue[i], apot_table.end[i], apot_table.values[i],
	       params);
      r_step = (apot_table.end[i] - r) / (NPLOT - 1);
    } else {
      x0 = 0;
      r_step = (apot_table.end[i] - r) / (NPLOT - 1);
    }
    for (l = 0; l < NPLOT; l++) {
      if (x0 == 0 || r <= x0)
	apot_table.fvalue[i] (r, apot_table.values[i], &temp);
      else
	temp =
	  params[0] * r * r * r + params[1] * r * r + params[2] * r +
	  params[3];
      fprintf(outfile, "%e %e\n", r, temp);
      r += r_step;
    }
    if (i != (apot_table.number - 1))
      fprintf(outfile, "\n\n");
  }

#endif
  fclose(outfile);
  printf("Potential plotting data written to %s\n", filename);
}

/*****************************************************************************
*
*  write alternate plot version of potential table
*  (same intervals for all pair and transfer functions)
*
******************************************************************************/

void write_altplot_pair(pot_table_t *pt, char *filename)
{
  FILE *outfile;
  char  msg[255];
  int   i, j, k, l, col;
  real  r, rmin = 100., rmax = 0., r_step, temp;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile) {
    sprintf(msg, "Could not open file %s\n", filename);
    error(msg);
  }

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
      col = i * ntypes + j;
      r = rmin;
      for (l = 0; l < NPLOT - 1; l++) {
#ifdef NEWSCALE
	fprintf(outfile, "%e %e\n", r,
		(r <= pt->end[k] ? splint_ne(pt, pt->table, k, r) : 0.)
		+ (r <= pt->end[paircol + i] ?
		   splint_ne(pt, pt->table, paircol + i, r) * lambda[j] : 0.)
		+ (r <=
		   pt->end[paircol + j] ? splint_ne(pt, pt->table,
						    paircol + j,
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
      fprintf(outfile, "%e %e\n", r,
	      r <= pt->end[i] ? splint_ne(pt, pt->table, i, r) : 0);
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
/********************************************************************
 *
 * write_pairdist(pot_table_t *pt, char *filename)
 *    - write distribution function of function access
 *
 *
 *******************************************************************/

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
  if (NULL == outfile) {
    sprintf(msg, "Could not open file %s\n", filename);
    error(msg);
  }


  /* Verteilungsfeld initialisieren */
  freq = (int *)malloc(ndimtot * sizeof(int));
  for (i = 0; i < ndimtot; i++)
    freq[i] = 0;

  for (h = firstconf; h < firstconf + myconf; h++) {
    for (i = 0; i < inconf[h]; i++) {
      atom = atoms + i + cnfstart[h];
      typ1 = atom->typ;

      /* Paarpotenzialfunktion */
      for (j = 0; j < atom->n_neigh; j++) {
	neigh = atom->neigh + j;
	typ2 = neigh->typ;
	col = (typ1 <= typ2) ?
	  typ1 * ntypes + typ2 - ((typ1 * (typ1 + 1)) / 2)
	  : typ2 * ntypes + typ1 - ((typ2 * (typ2 + 1)) / 2);
	/* Die Arbeit wurde bereits gemacht */
	if (neigh->r < pt->end[col])
	  freq[neigh->slot[0]]++;
#ifdef EAM
	/* Transferfunktion */
	col = paircol + typ2;
	if (neigh->r < pt->end[col])
	  freq[neigh->slot[1]]++;
#endif /* EAM */
      }
#ifdef EAM
      /* Finally: Einbettungsfunktion - hier muss Index festgestellt werden */
      col = paircol + ntypes + typ1;
      if (format == 3) {
	rr = atom->rho - pt->begin[col];
#ifdef NORESCALE
	if (rr < 0)
	  rr = 0;		/* Extrapolation */
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
  /* OK, jetzt haben wir die Daten - schreiben wir sie raus */
  j = 0;
//  rr=0.5*(pt->begin[0]+pt->xcoord[1]);
  for (col = 0; col < pt->ncols; col++) {
    for (i = pt->first[col]; i < pt->last[col]; i++) {
      rr = 0.5 * (pt->xcoord[i] + pt->xcoord[i + 1]);
      fprintf(outfile, "%f %d\n", rr, freq[i]);
    }
    fprintf(outfile, "\n\n");
  }
  fclose(outfile);
  printf("Distribution data written to %s\n", filename);
}
#endif
