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

/*added to use OptParamType*/
#include "kim/kim.h"
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


/*added (modified )*/
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


/* added */
#else /* !KIM */
			if(format != 5)
				error(1,"The potential file format in file (%s) is not compatiable with KIM. To use "
								"potfit-KIM, the format should be 5.", filename);
      /* although we require the user to use format 5, but internally, we still use
       * either format 0 (when ``nolimits'' is not in the compilation target), and
       * use format 3 (when it is in the compilation target). We do this because
       * format are used in ohter places, and we do not want to make that changes */
#ifndef NOLIMITS
      format = 0;
#else /* !NOLIMITS */
      format = 3;
#endif /* !NOLIMITS */

			if(size != 1) {
				warning("The number of potentials should always be 1 when KIM Model is used. "
								"You specifiy %d in file (%s), and it is reset to 1.\n", size,filename);
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


/*added */
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



#endif/* !KIM */
/*added ends*/



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
/*added modified (allocate more memory)*/
			/* apt->names[i] = (char *)malloc(20 * sizeof(char));*/
			apt->names[i] = (char *)malloc(255 * sizeof(char));
/*added ends*/

		if ((apt->n_par == NULL) || (apt->begin == NULL) || (apt->end == NULL)
				|| (apt->fvalue == NULL) || (apt->names == NULL) || (apt->pmin == NULL)
				|| (apt->pmax == NULL) || (apt->param_name == NULL)
				|| (apt->values == NULL))
			error(1, "Cannot allocate info block for analytic potential table %s", filename);
#endif /* APOT */


		switch (format) {
#ifdef APOT
			case 0:
				read_pot_table0(pt, apt, filename, infile);
				break;
#else
			case 3:
				read_pot_table3(pt, size, filename, infile);
				break;
/*added (we delete the read_pot_table4 in this file)*/
/*			case 4:
				read_pot_table4(pt, size, filename, infile);
*/
#endif /* APOT */
		}


		fclose(infile);

		printf("Reading potential file >> %s << ... done\n", filename);


/* added (modified)*/
		/* copy cutoff, cutoff equal pt->end[0] */
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
				rmin[i * ntypes + j] = 0.0;
				rcut[i * ntypes + j] = pt->end[0];
			}
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
	double  tmp_pmin, tmp_pmax;
	char tmp_value[255];
	int jj, kk, tmp_size;
	OptParamType OptParamSet;


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
	ReadPotentialKeywords(pt, filename, infile, &OptParamSet);

	/* There is only 1 potential for KIM Model. */
	i = 0;
	/* get the total size (some parameters may be array) of the optimizable 
	 * parameters and the `nestedvalue'. For analytic potential, apt->n_par[i]
	 * need to be equal to num_opt_param, since they should all be scalar. */
	apt->n_par[i] = get_OptimizableParamSize(&OptParamSet, name_opt_param, num_opt_param); 
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
	apt->param_name[i] = (char **)malloc(apt->n_par[i]* sizeof(char *));
	reg_for_free(apt->param_name[i], "apt->param_name[%d]", i);
	if (NULL == apt->values[i] || NULL == apt->invar_par[i] || NULL == apt->param_name[i])
		error(1, "Cannot allocate memory for potential paramters.\nAborting");
	for (j = 0; j < apt->n_par[i]; j++) {
		apt->param_name[i][j] = (char *)malloc(255 * sizeof(char));
		reg_for_free(apt->param_name[i][j], "apt->param_name[%d][%d]", i, j);
		if (NULL == apt->param_name[i][j]) {
			error(1, "Error in allocating memory for parameter name");
		}
	}

	/* allocate memory */
	apt->pmin[i] = (double *)malloc(apt->n_par[i] * sizeof(double));
	reg_for_free(apt->pmin[i], "apt->pmin[%d]", i);
	apt->pmax[i] = (double *)malloc(apt->n_par[i] * sizeof(double));
	reg_for_free(apt->pmax[i], "apt->pmax[%d]", i);
	if(NULL == apt->pmin[i] || NULL == apt->pmax[i])
		error(1, "Cannot allocate memory for potential paramters.\nAborting");


	/*  scan to read parameter values */
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
		/* how many parameter(s) are in this keyword? ( = 1 if scalar) */
		/* k is the palce(slot) that the keyword is in the OptParam struct */
		for ( k = 0; k < OptParamSet.Nparam; k++) {
			if (strcmp(OptParamSet.name[k], name_opt_param[j]) == 0)
				break;
		}
		tmp_size = 1;
		if( OptParamSet.rank[k] != 0) {
			for ( kk = 0; kk < OptParamSet.rank[k]; kk++) {
				tmp_size *= OptParamSet.shape[k][kk];
			}
		}
		/* read in the value and pmin and pmax */
		for (kk = 0; kk < tmp_size; kk++) {
			fgets(buffer, 255, infile);
			ret_val = sscanf(buffer, "%s %lf %lf", tmp_value, &tmp_pmin, &tmp_pmax);
			if (ret_val == 3) {
				apt->pmin[i][jj] = tmp_pmin;
				apt->pmax[i][jj] = tmp_pmax;
				/*if give `KIM', use KIM parameter, otherwise use readin */
				if (strncmp(tmp_value, "KIM", 3) == 0) {
					apt->values[i][jj] = *OptParamSet.nestedvalue[jj];
				} else if (1 > sscanf(tmp_value, "%lf", &apt->values[i][jj])) {
					error(1, " Value %d of (%s) is not float.\n", kk+1, OptParamSet.name[k]);
				}
			} else {
				error(0, "Not enough value(s) for (%s) are provided in %s. "
						"You listed %d value(s), but required are %d.\n", 
						OptParamSet.name[k], filename, kk, tmp_size);
				error(1, "Or line %d of (%s) are of wrong type. Only `KIM' and float "
						"data are acceptable.\n", kk+1,OptParamSet.name[k]);
			}
			/* If the parameter is array, the parameter name is given to the first value 
			 * and the other name will be given `\0'. */
			if(kk == 0) {
				strcpy(apt->param_name[i][jj], name_opt_param[j]);
			} else {
				strcpy(apt->param_name[i][jj], '\0');
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
	/* If it is anaytic potential (every PARAM_FREE* is scalar), keyword `cutoff' has
	 * to be specified in the input file. 
	 * You can give value directly, or give `KIM' to use the cutoff from KIM
	 * model. If it is tabualted potential, cutoff will read in from KIM model, even
	 * if it is given in the input file. */
	
	int have_cutoff;

	fsetpos(infile, &startpos);
	do {
		fgets(buffer, 255, infile);
		sscanf(buffer, "%s", name);
	} while (strcmp(name, "cutoff") != 0 && !feof(infile));

	if (apt->n_par[i] == num_opt_param) {
		printf("Each `PARAM_FREE_*' is scalar.\n");
		if (strcmp(name, "cutoff") == 0) {
			if(2 > sscanf(buffer, "%s %s", name, tmp_value)) {
				error(1,"Error reading in `cutoff'\n");
			} else {
				if (strcmp(tmp_value, "KIM") == 0) {
					have_cutoff = 0;
				} else {
					if(1 > sscanf(tmp_value,"%lf", &apt->end[i])){
						error(1,"Error reading in `cutoff'.\n");
					}
					have_cutoff = 1; 
				}
			}
		} else {
			error(1,"`cutoff' is missing in %s.\n", filename);
		}
	} else {
		printf("Some `PARAM_FREE_*' are array; (%s) is probably a tabulated potential.\n"
				,kim_model_name);
		if (strcmp(name, "cutoff") == 0) {
			printf("`cutoff' read in from %s is deprecated. Will use the `cutoff' in "
					" the KIM model throughout the fitting.\n", filename);
		}
		have_cutoff = 0;
	}

	/* get cutoff from KIM if cutoff has not been determined yet. */
	if (!have_cutoff) {
		for (j = 0; j < OptParamSet.Nparam; j++) { 
			if (strncmp(OptParamSet.name[j], "PARAM_FREE_cutoff", 17) == 0 ) {
				apt->end[i] = *OptParamSet.value[j];
				have_cutoff = 1;
				break;
			}
		}
		if(!have_cutoff) {
			error(1, "`PARAM_FREE_cutoff' cannot be found in `descriptor.kim'. Copy "
					"cutoff to potfit failed.");
		}
	}


	printf(" - Successfully read %d potential table(s)\n", apt->number);


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

	/*added don't need it any more; the calculation are done by KIM. But we still
		include it to make it possible to use the config.c file, where a lot of
		calc_table info are needed   */
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
	int   i, j, k, ret_val;
	char  buffer[255], name[255];
	fpos_t filepos, startpos; 
	OptParamType OptParamSet;

	/* save starting position */
	fgetpos(infile, &startpos);

	/* read the keywords and the names of parameters that will be optimized */
	ReadPotentialKeywords(pt, filename, infile, &OptParamSet);


	/* get the total size (some parameters may be array) of the optimizable 
	 * parameters and the `nestedvalue'. For analytic potential, apt->n_par[i]
	 * need to be equal to num_opt_param, since they should all be scalar. */
	get_OptimizableParamSize(&OptParamSet, name_opt_param, num_opt_param); 

	/* some potential table value */
	pt->first[0] = 0;
	pt->last[0] = pt->first[0] + OptParamSet.Nnestedvalue - 1;
	pt->len = OptParamSet.Nnestedvalue;
	pt->idxlen = OptParamSet.Nnestedvalue;

	/* allocate the function table */
	pt->table = (double *)malloc(pt->len * sizeof(double));
	reg_for_free(pt->table, "pt->table");
	pt->idx = (int *)malloc(pt->len * sizeof(int));
	reg_for_free(pt->idx, "pt->idx");
	if ((NULL == pt->table) || (NULL == pt->idx) )
		error(1, "Cannot allocate memory for potential table");

	/* copy parameters values to potfit */
	for (i = 0; i < OptParamSet.Nnestedvalue; i++) {
		pt->table[i] = *OptParamSet.nestedvalue[i];
		pt->idx[i] = i;
	}

	/* set end (end is cutoff) */
	int have_cutoff = 0;
	double tmp_cutoff;
	for (j = 0; j < OptParamSet.Nparam; j++) { 
		if (strncmp(OptParamSet.name[j], "PARAM_FREE_cutoff", 17) == 0 ) {
			pt->end[0] = *OptParamSet.value[j];
			have_cutoff = 1;
			break;
		}
	}
	if(!have_cutoff) {
		error(1, "`PARAM_FREE_cutoff' cannot be found in `descriptor.kim'. copy "
				"cutoff to potfit failed.");
	}

	printf(" - Successfully read potential parameters that will be optimized.\n");

	return;
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
						/*added function pointer is disabled, because we don't need it. The potential
							comes from KIM directly now. */
						/*	      apot_table.fvalue[i] (calct->xcoord[index], val, &f);
						 */
						/*added ends*/
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

		/*added disable it*/
		/*    (*val) = apot_grad(calc_pot.begin[i], val + 2, apot_table.fvalue[i]);
		 */
		/*added ends*/
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
				/*added  disable it*/
				/*  apot_table.fvalue[i] (calc_pot.xcoord[k], val, &f);
				 */
				/*added ends*/
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
