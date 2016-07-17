/****************************************************************
 *
 * potential_input_f5.c: Routines for reading a potential table
 *
 ****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
 *
 * http://potfit.sourceforge.net/
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

typedef struct {
  char const* filename;
  FILE* pfile;
  fpos_t startpos;
} apot_state;

void read_chemical_potentials(apot_state* pstate);
void read_pot_table5(char const* filename, FILE* pfile); 
#if !defined(NOLIMITS)
void read_pot_table5_limits(char const* filename, FILE* pfile);
#else
void read_pot_table5_nolimits(char const* filename, FILE* pfile);
#endif

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
  char name[255];

  fpos_t filepos;

  if (g_param.enable_cp) {
    /* search for cp */
    do {
      fgetpos(pstate->pfile, &filepos);
      fscanf(pstate->pfile, "%s", buffer);
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
    if (2 > fscanf(pfile, "%s %d", buffer, &g_param.compnodes)) {
      if (strcmp("type", buffer) == 0)
        g_param.compnodes = -1;
      else
        error(1,
              "Could not read number of composition nodes from potential "
              "file.\n");
    }
    if (strcmp(buffer, "cn") != 0 && g_param.ntypes > 1 &&
        g_param.g_param.compnodes != -1)
      error(1, "No composition nodes found in %s.\nUse \"cn 0\" for none.\n",
            pstate->filename);
    if (g_param.ntypes == 1) {
      g_param.compnodes = 0;
    }
    if (g_param.compnodes != -1) {
      apt->values[apt->number] = (double*)Realloc( apt->values[apt->number],
          (g_param.ntypes + g_param.g_param.compnodes) * sizeof(double));
      apt->pmin[apt->number] = (double*)Realloc( apt->pmin[apt->number],
          (g_param.ntypes + g_param.g_param.compnodes) * sizeof(double));
      apt->pmax[apt->number] = (double*)Realloc( apt->pmax[apt->number],
          (g_param.ntypes + g_param.g_param.compnodes) * sizeof(double));
      apt->chempot = apt->values[apt->number];
      compnodelist = (double*)Malloc((g_param.ntypes + g_param.g_param.compnodes) *
                                     sizeof(double));

      for (j = 0; j < g_param.compnodes; j++) {
        if (4 > fscanf(pfile, "%lf %lf %lf %lf", &compnodelist[j],
                       &apt->chempot[g_param.ntypes + j],
                       &apt->pmin[apt->number][g_param.ntypes + j],
                       &apt->pmax[apt->number][g_param.ntypes + j]))
          error(1, "Could not read composition node %d\n", j + 1);
        if (apt->pmin[apt->number][g_param.ntypes + j] > apt->chempot[g_param.ntypes + j] ||
            apt->pmax[apt->number][g_param.ntypes + j] < apt->chempot[g_param.ntypes + j])
          error(1, "composition node %d is out of bounds.\n", j + 1);
      }

      /* check g_param.compnodes for valid values */
      if (g_param.ntypes == 2) {
        for (j = 0; j < g_param.compnodes; j++)
          if (compnodelist[j] > 1 || compnodelist[j] < 0)
            error(1, "Composition node %d is %f but should be inside [0,1].\n",
                  j + 1, compnodelist[j]);
      }
    }
    if (g_param.compnodes != -1)
      printf("Enabled chemical potentials with %d extra composition node(s).\n",
             g_param.compnodes);
    if (g_param.compnodes == -1)
      g_param.compnodes = 0;
#endif  // CN
  }
#endif  // PAIR
}


void read_pot_table5(char const* filename, FILE* pfile)
{
#if !defined(NOLIMITS)
  read_pot_table5_limits(filename, pfile);  
#else
	read_pot_table5_nolimits(filename, pfile);
#endif

}


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
#if !defined(NOLIMITS)
void read_pot_table5_limits(char const* filename, FILE* pfile)
{
  apot_table_t* apt = &g_pot.apot_table;
  pot_table_t* pt = &g_pot.opt_pot;
  apot_state state;
  state.filename = filename;
  state.pfile = pfile;


	int   i, j, k, l, ret_val;
	char  buffer[255], name[255];
	double *val, *list;
	fpos_t startpos;
	double  tmp_pmin, tmp_pmax;
	char tmp_value[255];
	int jj, kk;
	void* pkim;
	int status;
	FreeParamType FreeParamSet;


	/* save starting position */
	fgetpos(pfile, &startpos);

  read_chemical_potentials(&state);

	/* read the keywords and the names of parameters that will be optimized */
	read_potential_keyword(pt, filename, pfile);

	/* write the descriptor.kim file for the test */
	write_temporary_descriptor_file(g_kim.kim_model_name);

	/* create KIM object with 1 atom and 1 species */
	status = setup_KIM_API_object(&pkim, 1, 1, g_kim.kim_model_name);
	if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__, "setup_KIM_API_object", status);
		exit(1);
	}

  get_compute_const(pkim);

  /* initialze the data struct for the free parameters with type double */
  get_free_param_double(pkim, &FreeParamSet);

  /* nest the optimizable params */
  nest_optimizable_param(pkim, &FreeParamSet, g_kim.name_opt_param, g_kim.num_opt_param);

	/* There is only 1 potential for KIM Model. */
	i = 0;
	apt->n_par[i] = FreeParamSet.Nnestedvalue; 
	apt->total_par += apt->n_par[i];
	
	/* copy name */
	strcpy(apt->names[i], g_kim.kim_model_name);

	/* set begin */
	apt->begin[i] = 0.0;

	/* allocate memory */
	apt->values[i] = (double *)Malloc(apt->n_par[i] * sizeof(double));
	/* last slot stores the total number of invariable parameters */
	apt->invar_par[i] = (int *)Malloc((apt->n_par[i] + 1) * sizeof(int));
	apt->pmin[i] = (double *)Malloc(apt->n_par[i] * sizeof(double));
	apt->pmax[i] = (double *)Malloc(apt->n_par[i] * sizeof(double));
	if(NULL == apt->pmin[i] || NULL == apt->pmax[i]
	   || NULL == apt->values[i] || NULL == apt->invar_par[i])
		error(1, "Cannot allocate memory for potential paramters.\nAborting");


	/* invariable parameters. Assume all are optimizable, and not optimiable ones
     will be marked later */
	/* the last slot is the number of invar */
	for (jj = 0; jj < apt->n_par[i]; jj++) {
		apt->invar_par[i][jj] = 0;
	}
	apt->invar_par[i][apt->n_par[i]] = 0;


	/* read parameter values */
  /* j: jth optimizable parameter 
     k: kth free parameter, index of jth optimizable param in the free param struct 
    jj: the current slot in the nested apot->value
    kk: the current slot of the jth (or kth) param 
  */
	jj = 0;  
	fsetpos(pfile, &startpos);
	for (j = 0; j < g_kim.num_opt_param; j++) {
		buffer[0] = '\0';
		name[0] = '\0';
		/* find the keyword `PARAM_FREE_*' */
		do {
			fgets(buffer, 255, pfile);
			sscanf(buffer, "%s", name);
		} while (strncmp(name, "PARAM_FREE", 10) != 0 && !feof(pfile));
		/* k is the palce(slot) that the param in the FreeParam struct */
		for ( k = 0; k < FreeParamSet.Nparam; k++) {
			if (strcmp(FreeParamSet.name[k], g_kim.name_opt_param[j]) == 0)
				break;
		}
		/* read in the value and pmin and pmax */
		for (kk = 0; kk < g_kim.size_opt_param[j]; kk++) {
			fgets(buffer, 255, pfile);
			ret_val = sscanf(buffer, "%s %lf %lf", tmp_value, &tmp_pmin, &tmp_pmax);
			if (ret_val == 3) {
				apt->pmin[i][jj] = tmp_pmin;
				apt->pmax[i][jj] = tmp_pmax;
				/*if give `KIM', use KIM parameter, otherwise use readin */
				if (strncmp(tmp_value, "KIM", 3) == 0) {
					apt->values[i][jj] = FreeParamSet.value[k][kk];
				} else if (1 > sscanf(tmp_value, "%lf", &apt->values[i][jj])) {
					error(1, "First data for parameter '%s' corrupted, it should be a type float or 'KIM'.\n",
                    g_kim.name_opt_param[j]);
				}
				if(apt->pmin[i][jj]>apt->values[i][jj] || apt->pmax[i][jj] < apt->values[i][jj])
				{
          error(1, "Value %d of '%s' is not within its limits in file '%s'.\n", kk+1, g_kim.name_opt_param[j], filename);
        }
        /* fixed param or not */
        ret_val = sscanf(buffer, "%*s %*s %*s %s", tmp_value);
        if (strncmp(tmp_value, "FIX", 3) == 0) {
          apt->invar_par[i][jj] = 1;
          apt->invar_par[i][apt->n_par[i]] += 1;
        }

			} else {
				error(0, "Data for parameter '%s' corrupted in file '%s'.\n",
                  g_kim.name_opt_param[j], filename);
        error(1, "Or you may want to try compiling with nolimits enabled.\n");
			}
			jj++;
		}  
	}


	printf("Successfully read potential parameters that will be optimized.\n\n");

	/* find cutoff */
	int have_cutoff;
  int have_PARAM_FREE_cutoff;
  have_cutoff = 0;	
  /* If every PARAM_FREE* is scalar keyword `cutoff' has
	 * to be specified in the input file. You can give value directly, or give `KIM' 
   * to use the cutoff from KIM model.
   * If it is tabualted potential, cutoff will read in from KIM model, even
	 * if it is given in the input file. */
	
  /* whether we have PARAM_FREE_cutoff in the KIM Model. 
  If not, cannot publish cutoff, then have_cutoff needs to be 0. */

  KIM_API_get_data(pkim, "PARAM_FREE_cutoff", &status);
   if (KIM_STATUS_OK == status)  have_PARAM_FREE_cutoff = 1;


	fsetpos(pfile, &startpos);
	do {
		fgets(buffer, 255, pfile);
		sscanf(buffer, "%s", name);
	} while (strcmp(name, "cutoff") != 0 && !feof(pfile));

	if (apt->n_par[i] == g_kim.num_opt_param) {
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
      exit(1);
    }
      apt->end[i] = *pcutoff;
  }

#ifdef PAIR
	if (g_param.enable_cp) {
		g_pot.cp_start = apt->total_par - apt->globals + g_param.ntypes * (g_param.ntypes + 1);
		apt->total_par += (g_param.ntypes + g_param.compnodes - apt->invar_par[apt->number][g_param.ntypes]);
	}
#endif /* PAIR */

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
	if (g_pot.have_globals)
		pt->len += apt->globals;

#ifdef PAIR
	if (g_param.enable_cp) {
		pt->len += (g_param.ntypes + g_param.compnodes);
	}
#endif /* PAIR */

	pt->table = (double *)Malloc(pt->len * sizeof(double));

	g_pot.calc_list = (double *)Malloc(pt->len * sizeof(double));

	pt->idx = (int *)Malloc(pt->len * sizeof(int));

	apt->idxpot = (int *)Malloc(apt->total_par * sizeof(int));

	apt->idxparam = (int *)Malloc(apt->total_par * sizeof(int));

	if ((NULL == pt->table) || (NULL == pt->idx) || (apt->idxpot == NULL)
			|| (apt->idxparam == NULL))
		error(1, "Cannot allocate memory for potential table.\n");
	for (i = 0; i < pt->len; i++) {
		pt->table[i] = 0;
		g_pot.calc_list[i] = 0;
		pt->idx[i] = 0;
	}


	/* no derivative slots are needed, val += 2 etc. is disabled )  */  
	/* this is the indirect index */
	k = 0;
	l = 0;
	val = pt->table;
	list = g_pot.calc_list;
	for (i = 0; i < apt->number; i++) {	/* loop over potentials */
		for (j = 0; j < apt->n_par[i]; j++) {	/* loop over parameters */
			*val = apt->values[i][j];
			*list = apt->values[i][j];
			val++;
			list++;
			if (!g_pot.invar_pot[i] && !apt->invar_par[i][j]) {
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
			} else {
				l++;
      }
		}
		if (!g_pot.invar_pot[i])
			pt->idxlen += apt->n_par[i] - apt->invar_par[i][apt->n_par[i]];
		apt->total_par -= apt->invar_par[i][apt->n_par[i]];
	}

#ifdef PAIR
	if (g_param.enable_cp) {
		init_chemical_potential(g_param.ntypes);
		i = apt->number;
		for (j = 0; j < (g_param.ntypes + g_param.compnodes); j++) {
			*val = apt->values[i][j];
			val++;
			if (!apt->invar_par[i][j]) {
				pt->idx[k] = l++;
				apt->idxpot[k] = i;
				apt->idxparam[k++] = j;
			}
		}
		pt->idxlen += (g_param.ntypes + g_param.compnodes - apt->invar_par[apt->number][g_param.ntypes]);
		g_pot.global_idx += (g_param.ntypes + g_param.compnodes - apt->invar_par[apt->number][g_param.ntypes]);
	}
#endif /* PAIR */

#ifdef NOPUNISH
	if (opt)
		warning("Gauge degrees of freedom are NOT fixed!\n");
#endif /* NOPUNISH */

	/* free memory */
	free_model_object(&pkim);

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
void read_pot_table5_nolimits(char const* filename, FILE *pfile)
{
	int   i, j, k, jj, kk;
	char  buffer[255], name[255];
	fpos_t startpos; 
	char tmp_value[255];
	void* pkim;
  int status;
	FreeParamType FreeParamSet;

	/* save starting position */
	fgetpos(pfile, &startpos);

	/* read the keywords and the names of parameters that will be optimized */
	read_potential_keyword(pt, filename, pfile);

	/* write the descriptor.kim file for the test */
	write_temporary_descriptor_file(g_kim.kim_model_name);

	/* create KIM object with 1 atom and 1 species */
	status = setup_KIM_API_object(&pkim, 1, 1, g_kim.kim_model_name);
	if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__, "setup_KIM_API_object", status);
		exit(1);
	}

  get_compute_const(pkim);

  /* initialze the data struct for the free parameters with type double */
  get_free_param_double(pkim, &FreeParamSet);

  /* nest the optimizable params */
  nest_optimizable_param(pkim, &FreeParamSet, g_kim.name_opt_param, g_kim.num_opt_param);

	/* some potential table value */
	pt->first[0] = 0;
	pt->last[0] = pt->first[0] + FreeParamSet.Nnestedvalue - 1;
	pt->len = FreeParamSet.Nnestedvalue;
	pt->idxlen = FreeParamSet.Nnestedvalue;

	/* allocate the function table */
	pt->table = (double *)Malloc(pt->len * sizeof(double));
	pt->idx = (int *)Malloc(pt->len * sizeof(int));
	if ((NULL == pt->table) || (NULL == pt->idx) )
		error(1, "Cannot allocate memory for potential table");

	/* read parameter values */
  /* j: jth optimizable parameter 
     k: kth free parameter, index of jth optimizable param in the free param struct 
    jj: the current slot in the nested apot->value
    kk: the current slot of the jth param
  */
	jj = 0;  
	fsetpos(pfile, &startpos);
	for (j = 0; j < g_kim.num_opt_param; j++) {
		buffer[0] = '\0';
		name[0] = '\0';
		/* find the keyword `PARAM_FREE_*' */
		do {
			fgets(buffer, 255, pfile);
			sscanf(buffer, "%s", name);
		} while (strncmp(name, "PARAM_FREE", 10) != 0 && !feof(pfile));
		/* k is the palce(slot) that the param in the FreeParam struct */
		for ( k = 0; k < FreeParamSet.Nparam; k++) {
			if (strcmp(FreeParamSet.name[k], g_kim.name_opt_param[j]) == 0)
				break;
		}
		for (kk = 0; kk < g_kim.size_opt_param[j]; kk++) {
			fgets(buffer, 255, pfile);
			if ( 1 == sscanf(buffer, "%s", tmp_value)) {
				/*if give `KIM', use KIM parameter, otherwise use readin */
				if (strncmp(tmp_value, "KIM", 3) == 0) {
					pt->table[jj] = *FreeParamSet.nestedvalue[jj];
				} else if (1 > sscanf(tmp_value, "%lf", &pt->table[jj])) {
					error(1, "First data for parameter '%s' corrupted, it should be a type float or 'KIM'.\n",
                    FreeParamSet.name[k]);
				}
			} else {
				error(1, "Data for parameter '%s' corrupted in file '%s'.\n", FreeParamSet.name[k], filename);
			}
		  pt->idx[jj] = jj;
			jj++;
		}  
	}

	/* set end (end is cutoff) */
  double* pcutoff;
  
  pcutoff = KIM_API_get_data(pkim, "cutoff", &status);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", status);
    exit(1);
  }
    pt->end[0] = *pcutoff;


	/* free memory */
	free_model_object(&pkim);
  
	printf("Successfully read potential parameters that will be optimized.\n\n");
	return;
}


#endif



