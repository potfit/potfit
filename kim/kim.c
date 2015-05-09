/*
	 Functions to connect to KIM.

	Contributor: Mingjian Wen
*/


#include "../potfit.h"				/* to use `nconf'	*/
#include "kim.h"				/* to use `NeighObjectType' */


void InitKIM() {

	/* local vars */
	int i;

/*NOTE(remove) the following two lines to create EAM model are in util now*/
	/* create EAM temp model in current dir */
/*	CreateModel();
*/
	/* system call to make to model just made */
/*	MakeModel();
*/

	printf("\nInitializing KIM ... started\n");
	/* create KIM objects and do the necessary initialization */
	InitObject();

	/*NOTE(change), repeated work, publish parameter for each config */
	OptParamAllConfig = (OptParamType*) malloc(nconf*sizeof(OptParamType));

	for (i = 0; i <nconf; i++) {
	/* nest optimizable parameters */
		get_OptimizableParamInfo(pkimObj[i], &OptParamAllConfig[i]);	
		nest_OptimizableParamValue(pkimObj[i], &OptParamAllConfig[i], 1);
	/* Publish cutoff (cutoff only need to be published once, so here) */
		PublishCutoff(pkimObj[i], rcutmax);
	}

	printf("Initializing KIM ... done\n");
	fflush(stdout);
}


/***************************************************************************
*
* Initialize KIM objects
*
*	create them once for all 
***************************************************************************/

int InitObject() {

	/* local variables */
	int status;
	int numOfconf = nconf;		/* number of configurations in reference data */
  int i;

	/* Allocate memory for KIM objects */
	pkimObj = (void**) malloc(numOfconf * sizeof(void *));
	if (NULL == pkimObj) {
	  KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}
		

	/* Checking whether .kim files in test and model are compatible or not */
	for (i = 0; i < numOfconf; i++) {  				
		status = KIM_API_file_init(&pkimObj[i], "descriptor.kim", kim_model_name);
 		if (KIM_STATUS_OK > status)
 		{
   		KIM_API_report_error(__LINE__, __FILE__, "KIM_API_file_init", status);
   		exit(1);
  	}
	}

	/* QUESTION: does each config has the same number of species? e.g. there are
	two types of species, but a specific config may one have one species). If
	not, then ntypes below need to be modified from config to config */

	/*  Initialize KIM objects */
	for (i = 0; i < numOfconf; i++) { 
		status = CreateKIMObj(pkimObj[i], inconf[i], ntypes, cnfstart[i]);
		if (KIM_STATUS_OK > status)	{
   		KIM_API_report_error(__LINE__, __FILE__,
														"KIM: initializing objects failed", status);
   		exit(1);
  	}
	}
	return 0;		
}


/***************************************************************************
*
* Get optimizable parameters info
*
*	nest all optimizable parameters together (including PARAM_FREE_cutoff)
* parameter names are nested into name
* parameter values pinter are nested into value
* parameter ranks are nested into rank  
* parameter shapes ...

* PARAM_FREE_cutoff will be in the last slot of ParamType->name. We do this
* because potfit will include cutoff as an optimizable paramter if _cp (smooth
* cutoff) is enabled, and the cutoff will be in the last slot in potfit's
* potential table. Here, we put it in the last slot for the convenience of
* later usage. 
***************************************************************************/
int get_OptimizableParamInfo(void* pkim, OptParamType* OptParam) {

	/*local vars*/
	int status;
	int NumFreeParam;
	int maxStringLength;
	char* pstr;
	char buffer[128];
  char name[64];
	char type[16];
	char* tmp_FreeParamName;
	int tmp_NumFreeParam;
	int tmp_size;
	int OptParamNum = 0;
	int NumFreeParamNoDouble = 0;
	int i, j, k;
  int have_PARAM_FREE_cutoff = 0;

/* get the number of free parameters */
	status = KIM_API_get_num_free_params(pkim, &NumFreeParam, &maxStringLength);
	if (KIM_STATUS_OK > status) {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_num_free_params", status);
  	return(status);
	}

	/* get the descriptor file, pointed by pstr, the type will be phrased from pstr */	
	status = KIM_API_get_model_kim_str(kim_model_name, &pstr);
	if (KIM_STATUS_OK > status) {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_kim_str", status);
  	return(status);
	}

	/* initialize name */
	OptParam->name = (char**) malloc(sizeof(char*));
	if (NULL==OptParam->name) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}

	/* infinite loop to find PARAM_FREE_* of type double, only they are optimizable. 
		 Note, ``PARAM_FREE_cutoff'' will be appended at the end */
	/* It's safe to do pstr = strstr(pstr+1,"PARAM_FREE") because the ``PARAM_FREE''
		will never ever occur at the beginning of the ``descriptor.kim'' file */	
	while (1) {
  	pstr = strstr(pstr+1,"PARAM_FREE");
		if (pstr == NULL) {
			break;
		} else {
			snprintf(buffer, sizeof(buffer), "%s", pstr);
			sscanf(buffer, "%s%s", name, type);
			if (strcmp(name, "PARAM_FREE_cutoff") == 0) {
        have_PARAM_FREE_cutoff = 1;
        continue;
      } else if (strcmp(type, "double") == 0) {
				OptParamNum++;					
				if (OptParamNum > 1) {
					OptParam->name = (char**) realloc(OptParam->name, (OptParamNum)*sizeof(char*));
					if (NULL==OptParam->name) {
						KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
						exit(1);
					}
				}
				OptParam->name[OptParamNum - 1] =  /*maxStringLength+1 to hold the `\0' at end*/
															(char*) malloc((maxStringLength+1)*sizeof(char));
				if (NULL==OptParam->name[OptParamNum - 1]) {
					KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
					exit(1);
				}		       			
				strcpy(OptParam->name[OptParamNum - 1], name);
			} else {
				NumFreeParamNoDouble++;
			}
		}
	}
  /*append `PARAM_FREE_cutoff' at the end */
  if (have_PARAM_FREE_cutoff) {
		OptParamNum++;					
		OptParam->name = (char**) realloc(OptParam->name, (OptParamNum)*sizeof(char*));
		if (NULL==OptParam->name) {
			KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
      exit(1);
    }
		OptParam->name[OptParamNum - 1] =  /*maxStringLength+1 to hold the `\0' at end*/
    														(char*) malloc((maxStringLength+1)*sizeof(char));
		if (NULL==OptParam->name[OptParamNum - 1]) {
			KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
			exit(1);
		}		       			
		strcpy(OptParam->name[OptParamNum - 1], "PARAM_FREE_cutoff");
  }


	OptParam->Nparam = OptParamNum;
	

	/* allocate memory for value */
	OptParam->value = (double**) malloc(OptParam->Nparam * sizeof(double*));
	if (NULL==OptParam->value) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}	

	/* get the pointer to parameter */
	for(i = 0; i < OptParam->Nparam; i++ ) {
 		OptParam->value[i] = KIM_API_get_data(pkim, OptParam->name[i], &status);
	}
	if (KIM_STATUS_OK > status) {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", status);
    return(status);
  }

  /* allocate memory for rank */
	OptParam->rank = (int*) malloc(OptParam->Nparam * sizeof(int));
	if (NULL==OptParam->rank) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}	

	/* get rank */
	for(i = 0; i < OptParam->Nparam; i++) {
		OptParam->rank[i] = KIM_API_get_rank(pkim, OptParam->name[i], &status);
	}
	if (KIM_STATUS_OK > status) {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_rank", status);
    return(status);
  }

	/* allocate memory for shape */
	OptParam->shape = (int**) malloc(OptParam->Nparam * sizeof(int*));
	if (NULL==OptParam->shape) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}	
	for (i = 0; i < OptParam->Nparam; i++) {
 		OptParam->shape[i] = (int*) malloc(OptParam->rank[i] * sizeof(int));
		/* should not do the check, since rank may be zero. Then in some implementation
			malloc zero would return NULL */	
		/*	if (NULL==OptParam->shape[i]) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
		}*/	
	}

	/* get shape */
	for(i = 0; i < OptParam->Nparam; i++) {
		KIM_API_get_shape(pkim, OptParam->name[i], OptParam->shape[i], &status);
	}
	if (KIM_STATUS_OK > status) {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_shape", status);
    return(status);
  }
	

	/* nestedvalue is not allocated here, give NULL pointer to it */
	OptParam->nestedvalue = NULL;



	return KIM_STATUS_OK;
}




/*****************************************************************************
* nest optimizable parameters values
* nest the values (pointer) obtained from get_OptParamInfo() to nestedvalue.
* nestedvalue would be equal to value if all parameters have rank 0. 
* include_cutoff: flag, whether cutoff will be incldued in the nested list
* return: the length of nestedvalue

* After calling get_OptParamInfo, the PARAM_FREE_cutoff will be in the last slot
* of name list. 
* this function will nest the FREE_PARAM_cutoff into nestedvaule if
* inlcude_cutoff is true, otherwise, do not incldue it. 
*****************************************************************************/
nest_OptimizableParamValue(void* pkim, OptParamType* OptParam, int include_cutoff)
{
  /* local variables */
  int tmp_size;
  int total_size;       /* the size of the nested values*/
  int i, j, k;
  
  /* nest values  */
  /*first compute the total size of the optimizable parameters */
  total_size = 0;
	for (i = 0; i < OptParam->Nparam; i++ ) {	
	  if (strcmp(OptParam->name[i], "PARAM_FREE_cutoff") == 0 && !include_cutoff) {
      continue;
    }
    tmp_size = 1;
	  if (OptParam->rank[i] == 0) {
	    total_size += 1;
	  } else {
	    for (j = 0; j <OptParam->rank[i]; j++) {
	      tmp_size *= OptParam->shape[i][j];
	    }
	    total_size += tmp_size;
	  }
	}
	
  /* allocate memory for nestedvalue*/
  OptParam->nestedvalue = (double**) malloc((total_size) * sizeof(double*));

  /*copy the values (pointers) to nestedvalue */
  k = 0;
  for (i = 0; i < OptParam->Nparam; i++ ) {
	  if (strcmp(OptParam->name[i], "PARAM_FREE_cutoff") == 0 && !include_cutoff) {
      continue;
    }
    tmp_size = 1;	
    if (OptParam->rank[i] == 0) {
      OptParam->nestedvalue[k] = OptParam->value[i];
      k++;
    } else {
      for (j = 0; j < OptParam->rank[i]; j++) {
        tmp_size *= OptParam->shape[i][j];
      }
      for (j = 0; j <tmp_size; j++) {
        OptParam->nestedvalue[k] = OptParam->value[i]+j;
        k++;
      }
    }
  }

  return total_size;
}














































/***************************************************************************
*
*	Function used to publsih cutoff
*
***************************************************************************/
PublishCutoff(void* pkim, double cutoff) {

  /* local variable */
  int status;
  double* pcutoff;
  
 	pcutoff = KIM_API_get_data(pkim, "PARAM_FREE_cutoff", &status);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", status);
    return(status);
  }

  *pcutoff = cutoff;
  
  return 0;
}




/* added */	
/***************************************************************************
* 
* Create KIM object	 
*
* transferred varialbes:
* pkim: KIM object pointer
* Natoms: number of atoms in this configuration  
* Nspecies: number of atom species in this configuraiton 
*	start: index of the first atom of this configuration in potfit atom array
*
* golbal potfit varialbes, used but not transferred: 
* atoms 
* elements
* box_side_len  (this is defined by us, not potfit) 
***************************************************************************/

int CreateKIMObj(void* pkim, int Natoms, int Nspecies, int start)
{
  /* model inputs */
  int* numberOfParticles;
  int* numberOfSpecies;
  int* particleSpecies;
  double* coords;	
	int* numberContrib;
	NeighObjectType* NeighObject;
	const char* NBCstr;
	int NBC;
  /* local variables */
  int status;	
	int species_code;	
	int halfflag;
	int	Ncontrib;	/* number of contriburting particles */	
	int i, j, k;
	/* We have to allocate additional memory for neighbor list, since potfit 
		assigns the neighbors of each particle to atoms[i].neigh[k].nr (nr is
		the index neighbor k of atom i.). The memory is not continuous.
		The following variable (neighListLength), together with `BeginIdx' in
		the NeighObject, are used to gather the neighbor info to continuous memory. 
	*/	
	int neighListLength; /* total length of neighList */


	/* set value */
	Ncontrib = Natoms;

	/* allocate memory for NeighObject ( simulator should be in charge of
	allocating memory for neighbor objects, not kim model ) */
	NeighObject = (NeighObjectType*) malloc(sizeof(NeighObjectType));
  if (NULL == NeighObject) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}	

	/* Allocate memory via the KIM system */
  KIM_API_allocate(pkim, Natoms, Nspecies, &status);
  if (KIM_STATUS_OK > status)
  {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_allocate", status);
    return(status);
  }

	/* register for neighObject */
	KIM_API_setm_data(pkim, &status, 1*4,
				   					 "neighObject",     1,   NeighObject,   1);
  if (KIM_STATUS_OK > status) {
   KIM_API_report_error(__LINE__, __FILE__,"KIM_API_setm_data",status);
    return(status);		
	 }

	/* register for get_neigh */
  status = KIM_API_set_method(pkim, "get_neigh", 1, (func_ptr) &get_neigh);
  if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__,"KIM_API_set_method",status);
    return(status);				
	}	
	/* call Model's init routine */
  status = KIM_API_model_init(pkim);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_init", status);
    return(status);		
  }

/* Determine which neighbor list type to use */
  halfflag = KIM_API_is_half_neighbors(pkim, &status);
  if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__,"is_half_neighbors", status);
    return(status);		
	}

  /* Unpack data from KIM object */
  KIM_API_getm_data(pkim, &status, 5*3,
                    "numberOfParticles",   &numberOfParticles,   1,
                    "numberOfSpecies",     &numberOfSpecies,     1,
                    "particleSpecies",     &particleSpecies,     1,
                    "coordinates",         &coords,              1,
					 "numberContributingParticles",  &numberContrib, (1==halfflag) );
  if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);
  	return status;
	}

  /* Set values */
  *numberOfParticles = Natoms;
  *numberOfSpecies 	 = Nspecies;	
	*numberContrib		 = Ncontrib;


	/* set boxSideLengths if MI_OPBC is used */
	/* determine which NBC is used */
  status = KIM_API_get_NBC_method(pkim, &NBCstr);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method", status);
    return status; }
  if ((!strcmp("NEIGH_RVEC_H",NBCstr)) || (!strcmp("NEIGH_RVEC_F",NBCstr))) {
    NBC = 0; }
  else if ((!strcmp("NEIGH_PURE_H",NBCstr)) || (!strcmp("NEIGH_PURE_F",NBCstr))) {
    NBC = 1; }
  else if ((!strcmp("MI_OPBC_H",NBCstr)) || (!strcmp("MI_OPBC_F",NBCstr))) {
    NBC = 2; }
  else if (!strcmp("CLUSTER",NBCstr)) {
    NBC = 3; }
  else {
    status = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "Unknown NBC method", status);
    return status; }

	if (NBC == 2) {
		/* define local varialbe */
		double* boxSideLen;
		int which_conf;		/* which config we are in? */
	
		which_conf = atoms[start].conf; 

  	/* Unpack data from KIM object */
  	KIM_API_getm_data(pkim, &status, 1*3,
					 					"boxSideLengths",      &boxSideLen,			1 );
  	if (KIM_STATUS_OK > status)
 	  {
			KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);
 		 	return status;
		}
 
	  /* Set values */
		boxSideLen[0]	= box_side_len[DIM*which_conf + 0];
		boxSideLen[1]	= box_side_len[DIM*which_conf + 1];
		boxSideLen[2] = box_side_len[DIM*which_conf + 2];
	}		


  /* set the species types */
	/* to use this, the #C Header line in configuration file has to be included.
		Also the order of elements name after #C should be in accordane	with the
		species code in the first column of configuration data, ranging from 0 to
		(ntypes-1).  
		Here, assume potfit ensures this.
		e.g. if we have `#C Ar Ne', then the code for Ar and Ne should be 0 and 1,
		respectively. */
		
	/* QUESTION: in KIM_API_get_species_code, doest the second argument has to
		 to be totally the same as that kim descriptor file? e.g. Upper case or lower
		 case matters or not? */

	for (i = 0; i < *numberOfParticles; i++) { 
		/* in potfit, atom types range from 0 to (ntype-1) */
		if (atoms[start+i].type < 0 || atoms[start+i].type >= *numberOfSpecies) {
 			status = KIM_STATUS_FAIL;			
			KIM_API_report_error(__LINE__, __FILE__,
                           "Unexpected species detected", status);
  		return status;
		}
		else {
			j = atoms[start+i].type;
			species_code =  KIM_API_get_species_code(pkim, elements[j], &status);	
		  if (KIM_STATUS_OK > status)
  		{
   			KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_species_code," 
														"the species names need to be exactly the same"
														"as that in KIM standard file", status);
 			 	return status;
 			}
 	    particleSpecies[i] = species_code;
		}
	}

	/* set coords values */
	for (i = 0; i < *numberOfParticles; i++) {
		coords[DIM*i]   = atoms[start+i].pos.x;
		coords[DIM*i+1] = atoms[start+i].pos.y;
		coords[DIM*i+2] = atoms[start+i].pos.z;
	}

	/* calcualte the length of neighborList */
	neighListLength = 0;
	for (i = 0; i < *numberOfParticles; i++) {
		neighListLength += atoms[start+i].num_neigh;	
	}

  /* allocate memory for NeighObject members */
  NeighObject->NNeighbors = (int*) malloc((*numberOfParticles) * sizeof(int));
  if (NULL==NeighObject->NNeighbors) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}	
  NeighObject->neighborList = (int*) malloc(neighListLength * sizeof(int));
  if (NULL==NeighObject->neighborList) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);		
	}	
 NeighObject->RijList = (double*) malloc((DIM*neighListLength) * sizeof(double));
  if (NULL==NeighObject->RijList) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);		
	}	
	NeighObject->BeginIdx	= (int*) malloc((*numberOfParticles) * sizeof(int));
	if (NULL==NeighObject->BeginIdx) {
		KIM_API_report_error(__LINE__,__FILE__,"malloc unsuccessful", -1);
		exit(1);		
	}

	/* copy the number of neighbors to NeighObject.NNeighbors */
	status = KIM_STATUS_FAIL; /* assume staus fails */
	for (i = 0; i < *numberOfParticles; i++) {
		NeighObject->NNeighbors[i] = atoms[start+i].num_neigh;	
	}
	if (i == *numberOfParticles){
		status = KIM_STATUS_OK;
	}
  if (KIM_STATUS_OK > status)	{
   	KIM_API_report_error(__LINE__, __FILE__,
												"copy number of neighbors to NNeighbors failed", status);
  return status;
 	}

	/* copy neighborlist from distributed memory locations to 
		continuous ones.
		copy atoms[i].neigh[j].nr to NeighObject.neighborList,
		copy atoms[i].neigh[k].dist.x ... to 	NeighObject.Rij[DIM*k] ... */	
	status = KIM_STATUS_FAIL; /* assume staus fails */	
	k = 0;
	for (i = 0; i < *numberOfParticles; i++) {
		NeighObject->BeginIdx[i] = k;
		for (j = 0; j < NeighObject->NNeighbors[i]; j++) {
			NeighObject->neighborList[k]  =	atoms[start+i].neigh[j].nr - start;
			NeighObject->RijList[DIM*k]   =	atoms[start+i].neigh[j].dist.x;
			NeighObject->RijList[DIM*k+1] =	atoms[start+i].neigh[j].dist.y;
			NeighObject->RijList[DIM*k+2] =	atoms[start+i].neigh[j].dist.z;
			k++;

/* unittest */
/* used together with to see that the stuff in the neighbor lista are correct
	set  start == 0 , will check for the first config. We set j==0, and j==last
	atom in the neighbor list of an atom, since we don't want that verbose info.
*/
/*
if (start != 0 && (j== 0 || j == atoms[start+i].num_neigh-1 )) {
printf("last neighbor: %d\n", atoms[start+ *numberOfParticles-1].num_neigh);
printf("which atom: %d\n",i);
printf("%d %f %f %f\n", atoms[start+i].neigh[j].nr,
												atoms[start+i].neigh[j].dist.x,
												atoms[start+i].neigh[j].dist.y,
												atoms[start+i].neigh[j].dist.z );
}
*/
/* unittest ends*/
		}
	}
	
	if (i == *numberOfParticles && k == neighListLength){
		status = KIM_STATUS_OK;
	}
  if (KIM_STATUS_OK > status)	{
   	KIM_API_report_error(__LINE__, __FILE__,
												"copy neighbor list failed", status);
  	return status; } 

	/* If the number of neighbors of an atom is zero, set the BeginIdx to the 
		 last position in neghbor list. Actually, the main purpose is to ensure 
		 that the BeginIdx of the last atom will not go beyond limit of 
		 neighborlist length.
		 e.g. there are 128 atoms in the config, and we use half neighbor list,
		 then the 128th atom will have no neighbors. Then the the begin index
		 for the last atom, BeginIdx[127] will go beyond the limit of Neighbor
		 list, which may result in segfault. So we need to do something to avoid
		 this.
		 I'm sure, there are better ways to do this.
	*/	
	for (i = 0; i < *numberOfParticles; i++) {
		if( NeighObject->NNeighbors[i] == 0) {
			NeighObject->BeginIdx[i] = k-1;
		}
	}

 	return KIM_STATUS_OK;
}
/* added ends */


/***************************************************************************
* 
* get_neigh 
*
***************************************************************************/

int get_neigh(void* kimmdl, int *mode, int *request, int* part,
                       int* numnei, int** nei1part, double** Rij)
{
   /* local variables */
  intptr_t* pkim = *((intptr_t**) kimmdl);
  int partToReturn;
  int status;
  int* numberOfParticles;
	int idx; /* index of first neighbor of each particle*/
  NeighObjectType* nl;

  /* initialize numnei */
  *numnei = 0;

	/* unpack neighbor list object */
  numberOfParticles = (int*) KIM_API_get_data(pkim, "numberOfParticles", &status);
  if (KIM_STATUS_OK > status) {
  KIM_API_report_error(__LINE__, __FILE__,"get_data", status);
	}

  nl = (NeighObjectType*) KIM_API_get_data(pkim, "neighObject", &status);
  if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__,"get_data", status);
	}

  /* check mode and request */
  if (0 == *mode) /* iterator mode */
  {
		if (0 == *request) /* reset iterator */
    {
	    (*nl).iteratorId = -1;
      return KIM_STATUS_NEIGH_ITER_INIT_OK;
		}
    else if (1 == *request) /* increment iterator */
    {
  	  (*nl).iteratorId++;
      if ((*nl).iteratorId >= *numberOfParticles)
      {
				return KIM_STATUS_NEIGH_ITER_PAST_END;
      }
      else
      {
				partToReturn = (*nl).iteratorId;
      }
    }
    else /* invalid request value */
    {
 	    KIM_API_report_error(__LINE__, __FILE__,"Invalid request in get_periodic_neigh",
													 KIM_STATUS_NEIGH_INVALID_REQUEST);
      return KIM_STATUS_NEIGH_INVALID_REQUEST;
		}
 	}
  else if (1 == *mode) /* locator mode */
  {
    if ((*request >= *numberOfParticles) || (*request < 0)) /* invalid id */
    {
			KIM_API_report_error(__LINE__, __FILE__,"Invalid part ID in get_periodic_neigh",
	 												KIM_STATUS_PARTICLE_INVALID_ID);
			return KIM_STATUS_PARTICLE_INVALID_ID;
    }
    else
    {
      partToReturn = *request;
    }
  }
  else /* invalid mode */
  {
		KIM_API_report_error(__LINE__, __FILE__,"Invalid mode in get_periodic_neigh",
		 											KIM_STATUS_NEIGH_INVALID_MODE);
    return KIM_STATUS_NEIGH_INVALID_MODE;
  }

	/* index of the first neigh of each particle */
	idx = (*nl).BeginIdx[partToReturn];

  /* set the returned part */
  *part = partToReturn;

  /* set the returned number of neighbors for the returned part */
  *numnei = (*nl).NNeighbors[partToReturn];

  /* set the location for the returned neighbor list */
  *nei1part = &((*nl).neighborList[idx]);

  /* set the pointer to Rij to appropriate value */
  *Rij = &((*nl).RijList[DIM*idx]);

	return KIM_STATUS_OK;
}

/* added ends */





/* added */
/***************************************************************************
* 
* Calculate force from KIM (general force, including force, virial, energy)
*
* transferred variaves 
*	Pkim: KIM object
* the following three are also model output 
*	energy;
* force;
*	virial;
*
***************************************************************************/

int CalcForce(void* pkim, double** energy, double** force, double** virial,
							int useforce, int usestress)
{	
	/* local variables */
  int status;
	/* Access to free parameters, also unpack other variables as needed */
	/* potfit will always compute energy, so is here */
	KIM_API_getm_data(pkim, &status, 1*3,
                    "energy",              energy,              1);
	if (useforce)
		KIM_API_getm_data(pkim, &status, 1*3,
                    "forces",              force,               1);
	if(usestress)
		KIM_API_getm_data(pkim, &status, 1*3,
										"virial",           	 virial,       			 1 );


  if (KIM_STATUS_OK > status)
  {
		KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);
    return status;
  }

  /* Call model compute */
  status = KIM_API_model_compute(pkim);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_compute", status);
    return status;
  }

	return KIM_STATUS_OK;
}
/* added ends */










/******************************************************************************
* get_OptimizableParam
* inquire KIM to get the number of optimizable parameters 
* and their names, which will be used later to nest the ParameterList


*******************************************************************************/

int get_OptimizableParamSize(OptParamType* OptParamSize, int include_cutoff) 
{
  /*local variables */
	void* pkim;
	int status;
	int size;			/* This would be equal to Nparam if all parameters have rank zero.*/	
	int NumFreeParamNoDouble = 0;  /*number of FREE_PARAM_* with type other than double*/
	char* pstr;
	char buffer[128];
	char name[64];
	char type[16];
	int tmp_size;
	int i, j;

	/* create a temporary KIM objects, in order to inquire KIM model for 
	 PARAM_FREE_* parameters info, delete at the end of the function */
	status = KIM_API_file_init(&pkim, "descriptor.kim", kim_model_name);
	if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__, "KIM_API_file_init", status);
   	exit(1);
  }
	
	/* Allocate memory via the KIM system */
  KIM_API_allocate(pkim, 2, 1, &status);
  if (KIM_STATUS_OK > status)
  {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_allocate", status);
    return(status);
  }

	/* call Model's init routine */
  status = KIM_API_model_init(pkim);
  if (KIM_STATUS_OK > status)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_init", status);
    return(status);		
  }
	
	/*get the info of the Optimizable parameters */
	get_OptimizableParamInfo(pkim, OptParamSize);

	/* compute the total size of the optimizable */
  size = nest_OptimizableParamValue(pkim, OptParamSize, include_cutoff);


	/*number of free parameters of type other than double */
		/* get the descriptor file, pointed by pstr, the type will be phrased from pstr */	
	status = KIM_API_get_model_kim_str(kim_model_name, &pstr);
	if (KIM_STATUS_OK > status) {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_kim_str", status);
  	return(status);
	}
		/* infinite loop to find PARAM_FREE_* of type other than double*/	
	while (1) {
  	pstr = strstr(pstr+1,"PARAM_FREE");
		if (pstr == NULL) {
			break;
		} else {
			snprintf(buffer, sizeof(buffer), "%s", pstr);
			sscanf(buffer, "%s%s", name, type);
			if ( strcmp(type, "double") != 0) {
				NumFreeParamNoDouble++;
			}
		}
	}

	printf("\nKIM Model being used: %s.\n\n", kim_model_name);
	if (NumFreeParamNoDouble != 0) {
    printf("-  There is (are) %d `PARAM_FREE_*' parameter(s) of type other than "
						"`double'.\n", NumFreeParamNoDouble);
	}


	/* deallocate */  
/*	
	status = KIM_API_model_destroy(pkim);
  if (KIM_STATUS_OK > status){ 
		KIM_API_report_error(__LINE__, __FILE__,"destroy", status);
		return status;
	}
*/

/*
	KIM_API_free(pkim, &status);
  if (KIM_STATUS_OK > status){ 
		KIM_API_report_error(__LINE__, __FILE__,"destroy", status);
		return status;
	}
*/

	/* return value */
	return size;
}

















































































































































































































































































/***************************************************************************
* 
* Create parameter file EAM
*	
* Number of R data, deltaR, number of rho data, deltaRho and cutoff will be 
* write here. Potential data will be filled with 0.0, and they will be altered
* later by calling the PublishParam function. 
*
***************************************************************************/
int CreateParamFileEAM() {	

	/* local varialbes */
	char tmpstring[80];
  FILE *pFile;
	int NumRho;	/* number of rho data for embedding function */
	int NumR;		/* nubmer of r data for density function and pair function */
		/* IdxRho IdxR which potfit potential we are visiting? ranging from 0 to ncols.
			 The data are arranged as \Phi(r),..,\rho(r),...,U(\rho) in opt_pot.table.
			 For example, if there are 2 species i j, the data will be arranged as 
			 \Phi(ii), \Phi(ij), \Phi(jj), \rho(i), \rho(j) U(i),U(j) */
	int IdxRho; 
	int IdxR;
	int i,j;

	sprintf(tmpstring, "%s/%s.params", kim_model_name, kim_model_name);	
	pFile = fopen(tmpstring,"w");
	if(pFile == NULL)
		error(1, "Error creating file: KIM Makefile. %s %d\n", __FILE__, __LINE__);
	else {
	
  	/* 3 lines of comments */
		fprintf(pFile, "# EAM potential in Dynamo `setfl' format, generated by potfit-kim.\n");		
		fprintf(pFile, "# For more info about KIM, visit: https://openkim.org\n");				
		fprintf(pFile, "# For more info about potfit, visit: http://potfit.sourceforge.net\n");		
		
    /* number of atoms types and atom species*/
		fprintf(pFile, "    %d", ntypes);	
		for (i = 0; i < ntypes; i++) {
			fprintf(pFile, "  %s", elements[i]);
		}
		fprintf(pFile, "\n");

		/* number of data points, delta  rcut */
		IdxRho = ntypes*(ntypes + 1)/2 + ntypes;	
	  NumRho = opt_pot.last[IdxRho] - opt_pot.first[IdxRho] + 1;
		IdxR = 0;
		NumR = opt_pot.last[IdxR] - opt_pot.first[IdxR] + 1;
		fprintf(pFile, "  %d %22.15e %d %22.15e %22.15e\n", NumRho, opt_pot.step[IdxRho], 
					  NumR, opt_pot.step[IdxR], rcutmax);

		/* embedding data, density data and their header */
		for (i = 0; i < ntypes; i++) {

/*NOTE(change)*/
			/* the header */
			fprintf(pFile, "  0  0.0  0.0  lattice_type\n");

			/* fill embedding data with 0.0 */
			for (j = 0; j < NumRho; j++) {
				fprintf(pFile, "  %2.1f", 0.0);
			}
			fprintf(pFile, "\n");
			/* fill density data with 0.0 */
			for (j = 0; j < NumR; j++) {
				fprintf(pFile, "  %2.1f", 0.0);
			}
			fprintf(pFile, "\n");
		}

		/* pair data */
		for (i = 0; i < ntypes*(ntypes + 1)/2; i++) {	
			IdxR = i;
			NumR = opt_pot.last[IdxR] - opt_pot.first[IdxR] + 1;
			/* fill pair data with 0.0 */
			for (j = 0; j < NumR; j++) {
				fprintf(pFile, "  %2.1f", 0.0);
			}
			fprintf(pFile, "\n");
		}
	}	
	fflush(pFile);
	fclose(pFile);
	return 0;
}


/***************************************************************************
*
* Create KIM Model
*
*	Create a directory for KIM model, and then create `Makefile' and '*.params'. 
***************************************************************************/

int CreateModel()
{
 
  char kim_model_driver_name[] = "EAM_CubicCompleteSpline__MD_000000111111_000" ;
	/*char kim_model_driver_name[] = "EAM_Dynamo__MD_120291908751_001";
*/
	/* local variables */
	char SysCmd[256];
	char tmpstring[80];
	int i,j;
  FILE *pFile;

  /* whether kim-api-build-config exists has been tested in the Makefile*/
	/* create Makefile.KIM_Config */
	system("kim-api-build-config --makefile-kim-config > Makefile.KIM_Config");

	/* create model directoy */
	/* if the model file name too long? */
	if(strlen(kim_model_name) > 249)
		error(1, "KIM Model name: '%s' is too long, use a shorter one.\n", kim_model_name);

/*NOTE(change) it would be better check whether the dir exits or not */
	strcpy(SysCmd, "mkdir ");
	strcat(SysCmd, kim_model_name);
	system(SysCmd);

	/* create Makefile */
	sprintf(tmpstring, "%s/Makefile", kim_model_name);
	pFile = fopen(tmpstring,"w");
	if(pFile == NULL)
		error(1, "Error creating file: KIM Makefile. %s %d\n", __FILE__, __LINE__);
	else {
		fprintf(pFile,
			"#\n"
			"# Copyright (c) 2013--2014, Regents of the University of Minnesota.\n"
			"# All rights reserved.\n"
			"#\n"
			"# Contributors:\n"
			"#    Ryan S. Elliott\n"
			"#    Ellad B. Tadmor\n"
			"#    Mingjian Wen\n"
			"#\n\n\n" );	

		fprintf(pFile,"# load all basic KIM make configuration\n" );	
		fprintf(pFile,"ifeq ($(wildcard ../Makefile.KIM_Config),)\n"
	  	"  $(error ../Makefile.KIM_Config does not exist. Something is wrong with your KIM API package setup)\n"
			"endif\n"
			"include ../Makefile.KIM_Config\n\n" );	

		fprintf(pFile,"# set model driver specific details\n" );
		fprintf(pFile,"MODEL_DRIVER_NAME    := %s\n", kim_model_driver_name);	
		fprintf(pFile,"MODEL_NAME           := %s\n", kim_model_name);	

		for (i = 0; i < ntypes; i++) {
			sprintf(tmpstring, "SPECIES_%03d_NAME", i+1);		
			fprintf(pFile,"%s     := %s\n", tmpstring, elements[i]);	
		}
		sprintf(tmpstring, "%s.params", kim_model_name);			
		fprintf(pFile,"PARAM_FILE_001_NAME  := %s\n\n", tmpstring);	
	
		fprintf(pFile,
			"# APPEND to compiler option flag lists\n"
			"#FFLAGS   +=\n"
			"#CFLAGS   +=\n"
			"#CXXFLAGS +=\n"
			"#LDFLAGS  +=\n"
			"#LDLIBS   +=\n\n");	

		fprintf(pFile, "# load remaining KIM make configuration\n"
		"include $(KIM_DIR)/$(builddir)/Makefile.ParameterizedModel\n");	
	}
	fflush(pFile);
	fclose(pFile);

	/* create EAM parameter file*/
#ifdef EAM
  CreateParamFileEAM();
#endif /* EAM*/

	return 0;
}


/***************************************************************************
*
* System call to make the model 
*
***************************************************************************/
int MakeModel() {

	/* local variables */
	char SysCmd[256];
	
	printf("\n\nMaking KIM model ... starting\n");
  sprintf(SysCmd,"cd %s; make", kim_model_name);
	system(SysCmd);
	printf("Making KIM model ... done\n");
	
	return 0;
}


