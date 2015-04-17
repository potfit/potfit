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
	/* create KIM objects and do the necessary initialization */
	InitObject();

	/*NOTE(change), repeated work, publish parameter for each config */
	FreeParamAllConfig = (FreeParamType*) malloc(nconf*sizeof(FreeParamType));

	for (i = 0; i <nconf; i++) {
	/* Create free parameters list */
		CreateFreeParamList(pkimObj[i], &FreeParamAllConfig[i]);	
	/* Publish cutoff (cutoff only need to be published once, so here) */
		PublishCutoff(pkimObj[i], &FreeParamAllConfig[i], rcutmax);
	}

	printf("\nInitializing KIM ... done\n");
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
		

get_OptimizableParam();

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
* Create free parameters list
*
*	nest all free parameters together
* parameter names are nested into name
* parameter values pinter are nested into value
* parameter ranks are nested into rank  
***************************************************************************/
int CreateFreeParamList(void* pkim, FreeParamType* FreeParam) {

	/*local vars */
	int status;
	int numberFreeParameters;
	int maxStringLength;
	int shape[3];
	int i,j;

	/* get the number of free parameters */
	status = KIM_API_get_num_free_params(pkim, &numberFreeParameters, &maxStringLength);
	if (KIM_STATUS_OK > status) {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_num_free_params", status);
    return(status);
  }
	FreeParam->Nparam = numberFreeParameters;


	/* allocate memory for name */
	FreeParam->name = (char**) malloc(numberFreeParameters * sizeof(char*));
  if (NULL==FreeParam->name) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}	
	for (i = 0; i < numberFreeParameters; i++) {
 		FreeParam->name[i] = (char*) malloc( maxStringLength * sizeof(char));
		if (NULL==FreeParam->name[i]) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
		}	
	}

	/* get name and record the position of cutoff in the list */
	for(i = 0; i < numberFreeParameters; i++ ) {
		status = KIM_API_get_free_parameter(pkim, i, &FreeParam->name[i]);
		if (!strcmp(FreeParam->name[i], "PARAM_FREE_cutoff")) {
			FreeParam->cutoffIdx = i;		
		}
	}
	if (KIM_STATUS_OK > status) {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_free_parameter", status);
    return(status);
  }

	/* allocate memory for value */
	FreeParam->value = (double**) malloc(numberFreeParameters * sizeof(double*));
	if (NULL==FreeParam->value) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}	

	/* get the pointer to parameter */
	for(i = 0; i < numberFreeParameters; i++ ) {
 		FreeParam->value[i] = KIM_API_get_data(pkim, FreeParam->name[i], &status);
	}
	if (KIM_STATUS_OK > status) {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", status);
    return(status);
  }

  /* allocate memory for rank */
	FreeParam->rank = (int*) malloc(numberFreeParameters * sizeof(int));
	if (NULL==FreeParam->rank) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}	

	/* get rank */
	for(i = 0; i < numberFreeParameters; i++) {
		FreeParam->rank[i] = KIM_API_get_rank(pkim, FreeParam->name[i], &status);
	}
	if (KIM_STATUS_OK > status) {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_rank", status);
    return(status);
  }

	/* allocate memory for shape */
	FreeParam->shape = (int**) malloc(numberFreeParameters * sizeof(int*));
	if (NULL==FreeParam->shape) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}	
	for (i = 0; i < numberFreeParameters; i++) {
 		FreeParam->shape[i] = (int*) malloc(FreeParam->rank[i] * sizeof(int));
		/* should not do the check, since rank may be zero. Then in some implementation
			malloc zero would return NULL */	
		/*	if (NULL==FreeParam->shape[i]) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
		}*/	
	}

	/* get shape */
	for(i = 0; i < numberFreeParameters; i++) {
		KIM_API_get_shape(pkim, FreeParam->name[i], FreeParam->shape[i], &status);
	}
	if (KIM_STATUS_OK > status) {
  	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_shape", status);
    return(status);
  }



	printf("\nThe following parameters should and will be pubished to KIM.\n"
					"(PARAM_FREE_cutoff is publsihed by PublishCutoff(), and the others\n"
					"are published by PublishParam().)\n");
	printf("      param name                  param shape\n");
	printf("     ############                ############\n");
	for(i = 0; i < numberFreeParameters; i++ ) {
		printf("  %-35s[ ", FreeParam->name[i]);
		for(j = 0; j < FreeParam->rank[i]; j++) {
			printf("%d ", FreeParam->shape[i][j]);
		}
		printf("]\n");
	}


	return 0;
}



/***************************************************************************
*
*	Function used to publsih cutoff
*
***************************************************************************/
PublishCutoff(void* pkim, FreeParamType* FreeParam, double cutoff) {
	*(FreeParam->value[FreeParam->cutoffIdx]) = cutoff;
	return 0;
}


/***************************************************************************
*
* Free KIM 
*
*	Free memory allocated by KIM. e.g., neighbor list, KIM object and KIM model
*
***************************************************************************/

void FreeKIM(void)
{
	/* local variables */
	int status;
	int i;
	NeighObjectType* NeighObject;
	void* pkim;

	for (i = 0; i < nconf; i++) {		
		pkim = pkimObj[i];
		
		/* call model destroy */
		status = KIM_API_model_destroy(pkim);
    if (KIM_STATUS_OK > status) 
			KIM_API_report_error(__LINE__, __FILE__,"destroy", status);

		/* free memory of neighbor list */
		/* first get neighbor object from KIM */
	  NeighObject = (NeighObjectType*) KIM_API_get_data(pkim, "neighObject", &status);
  	if (KIM_STATUS_OK > status) 
			KIM_API_report_error(__LINE__, __FILE__,"get_data", status);
		/* then free the stuff in neighbor list */	
		free(NeighObject->NNeighbors);	
		free(NeighObject->neighborList);	
		free(NeighObject->RijList);	
		free(NeighObject->BeginIdx);	
		/* free the neighbor list itself */
		free(NeighObject);	
	
		/* free memory of KIM objects */
  	KIM_API_free(pkim, &status);
 		if (KIM_STATUS_OK > status) 
			KIM_API_report_error(__LINE__, __FILE__,"destroy", status);
	}

  free(pkimObj);
	/* every thing is great */
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

* initialize the two golbal variabls: int OptParamNum; char** OptParamNames;

*******************************************************************************/

int get_OptimizableParam() {

	/*local variables */
	void* pkim;
	int status;
	int NumFreeParam;
	int maxStringLength;
	int OptParamSize = 0;			/* This would be equal to OptParamNum
																	if all parameters have rank zero.  */	
	int NumFreeParamNoDouble = 0;  /*number of FREE_PARAM_* with type other than double*/
	char* pstr;
	char buffer[128];
  char name[64];
	char type[16];
	char* tmp_FreeParamName;
	int tmp_NumFreeParam;
	int tmp_size;
	int* rank;
	int** shape;
	int i, j;
	/*initialize the global variable */
	OptParamNum = 0;

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


	/* infinite loop to find PARAM_FREE_* of type double, only they are optimizable. 
		 Note, ``PARAM_FREE_cutoff'' is not optimizable. */
	/* It's safe to do pstr = strstr(pstr+1,"PARAM_FREE") because the ``PARAM_FREE''
		will never ever occur at the beginning of the ``descriptor.kim'' file */	
	while (1) {
  	pstr = strstr(pstr+1,"PARAM_FREE");
		if (pstr == NULL) {
			break;
		} else {
			snprintf(buffer, sizeof(buffer), "%s", pstr);
			sscanf(buffer, "%s%s", name, type);

			if (strcmp(name, "PARAM_FREE_cutoff") != 0   /*cutoff is not a changable param*/
					&&				 strcmp(type, "double") == 0) {

				OptParamNum++;	

				OptParamNames = (char**) realloc(OptParamNames, (OptParamNum)*sizeof(char*));
				if (NULL==OptParamNames) {
					KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
					exit(1);
				}
	
				OptParamNames[OptParamNum - 1] =  /*maxStringLength+1 to hold the `\0' at end*/
															(char*) malloc((maxStringLength+1)*sizeof(char));
				if (NULL==OptParamNames[OptParamNum - 1]) {
					KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
					exit(1);
				}		       

				strcpy(OptParamNames[OptParamNum - 1], name);

			} else if ( strcmp(type, "double") != 0) {
				NumFreeParamNoDouble++;
			}
		}
	}

	/* Although the above should be correct, we double check it. This time, we 
		inquire OptParamNum through KIM call, not by phrasing the descriptor file. */
		/*first, check whether PARAM_FREE_cutoff is published */
	tmp_NumFreeParam = NumFreeParam;
	for(i = 0; i < NumFreeParam; i++ ) {
		status = KIM_API_get_free_parameter(pkim, i, &tmp_FreeParamName);
		if (KIM_STATUS_OK > status) {
  		KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_free_parameter", status);
    	return(status);
  	}
		if (!strcmp(tmp_FreeParamName, "PARAM_FREE_cutoff")) {
			tmp_NumFreeParam--;	
			break;
		}
	}
		/* then, check whether they match or not */
	if (tmp_NumFreeParam != (OptParamNum + NumFreeParamNoDouble)) {
		KIM_API_report_error(__LINE__, __FILE__, 
												"Number of optimizable parameters is incorrect", -1);
		exit(1);
	} 

	/* allocate memory for rank */
	rank = (int*) malloc(OptParamNum * sizeof(int));
	if (NULL==rank) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}	

	/* get rank */
	for(i = 0; i < OptParamNum; i++) {
		rank[i] = KIM_API_get_rank(pkim, OptParamNames[i], &status);
		if (KIM_STATUS_OK > status) {
		 	KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_rank", status);
	    return(status);
	  }
	}

	/* allocate memory for shape */
	shape = (int**) malloc(OptParamNum * sizeof(int*));
	if (NULL==shape) {
		KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
		exit(1);
	}	
	for (i = 0; i < OptParamNum; i++) {
		shape[i] = (int*) malloc(rank[i] * sizeof(int));
	/* should not do the check, since rank may be zero. Then in some implementation
		malloc zero would return NULL */	
	/*	if (NULL==shape[i]) {
			KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
			exit(1);
		}*/	
	}

	/* get shape */
	for(i = 0; i < OptParamNum; i++) {
		KIM_API_get_shape(pkim, OptParamNames[i], shape[i], &status);
		if (KIM_STATUS_OK > status) {
  		KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_shape", status);
    	return(status);
  	}	
	}

	/* compute the total size of the optimizable */
	for (i = 0; i < OptParamNum; i++ ) {	
		tmp_size = 1;
		if (rank[i] == 0) {
			OptParamSize += 1;
		} else {
			for (j = 0; j < rank[i]; j++) {
				tmp_size *= shape[i][j];
			}
			OptParamSize += tmp_size;
		}
	}



	printf("KIM Model being used: %s.\n", kim_model_name);
	printf("- There is (are) %d `PARAM_FREE_*' parameter(s) of type other than "
						"`double'.\n", NumFreeParamNoDouble);
	printf("- The following parameters should and will be pubished to KIM by "
					"`PublishParam()' each time `calc_force' is called.\n");
	printf("         param name                     param size\n");
	printf("        ############                   ############\n");
	for(i = 0; i < OptParamNum; i++ ) {
		printf("     %-35s[ ", OptParamNames[i]);
		for(j = 0; j < rank[i]; j++) {
			printf("%d ", shape[i][j]);
		}
		printf("]\n");
	}
	printf("- PARAM_FREE_cutoff is publsihed by `PublishCutoff()' once for all.\n");



	/* deallocate */  
	status = KIM_API_model_destroy(pkim);
  if (KIM_STATUS_OK > status){ 
		KIM_API_report_error(__LINE__, __FILE__,"destroy", status);
		return status;
	}
/*
	KIM_API_free(pkim, &status);
  if (KIM_STATUS_OK > status){ 
		KIM_API_report_error(__LINE__, __FILE__,"destroy", status);
		return status;
	}
*/
	free(shape);

	/* return value */
	return OptParamSize;
}






/* added */
/***************************************************************************
*
* Publish KIM parameter 
*
* transferred varialbes:
*	pkim: KIM ojbect
*	PotTable: potential table where the potential paramenters are stored
*	CutOff: cutoff radius   
*
***************************************************************************/
/*Don't forget to publish params for all KIM objest */
int PublishParam(void* pkim, FreeParamType* FreeParam, double* PotTable) {
	/*local variables*/  
	int status;
	double* param_epsilon;
	double* param_sigma;

	/* set convenient pointers to the parameters */
	param_epsilon = FreeParam->value[1];
	param_sigma   = FreeParam->value[2];
	
	/*publish parameters*/
 	*param_epsilon = PotTable[2]; 
	*param_sigma	 = PotTable[3];
 	status = KIM_API_model_reinit(pkim);
 	if (KIM_STATUS_OK > status)	
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_reinit", status);
    return status;
  }
	return KIM_STATUS_OK;
}
/* added ends */











































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


