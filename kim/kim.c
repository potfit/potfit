/*******************************************************************************
*
* kim.c 
*
* All functions all defined in this file
*
*******************************************************************************/

#define KIM_MAIN
#include "kim.h"
#undef KIM_MAIN

#include <stdio.h>
#include "../potfit.h" 


/*******************************************************************************
*
* Init KIM objects, each object for a reference configuration.
* Init optimizable parameters; nest them a single variable. 
*
*******************************************************************************/

void init_KIM() 
{
  int i;
  
  printf("\nInitializing KIM ... started\n");

  /* create KIM objects and do the necessary initialization */
  init_object();
  
  /* create free parameter data sturct and nest optimizable parameters */
  init_optimizable_param();

  for (i = 0; i <nconf; i++) {
    /* Publish cutoff (only needs to be published once, so here) */
    publish_cutoff(pkimObj[i], rcutmax);
  }

  printf("Initializing KIM ... done\n");
}


/*******************************************************************************
*
* Create KIM objects, and init the argument values
*
* potfit global variables:
* nconf
* inconf
* ntypes
* cnfstart
* useforce
* usestress
*
* potfit-KIM global variables:
* kim_model_name
*
*******************************************************************************/ 

void init_object()
{
  /* local variables */
  NeighObjectType* NeighObject;
  int status;
  int u_f;
  int u_s;
  int i;

  /* Allocate memory for KIM objects */
  pkimObj = (void**) malloc(nconf * sizeof(void *));
  if (NULL == pkimObj) {
    KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
    exit(1);
  }
	reg_for_free(pkimObj, "pkimObj");

  for (i = 0; i < nconf; i++) { 
     
    /* write descriptor file */
    u_f = useforce[i];
    u_s = 0;
#ifdef STRESS
    u_s = usestress[i];
#endif
    status = write_final_descriptor_file(u_f, u_s);
    if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__, "write_final_descriptor file", status);
      exit(1);
    }

    /* QUESTION: does each config has the same number of species? e.g. there are
     * two types of species, but a specific config may one have one species). If
     * not, then ntypes below need to be modified from config to config */
    /* Answer: actually this does not need any worry. If there is only one species in
     * the whole configuration, but we let ntypes = 2 in the following function call,
     * we just allocate more memory. But the species code for each atom would be
     * correct(see function init_KIM_API_argument). Then KIM Model would know how to
     * do the force calculation. */ 

    /* init KIM API object and allocate memory for data argument */
    status = setup_KIM_API_object(&pkimObj[i], inconf[i], ntypes, kim_model_name);
    if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__, "setup_KIM_API_object", status);
      exit(1);
    }

    /* init KIM API argument values */
    init_KIM_API_argument(pkimObj[i], inconf[i], ntypes, cnfstart[i]);
    if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__, "init_KIM_API_argument", status);
      exit(1);
    }

    /* allocate memory for NeighObject */ 
    NeighObject = (NeighObjectType*) malloc(sizeof(NeighObjectType));
    if (NULL == NeighObject) {
      KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
      exit(1);
    } 
		reg_for_free(NeighObject, "NeighObject");
    
		/* register access of neighborlist in KIM API */
    setup_neighborlist_KIM_access(pkimObj[i], NeighObject);
    if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__, "setup_neighborlist_KIM_access", status);
      exit(1);
    }

    /* initialize neighbor list */
    status = init_neighborlist(NeighObject, inconf[i], cnfstart[i]);
    if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__,"init_nieghborlist",status);
      exit(1);
    }
  }
}

/*******************************************************************************
*
* Create KIM API objects and allocate memory 
*
*******************************************************************************/ 

int setup_KIM_API_object(void** pkim, int Natoms, int Nspecies, char* modelname)
{
  /* local vars */
  int status;

  /* initialize KIM API object */
  status = KIM_API_file_init(pkim, "descriptor.kim", modelname);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_file_init", status);
    return(status);
  }

  /* Allocate memory for each data argument of initialized KIM object */
  KIM_API_allocate(*pkim, Natoms, Nspecies, &status);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_allocate", status);
    return(status);
  }

  /* call Model's init routine */
  status = KIM_API_model_init(*pkim);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_init", status);
    return(status);
  }
  return KIM_STATUS_OK;
}



/*******************************************************************************
*
* Init KIM API argument values 
*
* potfit global variables:
* atoms
*
* potfit-KIM global variables:
* box_side_len 
*
*******************************************************************************/ 

int init_KIM_API_argument(void* pkim, int Natoms, int Nspecies, int start)
{
  /* local vars */
  /* model inputs */
  int* numberOfParticles;
  int* numberOfSpecies;
  int* particleSpecies;
  double* coords; 
  int* numberContrib;
  /* other locals */
  const char* NBCstr;
  int NBC;
  int status; 
  int species_code; 
  int halfflag;
  int i, j;
  double* boxSideLen;
  int which_conf; /* which config we are in? */

  /* determine which neighbor list type to use */
  halfflag = KIM_API_is_half_neighbors(pkim, &status);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__,"is_half_neighbors", status);
    return(status);   
  }

  /* unpack data from KIM object */
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

  /* set various values */
  *numberOfParticles = Natoms;
  *numberOfSpecies   = Nspecies;  
	if (1==halfflag)
		*numberContrib   = Natoms;

  /* set coords values */
  for (i = 0; i < *numberOfParticles; i++) {
    coords[DIM*i]   = atoms[start+i].pos.x;
    coords[DIM*i+1] = atoms[start+i].pos.y;
    coords[DIM*i+2] = atoms[start+i].pos.z;
  }

  /* set species types */
  for (i = 0; i < *numberOfParticles; i++) { 
    j = atoms[start+i].type;
    species_code = KIM_API_get_species_code(pkim, elements[j], &status); 
    if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_species_code", status);
      return status;
    }
    particleSpecies[i] = species_code;
  }

  /* set boxSideLengths if MI_OPBC is used */
  /* determine which NBC is used */
  status = KIM_API_get_NBC_method(pkim, &NBCstr);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method", status);
    return status;
  }
  if ((!strcmp("NEIGH_RVEC_H",NBCstr)) || (!strcmp("NEIGH_RVEC_F",NBCstr))) {
    NBC = 0;
  }
  else if ((!strcmp("NEIGH_PURE_H",NBCstr)) || (!strcmp("NEIGH_PURE_F",NBCstr))) {
    NBC = 1;
  }
  else if ((!strcmp("MI_OPBC_H",NBCstr)) || (!strcmp("MI_OPBC_F",NBCstr))) {
    NBC = 2;
  }
  else if (!strcmp("CLUSTER",NBCstr)) {
    NBC = 3;
  }
  else {
    status = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "Unknown NBC method", status);
    return status;
  }

  if (NBC == 2) {
    which_conf = atoms[start].conf; 
    /* Unpack data from KIM object */
    KIM_API_getm_data(pkim, &status, 1*3, "boxSideLengths", &boxSideLen, 1);
    if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);
      return status;
    }
    /* set values */
    boxSideLen[0] = box_side_len[DIM*which_conf + 0];
    boxSideLen[1] = box_side_len[DIM*which_conf + 1];
    boxSideLen[2] = box_side_len[DIM*which_conf + 2];
  }   


  return KIM_STATUS_OK;
}


/*******************************************************************************
*
* Register neighborlist object and get_neigh fuction in KIM API object.
*
*******************************************************************************/ 

int setup_neighborlist_KIM_access(void* pkim, NeighObjectType* NeighObject)
{
  /* local variables */
  int status;
  /* register for neighObject */
  KIM_API_setm_data(pkim, &status, 1*4, "neighObject", 1, NeighObject, 1);
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
  return KIM_STATUS_OK;
}


/*******************************************************************************
*
* Create neighborlist and initialize 
*
* potfit global variables:
*
*******************************************************************************/

int init_neighborlist(NeighObjectType* NeighObject, int Natoms, int start)
{
  /* local variables */
  int i, j, k;
  int neighListLength = 0; /* total length of neighList */

  /* calcualte the length of neighborList */
  for (i = 0; i < Natoms; i++) {
    neighListLength += atoms[start+i].num_neigh;  
  }

  /* allocate memory for NeighObject members */
  NeighObject->NNeighbors = (int*) malloc(Natoms*sizeof(int));
  NeighObject->neighborList = (int*) malloc(neighListLength*sizeof(int));
  NeighObject->RijList = (double*) malloc((DIM*neighListLength)*sizeof(double));
  NeighObject->BeginIdx = (int*) malloc(Natoms*sizeof(int));
  if (NULL == NeighObject->NNeighbors || NULL == NeighObject-> neighborList
		||NULL == NeighObject->RijList || NULL == NeighObject->BeginIdx )
	{
    KIM_API_report_error(__LINE__,__FILE__,"malloc unsuccessful", -1);
    exit(1);
  }

	/* free memory */
	reg_for_free(NeighObject->NNeighbors, "NeighObject->NNeighbors");
	reg_for_free(NeighObject->neighborList, "NeighObject->neighborList");
	reg_for_free(NeighObject->RijList, "NeighObject->RijList");
	reg_for_free(NeighObject->BeginIdx, "NeighObject->BeginIdx");

  /* copy the number of neighbors to NeighObject.NNeighbors */
  for (i = 0; i < Natoms; i++) {
    NeighObject->NNeighbors[i] = atoms[start+i].num_neigh;  
  }
  
  /* copy neighborlist from distributed memory locations to continuous ones */
  k = 0;
  for (i = 0; i < Natoms; i++) {
    NeighObject->BeginIdx[i] = k;
    for (j = 0; j < NeighObject->NNeighbors[i]; j++) {
      NeighObject->neighborList[k]  = atoms[start+i].neigh[j].nr - start;
      NeighObject->RijList[DIM*k]   = atoms[start+i].neigh[j].dist.x;
      NeighObject->RijList[DIM*k+1] = atoms[start+i].neigh[j].dist.y;
      NeighObject->RijList[DIM*k+2] = atoms[start+i].neigh[j].dist.z;
      k++;
    }
  }

  /* If the number of neighbors of the last atom is zero, set the BeginIdx to the last
   * position in neghbor list. The purpose is to ensure that the BeginIdx of the
   * last atom will not go beyond the limit of neighborlist length.
   * e.g. say there are 128 atoms in a cluster, and we use half neighbor list,
   * then the 128th atom will have no neighbors. The begin index for the
   * last atom, BeginIdx[127] will go beyond the limit of Neighbor list, which may
   * result in segfault in `get_neigh'. */  
  if( NeighObject->NNeighbors[Natoms-1] == 0) {
    NeighObject->BeginIdx[Natoms-1] = k-1; 
  }

  return KIM_STATUS_OK;
}


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
  if (0 == *mode) { /* iterator mode */ 
    if (0 == *request) { /* reset iterator */
      (*nl).iteratorId = -1;
      return KIM_STATUS_NEIGH_ITER_INIT_OK;
    } else if (1 == *request) { /* increment iterator */
      (*nl).iteratorId++;
      if ((*nl).iteratorId >= *numberOfParticles) {
        return KIM_STATUS_NEIGH_ITER_PAST_END;
      } else {
        partToReturn = (*nl).iteratorId;
      }
    } else { /* invalid request value */
      KIM_API_report_error(__LINE__, __FILE__,"Invalid request in get_periodic_neigh",
                           KIM_STATUS_NEIGH_INVALID_REQUEST);
      return KIM_STATUS_NEIGH_INVALID_REQUEST;
    }
  } else if (1 == *mode) { /* locator mode */
    if ((*request >= *numberOfParticles) || (*request < 0)) { /* invalid id */
      KIM_API_report_error(__LINE__, __FILE__,"Invalid part ID in get_periodic_neigh",
                          KIM_STATUS_PARTICLE_INVALID_ID);
      return KIM_STATUS_PARTICLE_INVALID_ID;
    } else {
      partToReturn = *request;
    }
  } else { /* invalid mode */ 
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
  *nei1part = &(*nl).neighborList[idx];

  /* set the pointer to Rij to appropriate value */
  *Rij = &(*nl).RijList[DIM*idx];

  return KIM_STATUS_OK;
}


/*******************************************************************************
*
* Initialize optimizable parameter 
* allocate memory for free parameters, and nest the optimizable ones in to a 
* variable called `nestedvalue'. Also publish the parameter once for alll. 
*
* potfit global variable:
* rcutmax
*
* potfit-KIM global variable:
* num_opt_param
* name_opt_param
*
*******************************************************************************/

void init_optimizable_param()
{
  /* local vars */
  int i;

  /* allocate memory */
  FreeParamAllConfig = (FreeParamType*) malloc(nconf*sizeof(FreeParamType));
  if (NULL == FreeParamAllConfig) {
    KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
    exit(1);
  }
	reg_for_free(FreeParamAllConfig, "FreeParamAllConfig");	

  for (i = 0; i <nconf; i++) {
    /* nest optimizable parameters */
    get_free_param_double(pkimObj[i], &FreeParamAllConfig[i]);  
    nest_optimizable_param(pkimObj[i], &FreeParamAllConfig[i], 
                               name_opt_param, num_opt_param);
  }
}




/******************************************************************************* 
* 
*  Get NBC, is_half_neighbors and the computation flags of model
*
*******************************************************************************/
void get_compute_const(void* pkim) 
{
  /* local variables */
  int status;
  const char* NBCstr;
  
  /* NBC */
  status = KIM_API_get_NBC_method(pkim, &NBCstr);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method",status);
    exit(1);
  }
  strcpy(NBC_method, NBCstr);
  printf("KIM NBC: %s.\n", NBC_method);

  /* use half or full neighbor list?  1 = half, 0 =full; this will be used
   * in config.c to build up the neighbor list */
  is_half_neighbors = KIM_API_is_half_neighbors(pkim, &status);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_is_half_neighbors",status);
    exit(1);
  }

  /* model ability flag */ 
  get_KIM_model_has_flags();

  return;
}


/******************************************************************************* 
* 
*  get the computation flags of model
*
*******************************************************************************/

int get_KIM_model_has_flags()
{
  void* pkim;
  int status;

  /* get KIM API object representing the KIM Model only */
  status = KIM_API_model_info(&pkim, kim_model_name);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_info", status);
    return(status);
  }

  /* determine if the KIM Model can compute the total energy */
  KIM_API_get_index(pkim, (char*) "energy", &status);
  kim_model_has_energy = (status == KIM_STATUS_OK);

  /* determine if the KIM Model can compute the forces */
  KIM_API_get_index(pkim, (char*) "forces", &status);
  kim_model_has_forces = (status == KIM_STATUS_OK);

  /* determine if the KIM Model can compute the virial */
  KIM_API_get_index(pkim, (char*) "virial", &status);
  kim_model_has_virial = (status == KIM_STATUS_OK);
  KIM_API_get_index(pkim, (char*) "process_dEdr", &status);
  kim_model_has_virial = kim_model_has_virial || (status == KIM_STATUS_OK);

  /* tear down KIM API object */
  KIM_API_free(&pkim, &status);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_free", status);
    return(status);
  }

  return KIM_STATUS_OK;
}


/***************************************************************************
*
* Get free parameters of type double for the descriptor file of the model
*
* put all free parameters of type `double' into the data struct (only type
* `double' parameters are optimizable. `int' ones may be some flag.) However
* even if it is of type double, it is upto the user's desicion whether to 
* optimize it or not.
*
* `PARAM_FREE_cutoff' will also be included in ParamType->name.
*
***************************************************************************/

int get_free_param_double(void* pkim, FreeParamType* FreeParam) 
{
  /*local vars*/
  int status;
  int NumFreeParam;
  int maxStringLength;
  char* pstr;
  char buffer[128];
  char name[64];
  char type[16];
  int NumFreeParamDouble;
  int i;

  /* get the maxStringLength of free parameters */
  status = KIM_API_get_num_free_params(pkim, &NumFreeParam, &maxStringLength);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_num_free_params", status);
    return(status);
  }

  /* get the descriptor file, pointed by pstr. the type of data will be phrased from pstr*/  
  status = KIM_API_get_model_kim_str(kim_model_name, &pstr);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_kim_str", status);
    return(status);
  }

  /* infinite loop to find PARAM_FREE_* of type `double' */
  /* It's safe to do pstr = strstr(pstr+1,"PARAM_FREE") because the ``PARAM_FREE''
    will never ever occur at the beginning of the ``descriptor.kim'' file */  
  FreeParam->name = NULL;
  NumFreeParamDouble = 0;
  while (1) {
    pstr = strstr(pstr+1,"PARAM_FREE");
    if (pstr == NULL) {
      break;
    } else {
      snprintf(buffer, sizeof(buffer), "%s", pstr);
      sscanf(buffer, "%s%s", name, type);
      if (strcmp(type, "double") == 0) {
        NumFreeParamDouble++;          
        FreeParam->name = (char**) realloc(FreeParam->name, (NumFreeParamDouble)*sizeof(char*));
        if (NULL==FreeParam->name) {
          KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
          exit(1);
        }
        /*maxStringLength+1 to hold the `\0' at end*/
        FreeParam->name[NumFreeParamDouble - 1] = (char*) malloc((maxStringLength+1)*sizeof(char));
        if (NULL==FreeParam->name[NumFreeParamDouble - 1]) {
          KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
          exit(1);
        }               
				reg_for_free(FreeParam->name[NumFreeParamDouble-1], "FreeParam->name[NumFreeParamDouble-1]");
        strcpy(FreeParam->name[NumFreeParamDouble - 1], name);
      }
    }
  }
  
	reg_for_free(FreeParam->name, "FreeParam->name");

  FreeParam->Nparam = NumFreeParamDouble;

  /* allocate memory for value */
  FreeParam->value = (double**) malloc(FreeParam->Nparam * sizeof(double*));
  if (NULL==FreeParam->value) {
    KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
    exit(1);
  } 
	reg_for_free(FreeParam->value, "FreeParam->value");

  /* get the pointer to parameter */
  for(i = 0; i < FreeParam->Nparam; i++ ) {
    FreeParam->value[i] = KIM_API_get_data(pkim, FreeParam->name[i], &status);
    if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", status);
      return(status);
    }
  }
  
  /* allocate memory for rank */
  FreeParam->rank = (int*) malloc(FreeParam->Nparam * sizeof(int));
  if (NULL==FreeParam->rank) {
    KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
    exit(1);
  } 
	reg_for_free(FreeParam->rank, "FreeParam->rank");
  
	/* get rank */
  for(i = 0; i < FreeParam->Nparam; i++) {
    FreeParam->rank[i] = KIM_API_get_rank(pkim, FreeParam->name[i], &status);
    if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_rank", status);
      return(status);
    }
  }

  /* allocate memory for shape */
  FreeParam->shape = (int**) malloc(FreeParam->Nparam * sizeof(int*));
  if (NULL==FreeParam->shape) {
    KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
    exit(1);
  }
	reg_for_free(FreeParam->shape, "FreeParam->shape");
  
	for (i = 0; i < FreeParam->Nparam; i++) {
    FreeParam->shape[i] = (int*) malloc(FreeParam->rank[i] * sizeof(int));
    /* should not do the check, since rank may be zero. Then in some implementation
      malloc zero would return NULL */  
    /*  if (NULL==FreeParam->shape[i]) {
    KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
    exit(1);
    }*/ 
	reg_for_free(FreeParam->shape[i], "FreeParam->shape[i]");
  }

  /* get shape */
  for(i = 0; i < FreeParam->Nparam; i++) {
    KIM_API_get_shape(pkim, FreeParam->name[i], FreeParam->shape[i], &status);
    if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_shape", status);
      return(status);
    }
  }
  
  /* nestedvalue is not allocated here, give NULL pointer to it */
  FreeParam->nestedvalue = NULL;

  /* free the memory of model kim str */
  free(pstr);

  return KIM_STATUS_OK;
}


/*****************************************************************************
* 
* Nest optimizable parameters (pointer)
*
* Note that the optimizable parameters is a subset of the free parameters with 
* type `double'. Only the ones that the user list in the potential input file 
* will be optimized. This function should be called after `get_free_param_double'.
*
* input_param_name: the names of parameters that will be optimized 
* input_param_num: the number of parameters that will be optimized
*
*****************************************************************************/

int nest_optimizable_param(void* pkim, FreeParamType* FreeParam,
                       char** input_param_name, int input_param_num)
{
  /* local variables */
  int tmp_size;
  int total_size;           /* the size of the nested values*/
  int have_name;            /* flag, is the name in the data struct? */
  int idx[input_param_num]; /* the index in FreeParam */
  int i, j, k;

  /* nest values  */
  /*first compute the total size of the optimizable parameters */
  total_size = 0;
  for (i = 0; i < input_param_num; i++) { 
    have_name = 0;
    for (j = 0; j< FreeParam->Nparam; j++) {
      if (strcmp(FreeParam->name[j], input_param_name[i]) == 0) {
        idx[i] = j;
        have_name = 1;
        break;
      }
    }
    if (have_name) {
      tmp_size = 1;
      for (j = 0; j < FreeParam->rank[idx[i]]; j++) {
        tmp_size *= FreeParam->shape[idx[i]][j];
      }
      /* store tmp_size for later use in potential outpot */
      size_opt_param[i] = tmp_size; 
      total_size += tmp_size;
    } else {
      error(1, "The parameter '%s' is not optimizable, check spelling.\n",
          input_param_name[i]);
    }
  }
  
  /* allocate memory for nestedvalue*/
  FreeParam->nestedvalue = (double**) malloc((total_size) * sizeof(double*));
  if (NULL==FreeParam->nestedvalue) {
    KIM_API_report_error(__LINE__, __FILE__,"malloc unsuccessful", -1);
    exit(1);
  }
	reg_for_free(FreeParam->nestedvalue, "FreeParam->nestedvalue");
  
	/*copy the values (pointers) to nestedvalue */
  k = 0;
  for (i = 0; i < input_param_num; i++ ) {
    tmp_size = 1; 
    for (j = 0; j < FreeParam->rank[idx[i]]; j++) {
      tmp_size *= FreeParam->shape[idx[i]][j];
    }
   for (j = 0; j <tmp_size; j++) {
      FreeParam->nestedvalue[k] = FreeParam->value[idx[i]]+j;
      k++;
    }
  }

  /* store the number of total parameter values */
  FreeParam->Nnestedvalue = total_size;

  return total_size;
}


/***************************************************************************
* 
* calculate force from KIM (general force, including force, virial, energy)
*
* flags:
* usefore & usestress
*
* general force pointers
* energy;
* force;
* virial;
*
***************************************************************************/

int calc_force_KIM(void* pkim, double** energy, double** force, double** virial,
              int useforce, int usestress)
{ 
  /* local variables */
   int status;

  /* get data */
  KIM_API_getm_data(pkim, &status, 3*3,
                    "energy", energy, 1, 
                    "forces", force,  useforce,
                    "virial", virial, usestress);
  if (KIM_STATUS_OK > status) {
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


/***************************************************************************
*
* publsih cutoff
*
***************************************************************************/
int publish_cutoff(void* pkim, double cutoff) 
{
  /* local variable */
  int status;
  double* pcutoff;
 
  /* update cutoff */
  pcutoff = KIM_API_get_data(pkim, "PARAM_FREE_cutoff", &status);
  if (KIM_STATUS_OK != status) {
#ifndef NOLIMITS /* if defines nolimits; KIM cutoff is always used */
    warning("Publish cutoff to KIM failed. The cutoff in the KIM Model is used.\n"
            "Possibly the KIM Model has no 'PARAM_FREE_cutoff'.\n\n");
#endif 
  } else {
    *pcutoff = cutoff;
    /* reinit KIM model */
    status = KIM_API_model_reinit(pkim);
    if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_reinit", status);
      return status;
    }
  }

  return KIM_STATUS_OK;
}



/***************************************************************************
*
* publish KIM parameter 
*
* transferred varialbes:
* PotTable: potential table where the potential paramenters are stored
*
***************************************************************************/
int publish_param(void* pkim, FreeParamType* FreeParam, double* PotTable)
{ 
  /*local variables*/  
  int status;
  int i;

  /*publish parameters */ 
  for (i = 0; i < FreeParam->Nnestedvalue; i++) {
    *FreeParam->nestedvalue[i] = PotTable[i]; 
  }

  /* reinit KIM model */
  status = KIM_API_model_reinit(pkim);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_reinit", status);
    return status;
  }

  return KIM_STATUS_OK;
}


/******************************************************************************
*
* Read the keywords (`type', `cutoff', etc) and parameter names (beginning 
* with `PARAM_FREE') in the potential input file.
*
* In this function, the memory of three global variables: `num_opt_param',
* `name_opt_param', `size_opt_param' will be allocated and the first two will be
* initialized here (based on the infomation from the input), `size_opt_param' will
* be initialized in nest_optimizable_param. 
*
* pt: potential table
* filename: potential input file name 
* infile: potential input file
* FreeParam: data struct that contains the free parameters info of the KIM Model
* 
******************************************************************************/

int read_potential_keyword(pot_table_t* pt, char* filename, FILE* infile,
													FreeParamType* FreeParam)
{
	/* local variables */
  int   i, j, k, ret_val;
  char  buffer[255], name[255];
  fpos_t startpos; 
	void* pkim;
  int status;

  /* save starting position */
  fgetpos(infile, &startpos);

  /* scan for "type" keyword */
  buffer[0] = '\0';
  name[0] = '\0';
  do {
    fgets(buffer, 255, infile);
    sscanf(buffer, "%s", name);
  } while (strncmp(name, "type", 4) != 0 && !feof(infile));
  if (strncmp(name, "type", 4) != 0) {
    error(1, "Keyword 'type' is missing in file: %s.", filename);
  }
  if (1 > sscanf(buffer, "%*s %s", name))
    error(1, "KIM Model name missing in file: %s.", filename);

  /* copy name*/
  strcpy(kim_model_name, name);
  printf("\nKIM Model: %s.\n", kim_model_name);
 
  /* find `check_kim_opt_param' or `num_opt_param'. The two keywords are mutually
	 * exculsive, which comes first will be read, and the other one will be ignored. */
  fsetpos(infile, &startpos);

  do {
    fgets(buffer, 255, infile);
    sscanf(buffer, "%s", name);
  } while (strcmp(name, "check_kim_opt_param") != 0 
					 && strcmp(name, "num_opt_param") != 0 
           && !feof(infile));
  
	/* read `check_kim_opt_param' or `num_opt_param' */
	if (strncmp(buffer,"check_kim_opt_param", 19) == 0) {

  	/* create a temporary KIM objects to query the info */
		/* write temporary descriptor file */
		write_temporary_descriptor_file(kim_model_name);
		/* create KIM object with 1 atom and 1 species */
		status = setup_KIM_API_object(&pkim, 1, 1, kim_model_name);
		if (KIM_STATUS_OK > status) {
      KIM_API_report_error(__LINE__, __FILE__, "setup_KIM_API_object", status);
      exit(1);
    }
		/* initialze the data struct for the free parameters with type double */
		get_free_param_double(pkim, FreeParam);

		printf(" - The following potential parameters are available to fit. Include the "
        "name(s) (and the initial value(s) and lower and upper boundaries if "
        "`NOLIMITS' is not enabled while compilation) that you want to optimize "
        "in file: %s.\n",filename);
    printf("         param name                 param extent\n");
    printf("        ############               ##############\n");
    for(k = 0; k < FreeParam->Nparam; k++ ) {
      if (strncmp(FreeParam->name[k], "PARAM_FREE_cutoff", 17) == 0){
        continue;
      }
      printf("     %-35s[ ", FreeParam->name[k]);
      for(j = 0; j < FreeParam->rank[k]; j++) {
        printf("%d ", FreeParam->shape[k][j]);
      }
      printf("]\n");
   }
    printf("\n - Note that empty parameter extent (i.e. '[ ]') indicates that the "
           "parameter is a scalar.\n");
    printf(" - Note also that KIM array parameter is row based, while listing the "
        "initial values for such parameter, you should ensure that the sequence is "
        "correct. For example, if the extent of a parameter `PARAM_FREE_A' is "
        "[ 2 2 ], then you should list the initial values as: A[0 0], A[0 1], "
        "A[1 0], A[1 1].\n");

	 /* free the temporary kim model */
	  free_model_object(&pkim);
	  exit(1); 

  } else if (strncmp(buffer,"num_opt_param", 13) == 0) {
    if(1 != sscanf(buffer, "%*s%d", &num_opt_param)) {
      error(1, "Cannot read 'num_opt_param' in file: %s.", filename);  
    }
  } else {
    error(1, "Keyword 'num_opt_param' is missing in file: %s.", filename);
  }

  /* allocate memory */ 
  name_opt_param = (char**)malloc(num_opt_param*sizeof(char*)); 
  size_opt_param = (int*)malloc(num_opt_param*sizeof(int)); 
  reg_for_free(name_opt_param, "name_opt_param");
  reg_for_free(size_opt_param, "size_opt_param");
  if (NULL == name_opt_param || NULL == size_opt_param) {
    error(1, "Error in allocating memory for parameter name or size");
  }
  for (i = 0; i < num_opt_param; i++) {
    name_opt_param[i] = (char *)malloc(255 * sizeof(char));
    reg_for_free(name_opt_param[i], "name_opt_param[%d]", i);
    if (NULL == name_opt_param[i]) {
      error(1, "Error in allocating memory for parameter name");
    }
  }

	/* find parameter names beginning with `PARAM_FREE_*' */
	fsetpos(infile, &startpos);
	for (j = 0; j < num_opt_param; j++) {
		buffer[0] = '\0';
		name[0] = '\0';
		do {
			fgets(buffer, 255, infile);
			ret_val = sscanf(buffer, "%s", name);
		} while (strncmp(name, "PARAM_FREE", 10) != 0 && !feof(infile));
		if (feof(infile) ) {
			error(0, "Not enough parameter(s) 'PARAM_FREE_*' in file: %s.\n", filename);
			error(1, "You listed %d parameter(s), but required are %d.\n", j, num_opt_param);
		}
		if (ret_val == 1) {
			strcpy(name_opt_param[j], name); 
		} else {
			error(0, "parameter '%d' in file '%s' corrupted\n.", j + 1, filename);
			error(1, "Each parameter name should be in a single line.\n");
		}
	}

	return 0;
}


/*******************************************************************************
 *
 *  write potential table (format 5)
 *
 ******************************************************************************/

void write_pot_table5(pot_table_t* pt, char *filename)
{
  /*local variables */
  FILE *outfile = NULL;
  int i, j, k;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 5 1");
  if (have_elements) {
    fprintf(outfile, "\n#C");
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    fprintf(outfile, "\n##");
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
        fprintf(outfile, " %s-%s", elements[i], elements[j]);
  }
  fprintf(outfile, "\n#E");

  /* write KIM Model name */
  fprintf(outfile, "\n\n# KIM Model name");
  fprintf(outfile, "\ntype  %s", kim_model_name);

  /* write cutoff */
  fprintf(outfile, "\n\n# cutoff");
  fprintf(outfile, "\ncutoff  %24.16e", rcutmax);

  /* check KIM optimizable params */
  fprintf(outfile, "\n\n# uncomment the following line to check the optimizable "
                   "parameters of the KIM Model");
  fprintf(outfile, "\n#check_kim_opt_param");

  /* number of opt params */
  fprintf(outfile, "\n\n# the number of optimizable parameters that will be listed below");
  fprintf(outfile, "\nnum_opt_param  %d", num_opt_param);
  
  /* write data */
  k = 0; 
  fprintf(outfile, "\n\n# parameters");
  for (i = 0; i < num_opt_param; i++) {
    fprintf(outfile, "\n%s", name_opt_param[i]);
#ifndef NOLIMITS
    for (j = 0; j < size_opt_param[i]; j++) {
      fprintf(outfile, "\n%18.10e %18.10e %18.10e", pt->table[k], apot_table.pmin[0][k],
              apot_table.pmax[0][k]); 
      k++;
    }
    fprintf(outfile, "\n");
#endif  
  }

  fclose(outfile);
}


/******************************************************************************
*
* write final descriptor file
*
* arguments:
*
* u_f: flag; use force?
* u_s: flag; use stress? 
*
******************************************************************************/

int write_final_descriptor_file(int u_f, int u_s)
{
  int compute_energy;
  int compute_forces;
  int compute_virial;
  int status;


  /* set_compute flag */
  if (!kim_model_has_energy) {
    compute_energy = 0;
    error(1,"KIM Model does not provide `energy'.\n");
  } else {
    compute_energy = 1;
  }

  if (u_f) {
    if (kim_model_has_forces) {
      compute_forces = 1;
    } else {
      compute_forces = 0;
      error(1,"KIM Model does not provide `forces'.\n");
    } 
  }else {
    compute_forces = 0;
  }

  if (u_s) {
    if (kim_model_has_virial) {
      compute_virial = 1;
    } else {
      compute_virial = 0;
      error(1,"KIM Model does not provide `virial'.\n");
    } 
  }else {
    compute_virial = 0;
  }

  /* write the `descritor.kim' file for this test */
  status = write_descriptor_file(ntypes, elements, 
                              compute_energy, compute_forces, compute_virial); 
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "write_descriptor_file", status);
    return(status);
  }

  return KIM_STATUS_OK;
}


/******************************************************************************
 * 
 * write temporary `descriptor.kim' file for the test 
 *
 * arguments:
 * Nspecies: number of species that will be written into the descriptor file
 * species: the names of the species, e.g. Al, Cu 
 *
 *****************************************************************************/

int write_temporary_descriptor_file(char* modelname)
{
	/* local variables */
	void* pkim;
	char* species[1];
	int status;
	int Nspecies = 1;

  /* create a temporary object*/
  status = KIM_API_model_info(&pkim, modelname);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_info", status);
    return(status);
  }

  /* get the first species supported by the model */
  status = KIM_API_get_model_species(pkim, 0, &species[0]);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_species", status);
    return(status);
  }
 
  /* write a temporary `descriptor.kim' file, used only to query model info */
  /* we'll write `descriptor.kim' file with the species reading from potfit later. */
  status = write_descriptor_file(Nspecies, species, 0, 0, 0); 
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "write_descriptor_file", status);
    return(status);
  }

  /* free the temporary object */
  KIM_API_free(&pkim, &status);
  if (KIM_STATUS_OK > status) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_free", status);
    return(status);
  }

	return KIM_STATUS_OK;
}



/******************************************************************************
 * 
 * write descriptor.kim for the test 
 *
 * arguments:
 * Nspecies: number of species that will be written into the descriptor file
 * species: the names of the species, e.g. Al, Cu 
 * compute_erergy: flag; whether to write 'energy' in the Model Output section
 * compute_forces: flag; whether to write 'forces' in the Model Output section
 * compute_virial: flag; whether to write 'virial' in the Model Output section
 *
 *****************************************************************************/

int write_descriptor_file(int Nspecies, char** species, int compute_energy, 
                          int compute_forces, int compute_virial)
{
  /* local variables */
  int i;
  FILE* outfile;

  outfile = fopen("descriptor.kim", "w");
 
  /* header */
  fprintf(outfile,
    "#\n"
    "# CDDL HEADER START\n"
    "#\n"
    "# The contents of this file are subject to the terms of the Common Development\n"
    "# and Distribution License Version 1.0 (the \"License\").\n"
    "#\n"
    "# You can obtain a copy of the license at\n"
    "# http://www.opensource.org/licenses/CDDL-1.0.  See the License for the\n"
    "# specific language governing permissions and limitations under the License.\n"
    "#\n"
    "# When distributing Covered Code, include this CDDL HEADER in each file and\n"
    "# include the License file in a prominent location with the name LICENSE.CDDL.\n"
    "# If applicable, add the following below this CDDL HEADER, with the fields\n"
    "# enclosed by brackets \"[]\" replaced with your own identifying information:\n"
    "#\n"
    "# Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.\n"
    "#\n"
    "# CDDL HEADER END\n"
    "#\n\n"
  );

  /* copyright */
  fprintf(outfile, 
    "#\n"
    "# Copyright (c) 2013--2014, Regents of the University of Minnesota.\n"
    "# All rights reserved.\n"
    "#\n"
    "# Contributors:\n"
    "#    Ryan S. Elliott\n"
    "#    Ellad B. Tadmor\n"
    "#    Valeriu Smirichinski\n"
    "#    Stephen M. Whalen\n"
    "#\n\n"
  );

  /* versioni and units  */
  fprintf(outfile,
    "#######################################################################################################\n"
    "#\n"
    "# Release: This file is part of the kim-api-v1.6.3 package.\n"
    "#\n"
    "# See src/standard.kim for documentation about this file\n"
    "#\n"
    "#######################################################################################################\n\n"
    "KIM_API_Version := 1.6.3\n\n"
    "Unit_length      := A\n"
    "Unit_energy      := eV\n"
    "Unit_charge      := e\n"
    "Unit_temperature := K\n"
    "Unit_time        := ps\n\n\n"
  );

  /* particle species */
  /* code does not matter, so just give it 0 */
  fprintf(outfile, 
    "#######################################################################################################\n"
    "PARTICLE_SPECIES:\n"
    "# Symbol/name               Type                    code\n\n" 
  );
  for (i = 0; i < Nspecies; i++) {
    fprintf(outfile, "%s                          spec                    0\n\n", species[i]);
  }

  /* conversions */
  fprintf(outfile, 
    "\n#######################################################################################################\n"
    "CONVENTIONS:\n"
    "# Name                      Type\n\n"
    "ZeroBasedLists              flag\n\n"
    "Neigh_LocaAccess            flag\n\n"
    "MI_OPBC_H                   flag\n\n"
    "MI_OPBC_F                   flag\n\n"
    "NEIGH_RVEC_H                flag\n\n"
    "NEIGH_RVEC_F                flag\n\n\n"
    );

  /* Model output */
  fprintf(outfile, 
    "#######################################################################################################\n"
    "MODEL_INPUT:\n"
    "# Name                      Type         Unit                Shape             Requirements\n\n"
    "numberOfParticles           integer      none                []\n\n"
    "numberOfSpecies             integer      none                []\n\n"
    "particleSpecies             integer      none                [numberOfParticles]\n\n"
    "coordinates                 double       length              [numberOfParticles,3]\n\n"
    "boxSideLengths              double       length              [3]\n\n"
    "numberContributingParticles integer      none                []\n\n"
    "get_neigh                   method       none                []\n\n"
    "neighObject                 pointer      none                []\n\n\n"
  );

  /* Model output */
  fprintf(outfile,  
    "#######################################################################################################\n"
    "MODEL_OUTPUT:\n"
    "# Name                      Type         Unit                Shape\n\n"
    "# Requirements\n\n" 
    "destroy                     method       none                []\n\n"
    "compute                     method       none                []\n\n"
    "reinit                      method       none                []\n\n"
    "cutoff                      double       length              []\n\n"
  );
 
  if(compute_energy)
    fprintf(outfile,
      "energy                      double       energy              []\n\n"
    ); 
  if(compute_forces)
    fprintf(outfile,
      "forces                      double       force               [numberOfParticles,3]\n\n"
    ); 
  if(compute_virial)
    fprintf(outfile,
      "virial                      double       energy              [6]" 
  );
  
  fflush(outfile);
  fclose(outfile);

  return KIM_STATUS_OK;
}


/******************************************************************************
 * free KIM model and object 
 ******************************************************************************/
int free_model_object(void** pkim) 
{
	/* local variables */
	int status;

	/* call model destroy */
	status = KIM_API_model_destroy(*pkim);
	if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__,"KIM_API_model_destroy", status);
		return status;
	}
	/* free KIM objects */
	KIM_API_free(pkim, &status);
	if (KIM_STATUS_OK > status) {
		KIM_API_report_error(__LINE__, __FILE__,"KIM_API_free", status);
		return status;
	}
	return KIM_STATUS_OK;
}


/******************************************************************************
 * free KIM model and object 
 ******************************************************************************/
void free_KIM()
{
	/* local variables */
	int i;
	
	for(i = 0; i < nconf; i++) {
		free_model_object(&pkimObj[i]);
	}	
}



