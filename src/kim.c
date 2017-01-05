/****************************************************************
 *
 * kim.h: KIM interface
 *
 ****************************************************************
 *
 * Copyright 2002-2017 - the potfit development team
 *
 * https://www.potfit.net/
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

#include <stdio.h>

#include "potfit.h"

#include "kim.h"
#include "memory.h"
#include "utils.h"

typedef struct {
  int iteratorId;
  int* NNeighbors;
  int* neighborList;
  double* RijList;
  int* BeginIdx;  /* The position of first neighbor of each atom in the neighbor list */
} neigh_obj_t;

#define KIM_INPUT_TEMPLATE_HEADER \
    "KIM_API_Version := 1.7.3\n" \
    "Unit_length      := A\n" \
    "Unit_energy      := eV\n" \
    "Unit_charge      := e\n" \
    "Unit_temperature := K\n" \
    "Unit_time        := ps\n" \
    "PARTICLE_SPECIES:\n"

#define KIM_INPUT_TEMPLATE_CONVENTIONS \
    "CONVENTIONS:\n" \
    "ZeroBasedLists              flag\n" \
    "Neigh_LocaAccess            flag\n" \
    "NEIGH_RVEC_H                flag\n" \
    "NEIGH_RVEC_F                flag\n" \
    "MI_OPBC_H                   flag\n" \
    "MI_OPBC_F                   flag\n"

#define KIM_INPUT_TEMPLATE_INPUT \
    "MODEL_INPUT:\n" \
    "numberOfParticles           integer    none     []\n" \
    "numberOfSpecies             integer    none     []\n" \
    "particleSpecies             integer    none     [numberOfParticles]\n" \
    "coordinates                 double     length   [numberOfParticles,3]\n" \
    "boxSideLengths              double     length   [3]\n" \
    "numberContributingParticles integer    none     []\n" \
    "get_neigh                   method     none     []\n" \
    "neighObject                 pointer    none     []\n" \

#define KIM_INPUT_TEMPLATE_OUTPUT \
    "MODEL_OUTPUT:\n" \
    "destroy                     method     none     []\n" \
    "compute                     method     none     []\n" \
    "reinit                      method     none     [] optional\n" \
    "cutoff                      double     length   []\n"

void create_descriptor_for_config(int config_index, char** descriptor)
{
  *descriptor = (char*)malloc(8 * 1024 * sizeof(char));
  if (!descriptor)
    error(1, "Error allocating descriptor string!\n");
  memset(*descriptor, 0, 8 * 1024 * sizeof(char));

  // write header

  char* p = *descriptor;

  strncpy(p, KIM_INPUT_TEMPLATE_HEADER, strlen(KIM_INPUT_TEMPLATE_HEADER));
  p += strlen(KIM_INPUT_TEMPLATE_HEADER);

  // write species data

  for (int i = 0; i < g_param.ntypes; ++i) {
    char temp[255];
    snprintf(temp, 255, "%s                      spec                0\n", g_config.elements[i]);
    strncpy(p, temp, strlen(temp));
    p += strlen(temp);
  }

  // write conventions

  strncpy(p, KIM_INPUT_TEMPLATE_CONVENTIONS, strlen(KIM_INPUT_TEMPLATE_CONVENTIONS));
  p += strlen(KIM_INPUT_TEMPLATE_CONVENTIONS);

  // write input

  strncpy(p, KIM_INPUT_TEMPLATE_INPUT, strlen(KIM_INPUT_TEMPLATE_INPUT));
  p += strlen(KIM_INPUT_TEMPLATE_INPUT);

  // write output

  strncpy(p, KIM_INPUT_TEMPLATE_OUTPUT, strlen(KIM_INPUT_TEMPLATE_OUTPUT));
  p += strlen(KIM_INPUT_TEMPLATE_OUTPUT);

  // write energe, forces and virial

  if (!g_kim.model_has_energy) {
    error(1,"KIM Model does not provide 'energy'.\n");
  } else {
    char temp[255];
    snprintf(temp, 255, "energy                     double     energy      []\n");
    strncpy(p, temp, strlen(temp));
    p += strlen(temp);
  }

  if (g_config.useforce[config_index]) {
    if (!g_kim.model_has_forces)
      error(1, "KIM model does not provide 'forces'.\n");
    char temp[255];
    snprintf(temp, 255, "forces                     double     force       [numberOfParticles,3]\n");
    strncpy(p, temp, strlen(temp));
    p += strlen(temp);
  }

#if defined(STRESS)
  if (g_config.usestress[config_index]) {
    if (!g_kim.model_has_virial)
      error(1, "KIM model does not provide 'virial'.\n");
    char temp[255];
    snprintf(temp, 255, "virial                     double     energy      [6]");
    strncpy(p, temp, strlen(temp));
    p += strlen(temp);
  }
#endif
}

void set_kim_model_data(int config_index)
{
  int status;
  int* numberOfParticles;
  int* numberOfSpecies;
  int* particleSpecies;
  double* coords;
  int* numberContrib;

  // unpack data from KIM object
  KIM_API_getm_data(g_kim.pkim[config_index],       &status,            5 * 3,
                    "numberOfParticles",            &numberOfParticles, 1,
                    "numberOfSpecies",              &numberOfSpecies,   1,
                    "particleSpecies",              &particleSpecies,   1,
                    "coordinates",                  &coords,            1,
                    "numberContributingParticles",  &numberContrib,     g_kim.is_half_neighbors);
  if (KIM_STATUS_OK > status)
    error(1, "KIM_API_getm_data failed!\n");

  // set various values
  *numberOfParticles = g_config.inconf[config_index];
  *numberOfSpecies   = g_param.ntypes;
  if (g_kim.is_half_neighbors)
    *numberContrib   = g_config.inconf[config_index];

  // set coords values
  for (int i = 0; i < *numberOfParticles; ++i) {
    coords[DIM * i + 0] = g_config.atoms[g_config.cnfstart[config_index] + i].pos.x;
    coords[DIM * i + 1] = g_config.atoms[g_config.cnfstart[config_index] + i].pos.y;
    coords[DIM * i + 2] = g_config.atoms[g_config.cnfstart[config_index] + i].pos.z;
  }

  // set species types
  for (int i = 0; i < *numberOfParticles; ++i) {
    const int j = g_config.atoms[g_config.cnfstart[config_index] + i].type;
    int species_code = KIM_API_get_species_code(g_kim.pkim[config_index], g_config.elements[j], &status);
    if (KIM_STATUS_OK > status)
      error(1, "");
    particleSpecies[i] = species_code;
  }

  // handle certain PBC
  if (g_kim.NBC == KIM_NEIGHBOR_TYPE_OPBC) {
    int which_conf = DIM * g_config.atoms[g_config.cnfstart[config_index]].conf;
    double* boxSideLen = NULL;
    // Unpack data from KIM object
    KIM_API_getm_data(g_kim.pkim[config_index],   &status,      1 * 3,
                      "boxSideLengths",           &boxSideLen,  1);
    if (KIM_STATUS_OK > status)
      error(1, "\n");

    // set values
    boxSideLen[0] = g_kim.box_vectors[which_conf + 0];
    boxSideLen[1] = g_kim.box_vectors[which_conf + 1];
    boxSideLen[2] = g_kim.box_vectors[which_conf + 2];
  }
}

/*******************************************************************************
*
* Create neighborlist and initialize
*
* potfit global variables:
*
*******************************************************************************/

void init_neighborlist(neigh_obj_t* pbuf, int config_index)
{
  const int start = g_config.cnfstart[config_index];
  int neighListLength = 0;

  // calcualte the length of neighborList
  for (int i = 0; i < g_config.inconf[config_index]; ++i)
    neighListLength += g_config.atoms[start + i].num_neigh;

  // allocate memory for pbuf members
  pbuf->NNeighbors = (int*)Malloc(g_config.inconf[config_index] * sizeof(int));
  pbuf->neighborList = (int*)Malloc(neighListLength * sizeof(int));
  pbuf->RijList = (double*)Malloc(DIM * neighListLength * sizeof(double));
  pbuf->BeginIdx = (int*)Malloc(g_config.inconf[config_index] * sizeof(int));

  // copy the number of neighbors to pbuf.NNeighbors
  for (int i = 0; i < g_config.inconf[config_index]; ++i)
    pbuf->NNeighbors[i] = g_config.atoms[start+i].num_neigh;

  // copy neighborlist from distributed memory locations to continuous ones
  int k = 0;
  for (int i = 0; i < g_config.inconf[config_index]; ++i) {
    pbuf->BeginIdx[i] = k;
    for (int j = 0; j < pbuf->NNeighbors[i]; ++j) {
      pbuf->neighborList[k]      = g_config.atoms[start + i].neigh[j].nr - start;
      pbuf->RijList[DIM * k + 0] = g_config.atoms[start + i].neigh[j].dist.x;
      pbuf->RijList[DIM * k + 1] = g_config.atoms[start + i].neigh[j].dist.y;
      pbuf->RijList[DIM * k + 2] = g_config.atoms[start + i].neigh[j].dist.z;
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

  if (pbuf->NNeighbors[g_config.inconf[config_index] - 1] == 0)
    pbuf->BeginIdx[g_config.inconf[config_index] - 1] = k - 1;
}

/***************************************************************************
*
* get_neigh
*
***************************************************************************/

int get_neigh(void* kimmdl, int* mode, int* request, int* part,
                       int* numnei, int** nei1part, double** Rij)
{
   /* local variables */
  intptr_t* pkim = *((intptr_t**) kimmdl);
  int partToReturn;
  int status;
  int* numberOfParticles;
  int idx; /* index of first neighbor of each particle*/

  // reset numnei
  *numnei = 0;

  neigh_obj_t* pbuf = (neigh_obj_t*)KIM_API_get_data(pkim, "neighObject", &status);
  if (KIM_STATUS_OK > status)
    error(1, "KIM_API_get_sim_buffer failed!\n");

  // unpack neighbor list object
  numberOfParticles = (int*) KIM_API_get_data(pkim, "numberOfParticles", &status);
  if (KIM_STATUS_OK > status)
    error(1, "KIM_API_get_data failed!\n");

  /* check mode and request */
  if (0 == *mode) { /* iterator mode */
    if (0 == *request) { /* reset iterator */
      pbuf->iteratorId = -1;
      return KIM_STATUS_NEIGH_ITER_INIT_OK;
    } else if (1 == *request) { /* increment iterator */
      pbuf->iteratorId++;
      if (pbuf->iteratorId >= *numberOfParticles) {
        return KIM_STATUS_NEIGH_ITER_PAST_END;
      } else {
        partToReturn = pbuf->iteratorId;
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
  idx = pbuf->BeginIdx[partToReturn];

  /* set the returned part */
  *part = partToReturn;

  /* set the returned number of neighbors for the returned part */
  *numnei = pbuf->NNeighbors[partToReturn];

  /* set the location for the returned neighbor list */
  *nei1part = &pbuf->neighborList[idx];

  /* set the pointer to Rij to appropriate value */
  *Rij = &pbuf->RijList[DIM*idx];

  return KIM_STATUS_OK;
}

/*******************************************************************************
*
* Create KIM objects, and init the argument values
*
*******************************************************************************/

void create_KIM_objects()
{
  // Allocate memory for KIM objects
  g_kim.pkim = (void**)Malloc(g_config.nconf * sizeof(void *));
  g_kim.param_value = (void***)Malloc(g_config.nconf * sizeof(void**));

  for (int i = 0; i < g_config.nconf; i++) {
    char* descriptor = NULL;
    create_descriptor_for_config(i, &descriptor);

    // initialize KIM API object
    int status = KIM_API_string_init(&g_kim.pkim[i], descriptor, g_kim.model_name);
    if (KIM_STATUS_OK > status)
      error(1, "KIM_API_string_init failed!\n");

    free(descriptor);

    // Allocate memory for each data argument of initialized KIM object
    KIM_API_allocate(g_kim.pkim[i], g_config.inconf[i], g_param.ntypes, &status);
    if (KIM_STATUS_OK > status)
      error(1, "KIM_API_allocate failed!\n");

    // call Model's init routine
    status = KIM_API_model_init(g_kim.pkim[i]);
    if (KIM_STATUS_OK > status)
      error(1, "KIM_API_model_init failed!\n");

    set_kim_model_data(i);

    neigh_obj_t* neigh = (neigh_obj_t*)Malloc(sizeof(neigh_obj_t));

    // store the neighbor object
    status = KIM_API_set_data(g_kim.pkim[i], "neighObject", 1, neigh);
    if (KIM_STATUS_OK > status)
      error(1, "KIM_API_set_data failed!\n");

    // initialize neighbor list
    init_neighborlist(neigh, i);

    // register for get_neigh
    status = KIM_API_set_method(g_kim.pkim[i], "get_neigh", 1, (func_ptr) &get_neigh);
    if (KIM_STATUS_OK > status)
      error(1, "KIM_API_set_method failed!\n");

    // store value pointers
    g_kim.param_value[i] = (void**)Malloc(g_kim.num_opt_param * sizeof(void*));
    for (int j = 0; j < g_kim.num_opt_param; ++j) {
      const int idx = g_kim.idx_opt_param[j];
      int status = KIM_STATUS_FAIL;
      void* data = KIM_API_get_data(g_kim.pkim[i], g_kim.freeparams.name[idx], &status);
      if (status != KIM_STATUS_OK)
        error(1, "KIM_API_get_data failed!\n");
      g_kim.param_value[i][j] = data;
    }
  }
}

/*******************************************************************************
*
* Init KIM objects, each object for a reference configuration.
* Init optimizable parameters; nest them a single variable.
*
*******************************************************************************/

void initialize_KIM()
{
  printf("\nInitializing KIM interface ...\n");

  // create KIM objects and do the necessary initialization
  create_KIM_objects();

  for (int i = 0; i <g_config.nconf; i++) {
    if (g_kim.freeparams.cutoff_is_free_param) {
      int status;
      double* pcutoff = KIM_API_get_data(g_kim.pkim[i], "PARAM_FREE_cutoff", &status);
      if (KIM_STATUS_OK != status)
        error(1, "KIM_API_get_data failed\n");

      *pcutoff = g_config.rcutmax;

      // reinit KIM model
      status = KIM_API_model_reinit(g_kim.pkim[i]);
      if (KIM_STATUS_OK > status)
        error(1, "KIM_API_model_reinit failed!\n");
    }
  }

  printf("\nInitializing KIM interface ... done\n");
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

// int nest_optimizable_param(int config_index)
// {
//   /* local variables */
//   int tmp_size;
//   int total_size;           /* the size of the nested values*/
//   int have_name;            /* flag, is the name in the data struct? */
//   int idx[input_param_num]; /* the index in FreeParam */
//   int i, j, k;
//
//   /* nest values  */
//   /*first compute the total size of the optimizable parameters */
//   total_size = 0;
//   for (i = 0; i < input_param_num; i++) {
//     have_name = 0;
//     for (j = 0; j< FreeParam->Nparam; j++) {
//       if (strcmp(FreeParam->name[j], input_param_name[i]) == 0) {
//         idx[i] = j;
//         have_name = 1;
//         break;
//       }
//     }
//     if (have_name) {
//       tmp_size = 1;
//       for (j = 0; j < FreeParam->rank[idx[i]]; j++) {
//         tmp_size *= FreeParam->shape[idx[i]][j];
//       }
//       /* store tmp_size for later use in potential outpot */
//       g_kim.size_opt_param[i] = tmp_size;
//       total_size += tmp_size;
//     } else {
//       error(1, "The parameter '%s' is not optimizable, check spelling.\n",
//           input_param_name[i]);
//     }
//   }
//
//   /* allocate memory for nestedvalue*/
//   FreeParam->nestedvalue = (double**) Malloc((total_size) * sizeof(double*));
//
// 	/*copy the values (pointers) to nestedvalue */
//   k = 0;
//   for (i = 0; i < input_param_num; i++ ) {
// /*    tmp_size = 1;
//     for (j = 0; j < FreeParam->rank[idx[i]]; j++) {
//       tmp_size *= FreeParam->shape[idx[i]][j];
//     }
//     for (j = 0; j <tmp_size; j++) {
// */
//     for (j = 0; j < g_kim.size_opt_param[i]; j++) {
//       FreeParam->nestedvalue[k] = FreeParam->value[idx[i]]+j;
//       k++;
//     }
//   }
//
//   /* store the number of total parameter values */
//   FreeParam->Nnestedvalue = total_size;
//
//   return total_size;
// }





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

void shutdown_KIM()
{
	/* local variables */
	int i;

	for(i = 0; i < g_config.nconf; i++) {
		free_model_object(&g_kim.pkim[i]);
	}
}
