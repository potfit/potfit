/****************************************************************
 *
 * kim.h: KIM interface
 *
 ****************************************************************
 *
 * Copyright 2002-2018 - the potfit development team
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

#include "potfit.h"

#include "force.h"
#include "kim.h"
#include "memory.h"

void check_model_routines()
{
  int present = 0;
  int required = 0;

  int res = KIM_Model_IsRoutinePresent(g_kim.model, KIM_MODEL_ROUTINE_NAME_Extension, &present, &required);
  if (res)
    error(1, "KIM_Model_IsRoutinePresent failed: %d\n", res);

  if (present)
    g_kim.supported_routines |= KIM_MODEL_ROUTINE_EXTENSION_SUPPORTED;

  res = KIM_Model_IsRoutinePresent(g_kim.model, KIM_MODEL_ROUTINE_NAME_Refresh, &present, &required);
  if (res)
    error(1, "KIM_Model_IsRoutinePresent failed: %d\n", res);

  if (present)
    g_kim.supported_routines |= KIM_MODEL_ROUTINE_REFRESH_SUPPORTED;
  else
    error(1, "KIM Model %s does not support Refresh routine", g_kim.model_name);

  res = KIM_Model_IsRoutinePresent(g_kim.model, KIM_MODEL_ROUTINE_NAME_WriteParameterizedModel, &present, &required);
  if (res)
    error(1, "KIM_Model_IsRoutinePresent failed: %d\n", res);

  if (present)
    g_kim.supported_routines |= KIM_MODEL_ROUTINE_WRITEPARAMS_SUPPORTED;
}

/****************************************************************
  init_kim_model
****************************************************************/

void init_kim_model()
{
  int requestedUnitsAccepted = 0;

  int res = KIM_Model_Create(KIM_NUMBERING_zeroBased, KIM_LENGTH_UNIT_A, KIM_ENERGY_UNIT_eV,
                             KIM_CHARGE_UNIT_e, KIM_TEMPERATURE_UNIT_unused, KIM_TIME_UNIT_unused,
                             g_kim.model_name, &requestedUnitsAccepted, &g_kim.model);

  if (res)
    error(1, "KIM_Model_Create failed: %d\n", res);

  if (!requestedUnitsAccepted)
    error(1, "Model %s does not support the required units!\n", g_kim.model_name);

  int num_supported_species = 0;
  g_kim.nspecies = 0;

  check_model_routines();

  KIM_SPECIES_NAME_GetNumberOfSpeciesNames(&num_supported_species);

  for (int i = 0; i < num_supported_species; ++i) {
    KIM_SpeciesName species_name;
    res = KIM_SPECIES_NAME_GetSpeciesName(i, &species_name);
    if (res)
      error(1, "Cannot get species name %d\n", i);
    int supported = 0;
    int code = 0;
    res = KIM_Model_GetSpeciesSupportAndCode(g_kim.model, species_name, &supported, &code);
    if (res)
      error(1, "Cannot get species support for species %d\n", i);
    if (!supported)
      continue;
    g_kim.nspecies++;
    g_kim.species = (KIM_SpeciesName*)Realloc(g_kim.species, g_kim.nspecies * sizeof(KIM_SpeciesName));
    g_kim.species[g_kim.nspecies - 1] = species_name;
  }

  KIM_Model_GetNumberOfParameters(g_kim.model, &g_kim.nparams);
  if (g_kim.nparams < 1)
    error(1, "Selected KIM model \"%s\" does not support any parameters for optimization!\n", g_kim.model_name);

  printf("\nKIM model \"%s\" supports %d species and %d parameters.\n\n", g_kim.model_name, g_kim.nspecies, g_kim.nparams);

  g_kim.params = (kim_parameter_t*)Malloc(g_kim.nparams * sizeof(kim_parameter_t));

  for (int i = 0; i < g_kim.nparams; ++i) {
    res = KIM_Model_GetParameterMetadata(g_kim.model, i, &g_kim.params[i].type, &g_kim.params[i].extent, &g_kim.params[i].name, &g_kim.params[i].desc);
    if (res)
      error(1, "Error getting parameter type and description for parameter %d\n", i);
    g_kim.total_params += g_kim.params[i].extent;
    g_kim.params[i].values = (value_t*)Malloc(g_kim.params[i].extent * sizeof(value_t));
    for (int j = 0; j < g_kim.params[i].extent; ++j) {
      if (KIM_DataType_Equal(g_kim.params[i].type, KIM_DATA_TYPE_Integer))
        res = KIM_Model_GetParameterInteger(g_kim.model, i, j, &g_kim.params[i].values[j].i);
      else if (KIM_DataType_Equal(g_kim.params[i].type, KIM_DATA_TYPE_Double))
        res = KIM_Model_GetParameterDouble(g_kim.model, i, j, &g_kim.params[i].values[j].d);
      else
        error(1, "Unknown parameter type for parameter %d: %s", i, KIM_DataType_ToString(g_kim.params[i].type));
    }
  }
}

/*******************************************************************************
    check_KIM_model_support
*******************************************************************************/

void check_KIM_model_support()
{
  int res = 0;

  g_kim.species_map = (int*)Malloc(g_param.ntypes * sizeof(int));

  for (int i = 0; i < g_param.ntypes; ++i) {
    KIM_SpeciesName name = KIM_SpeciesName_FromString(g_config.elements[i]);
    if (name.speciesNameID == -1)
      error(1, "Unknown species found: %s\n", g_config.elements[i]);

    // check species
    int speciesIsSupported = 0;
    int modelCode = -1;
    res = KIM_Model_GetSpeciesSupportAndCode(g_kim.model, name, &speciesIsSupported, &modelCode);
    if (res || !speciesIsSupported)
      error(1, "Species %s not supported", g_config.elements[i]);
    g_kim.species_map[i] = modelCode;
  }

  // check if forces and stresses are actually required
  int use_forces = 0;
#if defined(STRESS)
  int use_stress = 0;
#endif
  for (int j = 0; j < g_config.nconf; ++j) {
    if (g_config.useforce[j])
      use_forces = 1;
#if defined(STRESS)
    if (g_config.usestress[j])
      use_stress = 1;
    if (use_forces && use_stress)
#else
    if (use_forces)
#endif
      break;
  }

  int numberOfComputeArgumentNames = 0;
  KIM_ComputeArgumentName computeArgumentName;
  KIM_SupportStatus supportStatus;

  KIM_COMPUTE_ARGUMENT_NAME_GetNumberOfComputeArgumentNames(&numberOfComputeArgumentNames);
  for (int i = 0; i < numberOfComputeArgumentNames; ++i) {
    res = KIM_COMPUTE_ARGUMENT_NAME_GetComputeArgumentName(i, &computeArgumentName);
    if (res)
      error(1, "Error getting argument name %d", i);
    res = KIM_ComputeArguments_GetArgumentSupportStatus(g_kim.arguments[0], computeArgumentName, &supportStatus);
    if (res)
      error(1, "Error getting argument supportStatus %d", i);

    // potfit requires energy
    if (KIM_ComputeArgumentName_Equal(computeArgumentName, KIM_COMPUTE_ARGUMENT_NAME_partialEnergy)) {
      if (KIM_SupportStatus_Equal(supportStatus, KIM_SUPPORT_STATUS_notSupported))
        error(1, "Selected KIM model does not support energy calculation");
      continue;
    }

    if (use_forces) {
      if (KIM_ComputeArgumentName_Equal(computeArgumentName, KIM_COMPUTE_ARGUMENT_NAME_partialForces)) {
        if (KIM_SupportStatus_Equal(supportStatus, KIM_SUPPORT_STATUS_notSupported))
          error(1, "Selected KIM model does not support force calculation");
        continue;
      }
    }

#if defined(STRESS)
    if (use_stress) {
      if (KIM_ComputeArgumentName_Equal(computeArgumentName, KIM_COMPUTE_ARGUMENT_NAME_partialVirial)) {
        if (KIM_SupportStatus_Equal(supportStatus, KIM_SUPPORT_STATUS_notSupported))
          error(1, "Selected KIM model does not support stress calculation");
        continue;
      }
    }
#endif

    if (KIM_SupportStatus_Equal(supportStatus, KIM_SUPPORT_STATUS_required))
     error(1, "Selected KIM model has an unsupported required argument %s", KIM_ComputeArgumentName_ToString(computeArgumentName));
  }
}

/*******************************************************************************
    initialize_KIM
*******************************************************************************/

void initialize_KIM()
{
  printf("\nInitializing KIM interface ...\n");

  int res = 0;

  // Allocate memory for KIM compute arguments
  g_kim.arguments = (KIM_ComputeArguments**)Malloc(g_config.nconf * sizeof(KIM_ComputeArguments*));
  g_kim.helpers = (potfit_compute_helper_t*)Malloc(g_config.nconf * sizeof(potfit_compute_helper_t));

  // TODO check for MPI compatibility
  for (int i = 0; i < g_config.nconf; i++) {
    res = KIM_Model_ComputeArgumentsCreate(g_kim.model, &g_kim.arguments[i]);
    if (res)
      error(1, "Error calling KIM_Model_ComputeArgumentsCreate: %d\n", res);

    if (i == 0)
      check_KIM_model_support();

    // check arguments
    int numberOfComputeCallbackNames;
    KIM_ComputeCallbackName computeCallbackName;
    KIM_SupportStatus supportStatus;

    // check call backs
    KIM_COMPUTE_CALLBACK_NAME_GetNumberOfComputeCallbackNames(&numberOfComputeCallbackNames);
    for (int j = 0; j < numberOfComputeCallbackNames; ++j) {
      res = KIM_COMPUTE_CALLBACK_NAME_GetComputeCallbackName(j, &computeCallbackName);
      if (res)
        error(1, "Error getting ComputeCallbackName %d: %d\n", j, res);

      res = KIM_ComputeArguments_GetCallbackSupportStatus(g_kim.arguments[i], computeCallbackName, &supportStatus);
      if (res)
        error(1, "Error getting CallbackSupportStatus for callback %d: %d\n", j + 1, res);

      if (KIM_ComputeCallbackName_Equal(computeCallbackName, KIM_COMPUTE_CALLBACK_NAME_GetNeighborList)) {
        res = KIM_ComputeArguments_SetCallbackPointer(g_kim.arguments[i], KIM_COMPUTE_CALLBACK_NAME_GetNeighborList, KIM_LANGUAGE_NAME_c, (KIM_Function*)get_neigh, &g_kim.helpers[i]);
        if (res)
          error(1, "Error setting GetNeighborList callback function: %d\n", res);
        continue;
      }

      if (KIM_ComputeCallbackName_Equal(computeCallbackName, KIM_COMPUTE_CALLBACK_NAME_ProcessDEDrTerm)) {
        if (KIM_SupportStatus_Equal(supportStatus, KIM_SUPPORT_STATUS_required) || KIM_SupportStatus_Equal(supportStatus, KIM_SUPPORT_STATUS_optional)) {
          res = KIM_ComputeArguments_SetCallbackPointer(g_kim.arguments[i], KIM_COMPUTE_CALLBACK_NAME_ProcessDEDrTerm, KIM_LANGUAGE_NAME_c, (KIM_Function*)process_DEDr, NULL);
          if (res)
            error(1, "Error setting ProcessDEDrTerm callback function: %d\n", res);
        }
        continue;
      }

      if (KIM_ComputeCallbackName_Equal(computeCallbackName, KIM_COMPUTE_CALLBACK_NAME_ProcessD2EDr2Term)) {
        if (KIM_SupportStatus_Equal(supportStatus, KIM_SUPPORT_STATUS_required) || KIM_SupportStatus_Equal(supportStatus, KIM_SUPPORT_STATUS_optional)) {
          res = KIM_ComputeArguments_SetCallbackPointer(g_kim.arguments[i], KIM_COMPUTE_CALLBACK_NAME_ProcessD2EDr2Term, KIM_LANGUAGE_NAME_c, (KIM_Function*)process_D2EDr2, NULL);
          if (res)
            error(1, "Error setting ProcessD2EDr2Term callback function: %d\n", res);
        }
        continue;
      }

      // cannot handle any other "required" call backs
      if (KIM_SupportStatus_Equal(supportStatus, KIM_SUPPORT_STATUS_required))
        error(1, "KIM model requires the following unsupported callback: %s\n", KIM_ComputeCallbackName_ToString(computeCallbackName));
    }

    res = KIM_ComputeArguments_SetArgumentPointerInteger(g_kim.arguments[i], KIM_COMPUTE_ARGUMENT_NAME_numberOfParticles, g_config.number_of_particles + i);
    if (res)
      error(1, "Error setting numberOfParticles");

    res = KIM_ComputeArguments_SetArgumentPointerInteger(g_kim.arguments[i], KIM_COMPUTE_ARGUMENT_NAME_particleSpeciesCodes, g_config.species_codes[i]);
    if (res)
      error(1, "Error setting particleSpeciesCodes");

    res = KIM_ComputeArguments_SetArgumentPointerInteger(g_kim.arguments[i], KIM_COMPUTE_ARGUMENT_NAME_particleContributing, g_config.particle_contributing[i]);
    if (res)
      error(1, "Error setting particleContributing");

    res = KIM_ComputeArguments_SetArgumentPointerDouble(g_kim.arguments[i], KIM_COMPUTE_ARGUMENT_NAME_coordinates, g_config.coordinates[i]);
    if (res)
      error(1, "Error setting coordinates");

    g_kim.helpers[i].config = i;

    g_kim.helpers[i].forces = (double*)Malloc(3 * g_config.number_of_particles[i] * sizeof(double));

    res = KIM_ComputeArguments_SetArgumentPointerDouble(g_kim.arguments[i], KIM_COMPUTE_ARGUMENT_NAME_partialForces, g_kim.helpers[i].forces);
    if (res)
      error(1, "Error setting forces!\n");
  }

  printf("Initializing KIM interface ... done\n");
}

/******************************************************************************
    shutdown_KIM
******************************************************************************/

void shutdown_KIM()
{
  if (g_kim.arguments) {
    for(int i = 0; i < g_config.nconf; i++) {
      if (g_kim.arguments[i])
        KIM_Model_ComputeArgumentsDestroy(g_kim.model, &g_kim.arguments[i]);
    }
  }
  if (g_kim.model)
    KIM_Model_Destroy(&g_kim.model);
}
