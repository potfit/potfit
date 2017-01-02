/****************************************************************
 *
 * kim.h: header file for KIM interface
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

#ifndef KIM_H_INCLUDED
#define KIM_H_INCLUDED

#if defined(KIM)

#include "KIM_API_C.h"
#include "KIM_API_status.h"

#include "kim.h"

#if !defined(DIM)
#define DIM 3
#else
STATIC_ASSERT(DIM==3, potfit_kim_support_requires_DIM_eq_3);
#endif

// called from potfit.c
void initialize_KIM();
void shutdown_KIM();

/* called in `potential_input.c' */
int read_potential_keyword(pot_table_t* pt, char const* filename, FILE* infile);

/* called in `potential_input.c' and by `read_potential_keyword' */
int write_temporary_descriptor_file(char* modelname);

/* called by `init_object' */
int setup_KIM_API_object(void** pkim, int Natoms, int Nspecies, char* modelname);

/* set compute flags and get NBC method, called in `potential_input.c' */
void get_compute_const(void* pkim);

/* free memory */
int free_model_object(void** pkim);

int nest_optimizable_param(void* pkim, FreeParamType* FreeParam,
                              char** input_param_name, int input_param_num);


int get_free_param_double(void* pkim, FreeParamType* FreeParam);

int publish_param(void* pkim, FreeParamType* FreeParam, double* PotTable);

/* called in `force_kim.c' */
int calc_force_KIM(void* pkim, double** energy, double** force, double** virial,
							int useforce, int usestress);

#endif // KIM

#endif // KIM_H_INCLUDED
