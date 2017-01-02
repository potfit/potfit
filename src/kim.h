/******************************************************************************
*
* potfit-KIM
*
* kim.h
*
* header file for all KIM stuff
*
* Author: Mingjian Wen (wenxx151@umn.edu), University of Minnesota
******************************************************************************/

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

/* function prototypes */
/******************************************************************************/

// called from potfit.c
void init_KIM();
void free_KIM();

/* called by `init_KIM' */
int write_descriptor_file(int Nspecies, const char** species, int compute_energy,
                          int compute_forces, int compute_virial);

void init_object();

void init_optimizable_param();

/* called by `init_object' */
int setup_KIM_API_object(void** pkim, int Natoms, int Nspecies, char* modelname);

int init_KIM_API_argument(void* pkim, int Natoms, int Nspecies, int start);

/* function pointer assigned in `setup_neighborlist_KIM_access' */
int get_neigh(void* kimmdl, int* mode, int *request, int* part,
              int* numnei, int** nei1part, double** Rij);

/* called in `potential_input.c' */
int read_potential_keyword(pot_table_t* pt, char const* filename, FILE* infile);

/* called in `potential_input.c' and by `read_potential_keyword' */
int write_temporary_descriptor_file(char* modelname);

/* called by `init_obj' */
int write_final_descriptor_file(int u_f, int u_s);

int get_free_param_double(void* pkim, FreeParamType* FreeParam);

int nest_optimizable_param(void* pkim, FreeParamType* FreeParam,
                              char** input_param_name, int input_param_num);

/* called by `init_KIM' */
int publish_cutoff(void* pkim, double cutoff);

/* called in `force_kim.c' */
int calc_force_KIM(void* pkim, double** energy, double** force, double** virial,
							int useforce, int usestress);

int publish_param(void* pkim, FreeParamType* FreeParam, double* PotTable);

/* set compute flags and get NBC method, called in `potential_input.c' */
void get_compute_const(void* pkim);

int get_KIM_model_has_flags();

/* free memory */
int free_model_object(void** pkim);

#endif // KIM

#endif // KIM_H_INCLUDED
