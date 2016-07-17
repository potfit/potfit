/******************************************************************************
*
* potfit-KIM
* 
* kim.h
*
* header file for all KIM stuff  
******************************************************************************/ 

#ifndef KIM_H
#define KIM_H

#include <stdio.h> 
#include "KIM_API_C.h"
#include "KIM_API_status.h"
/*#include "potfit.h"
*/


#ifndef DIM
#define DIM 3
#endif


/* types */
/******************************************************************************/ 
/* neighborlist struct */
typedef struct
{
   int iteratorId;
   int* NNeighbors;
   int* neighborList;
   double* RijList;
	 int* BeginIdx;       /* The position of first neighbor of each atom in 
	 									     * the neighbor list */
} NeighObjectType;

/* sturct of all free parameters of data type `double' */
typedef struct
{
  char** name;
 	double** value;       /* the pointers to parameters */
  int* rank;
	int** shape;
	int Nparam;				    /* number of free parameters */
	double** nestedvalue; /* nest all values(pointers) of optimizable parameters
                         * into this long list */
  int Nnestedvalue;     /* the length of nestedvalue */
} FreeParamType;

/* varialbes */
typedef struct
{
  void** pkimObj;                    /* pointers to kim objects */
  FreeParamType* FreeParamAllConfig; /* free param struct for all configurations */
  char kim_model_name[255]; /* kim model name (read in from input )*/
  char NBC_method[64];      /* neighbor list and boundary conditions */
  char** name_opt_param;    /* optimizable parameter names (read in from input) */
  int kim_model_has_forces;
  int kim_model_has_virial;
  int num_elements;         /* number of elements read in from config file */
  int is_half_neighbors;    /* using half neighbor list? 1 = half, 0 = full */
  int kim_model_has_energy; /* flag, whether KIM model has the routine to 
                             * compute energy, forces, and virial */
  int num_opt_param;        /* number of optimizalbe params( read in from input) */
  int* size_opt_param;      /* size of each parameter (number of values each
                             * parameter name represents) */

  double* box_side_len;     /* box_side_len is used to enable MI_OPBC in KIM. */ 

} potfit_kim;

/* function prototypes */
/******************************************************************************/ 

/* called in `potfit.c' */
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

int setup_neighborlist_KIM_access(void* pkim, NeighObjectType* NeighObject); 

int init_neighborlist(NeighObjectType* NeighObject, int Natoms, int start);

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

/* used in potential.c */
/*#ifndef NOLIMITS
void read_pot_table5_no_nolimits(pot_table_t *pt, apot_table_t *apt, char *filename, FILE *infile);
#else
void read_pot_table5_with_nolimits(pot_table_t *pt, int size, char *filename, FILE *infile);
#endif
*/



/* assigned to function pointer `write_pot_table' in `potfit.c' */
/*void write_pot_table5(pot_table_t *pt, char *filename);
*/


/* set compute flags and get NBC method, called in `potential_input.c' */
void get_compute_const(void* pkim); 

int get_KIM_model_has_flags();


/* free memory */
int free_model_object(void** pkim);

void free_KIM();


#endif /* KIM_H */
