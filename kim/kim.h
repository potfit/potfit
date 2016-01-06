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

#include "KIM_API_C.h"
#include "KIM_API_status.h"
#include "../potfit.h"

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


/* global variables */
/******************************************************************************/ 
#ifndef EXTERN
#ifdef KIM_MAIN
#define EXTERN          /* define variables in main */
#else
#define EXTERN extern   /* declare them extern otherwise */
#endif /* KIM_MAIN */
#endif /* EXTERN */

EXTERN void** pkimObj;  /* pointers to kim objects */

FreeParamType* FreeParamAllConfig; /* free param struct for all configurations */

EXTERN char kim_model_name[255];  /* kim model name (read in from input )*/

EXTERN char** name_opt_param;  /* optimizable parameter names (read in from input) */

EXTERN int num_opt_param;  /* number of optimizalbe params( read in from input) */

EXTERN int* size_opt_param;  /* size of each parameter (number of values each
                              * parameter name represents) */

EXTERN int is_half_neighbors;  /* using half neighbor list? 1 = half, 0 = full */

EXTERN double* box_side_len;  /*box_side_len is used to enable MI_OPBC_H in KIM. 
                               * box_side_len is used to store the info of box
                               * size of each configuration. box_side_len[0], 
                               * box_side_len[1], box_side_len[2] store box_x.x,
                               * box_y.y box_z.z of the first configuration,
                               * respectively. box_side_len[3] stores box_x.x
                               * of the second config...  Note this is only 
                               * possible when box_x.y = 0, and similar for other
                               * components of the box. If not, we need to come 
                               * up with other methods.*/ 


/* function prototypes */
/******************************************************************************/ 

/* called in `potfit.c' */
void init_KIM();

void free_KIM();

/* called by `init_KIM' */
int write_descriptor_file(int Nspecies, char** species);

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
int read_potential_keyword(pot_table_t* pt, char* filename, FILE* infile,
								    			FreeParamType* FreeParam);

/* called in `potential_input.c' and by `read_potential_keyword'  `init_optimizable_param' */
int write_temporary_descriptor_file(char* modelname);

int get_free_param_double(void* pkim, FreeParamType* FreeParam);

int nest_optimizable_param(void* pkim, FreeParamType* FreeParam, 
                              char** input_param_name, int input_param_num);

/* called by `init_optimizable_param' */
int publish_cutoff(void* pkim, double cutoff);

/* called in `force_kim.c' */
int calc_force_KIM(void* pkim, double** energy, double** force, double** virial,
							int useforce, int usestress);

int publish_param(void* pkim, FreeParamType* FreeParam, double* PotTable);

/* used in potential.c */
void read_pot_table5_no_nolimits(pot_table_t *pt, apot_table_t *apt, char *filename, FILE *infile);
void read_pot_table5_with_nolimits(pot_table_t *pt, int size, char *filename, FILE *infile);

/* assigned to function pointer `write_pot_table' in `potfit.c' */
void write_pot_table5(pot_table_t *pt, char *filename);


/* free memory */
int free_model_object(void** pkim);

void free_KIM();

#ifdef PAIR
/* used only for check purpose (could be deleted ) */
int AnalyticForce(double epsilon, double sigma, double cutoff,
									double r, double* phi_val, double* phi_grad); 
#endif /* PAIR */


#endif /* KIM_H */
