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
void InitKIM();

void FreeKIM();

/* called by `InitKIM' */
void InitObject();

void InitOptimizableParam();

int PublishCutoff(void* pkim, double cutoff);

/* called by `InitObject' */
int CreateKIMObj(void* pkim, int Natoms, int Nspecies, int start);

/* called in `potential_input.c' */
int get_OptimizableParamSize(FreeParamType* FreeParam,
                              char** input_param_name, int input_param_num);

int ReadPotentialKeywords(pot_table_t* pt, char* filename, FILE* infile,
								    			FreeParamType* FreeParam);

/* called by `get_OptimizableParamSize' and `InitOptimizableParam' */
int get_FreeParamDouble(void* pkim, FreeParamType* FreeParam);

int nest_OptimizableParamValue(void* pkim, FreeParamType* FreeParam, 
                              char** input_param_name, int input_param_num);

/* called in `force_kim.c' */
int get_neigh(void* kimmdl, int* mode, int *request, int* part,
              int* numnei, int** nei1part, double** Rij);					

int CalcForce(void* pkim, double** energy, double** force, double** virial,
							int useforce, int usestress);

int PublishParam(void* pkim, FreeParamType* FreeParam, double* PotTable);

/* assigned to function pointer `write_pot_table' in `potfit.c' */
void write_pot_table5(pot_table_t *pt, char *filename);

#ifdef PAIR
/* used only for check purpose (could be deleted ) */
int AnalyticForce(double epsilon, double sigma, double cutoff,
									double r, double* phi_val, double* phi_grad); 
#endif /* PAIR */


#endif /* KIM_H */
