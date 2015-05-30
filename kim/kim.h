/****************************************************************
*
* KIM-potfit
* 
* force_kim.h
*
* header file for calculating force via KIM
****************************************************************/ 

#ifndef KIM_H
#define KIM_H

#include "KIM_API_C.h"
#include "KIM_API_status.h"

#ifndef DIM
#define DIM 3
#endif



/* Define neighborlist structure */
typedef struct
{
   int iteratorId;
   int* NNeighbors;
   int* neighborList;
   double* RijList;
	 int* BeginIdx;   /* The position of first neighbor of each atom in 
	 										the neighbor list */
} NeighObjectType;

typedef struct
{
  char** name;
 	double** value;       /* the pointer to parameters */
	double** nestedvalue; /* nest all values(pointer) into one long variable */
  int* rank;
	int** shape;
	int Nparam;				    /* number of optimizable parameters */
  int Nnestedvalue;     /* the length of nestedvalue, it would be equal to
                         * (Nparam - 1) if all the parameters are scalar.
                         * Otherwise, it is larger. (Nparam - 1), because
                         * PARAM_FREE_cutoff is counted in Nparam, but not in
                         * Nnestedvalue. */ 
} FreeParamType;

FreeParamType* FreeParamAllConfig;


/* function prototypes */

/* functions defined in `kim.c' */
void InitKIM();

void FreeKIM();

int InitObject();

int get_FreeParamDouble(void* pkim, FreeParamType* FreeParam);

int nest_OptimizableParamValue(void* pkim, FreeParamType* FreeParam, 
                              char** input_param_name, int input_param_num);

int get_OptimizableParamSize(FreeParamType* FreeParam,
                              char** input_param_name, int input_param_num);

int PublishCutoff(void* pkim, double cutoff);

int CreateKIMObj(void* pkim, int Natoms, int Nspecies, int start);

int get_neigh(void* kimmdl, int* mode, int *request, int* part,
              int* numnei, int** nei1part, double** Rij);					

int CalcForce(void* pkim, double** energy, double** force, double** virial,
							int useforce, int usestress);

/* functions in force_[interaction]_kim.c */
int PublishParam(void* pkim, FreeParamType* FreeParam, double* PotTable);

/* read potential input file */
int ReadPotentialKeywords(pot_table_t* pt, char* filename, FILE* infile,
								    			FreeParamType* FreeParam);
/* write potential file */
void write_pot_table5(pot_table_t *pt, char *filename);

#ifdef PAIR
/* used only for check purpose (could be deleted ) */
int AnalyticForce(double epsilon, double sigma, double cutoff,
									double r, double* phi_val, double* phi_grad); 
#endif /* PAIR */

#ifdef EAM
int CreateModel();
int MakeModel();
int CreateParamFileEAM();
#endif /*EAM*/

#endif /* KIM_H */
