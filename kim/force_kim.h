/****************************************************************
*
* KIM-potfit
* 
* force_kim.h
*
* header file for calculating force via KIM
****************************************************************/ 

#ifndef FORCEKIM_H
#define FORCEKIM_H

#include "KIM_API_C.h"
#include "KIM_API_status.h"

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

/* function prototypes */
int CreateKIMObj(void* pkim, int Natoms, int Nspecies, int start);

int get_neigh(void* kimmdl, int *mode, int *request, int* part,
              int* numnei, int** nei1part, double** Rij);					

int CalcForce(void* pkim, double** energy, double** force, double** virial);

#ifdef PAIR
int PublishParam(void* km, double* PotTable, double CutOff);
/* used only for check purpose */
int AnalyticForce(double epsilon, double sigma, double cutoff,
									double r, double* phi_val, double* phi_grad); 
#endif /* PAIR */

#ifdef EAM
int PublishParam(void* km, double* PotTable);
void* GetShapeAndPointer(void* const pkim, char const* const argName, 
												 int* const shape, int* const status);
#endif /*EAM*/

#endif /* FORCEKIM_H */
