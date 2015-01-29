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
int CreateKIMObj(void* km, int Natoms, int Nspecies, int start);

int get_neigh(void* kimmdl, int *mode, int *request, int* part,
              int* numnei, int** nei1part, double** Rij);

int PublishParam(void* km, double* PotTable, double CutOff);							

int CalcForce(void* km, double** energy, double** force, double** virial);

/* used only for check purpose */
int AnalyticForce(double epsilon, double sigma, double cutoff,
									double r, double* phi_val, double* phi_grad); 
#endif /* FORCEKIM_H */
