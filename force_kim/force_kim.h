/****************************************************************
 *
 * config.h: header file for reading atomic configurations
 *
 ****************************************************************
 *
 * Copyright 2002-2014
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.sourceforge.net/
 *
 ****************************************************************
 *
 *   This file is part of potfit.
 *
 *   potfit is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   potfit is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
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


#endif /* FORCEKIM_H */
