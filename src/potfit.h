/****************************************************************
 *
 * potfit.h: main potfit header file
 *
 ****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
 *
 * http://potfit.sourceforge.net/
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

#ifndef POTFIT_H_INCLUDED
#define POTFIT_H_INCLUDED

#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(MPI)
#include <mpi.h>
#endif  // MPI

#define POTFIT_VERSION "potfit-git"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define SPROD(a, b) (((a).x * (b).x) + ((a).y * (b).y) + ((a).z * (b).z))
#define SWAP(A, B, C) \
  (C) = (A);          \
  (A) = (B);          \
  (B) = (C);

// include preprocessor flags and certain compile time constants
#include "defines.h"

// include type definitions after all preprocessor flags are properly set
#include "types.h"

// include kim types
#if defined(KIM)
#include "kim.h"
#endif

// global variables (defined in potfit.c)

extern potfit_calculation g_calc;
extern potfit_configurations g_config;
extern potfit_filenames g_files;
extern potfit_mpi_config g_mpi;
extern potfit_parameters g_param;
extern potfit_potentials g_pot;
#if defined(KIM)
extern potfit_kim g_kim;
#endif


// general functions for warning and error output

void error(int done, const char* msg, ...);
void warning(const char* msg, ...);

#endif  // POTFIT_H_INCLUDED
