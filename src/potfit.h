/****************************************************************
 *
 * potfit.h: main potfit header file
 *
 ****************************************************************
 *
 * Copyright 2002-2015
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

#ifndef POTFIT_H
#define POTFIT_H

#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(MPI)
#include "mpi_utils.h"
#endif /* MPI */

#include "random.h"

#define POTFIT_VERSION "potfit-git"

#define MAX(a,b)   ((a) > (b) ? (a) : (b))
#define MIN(a,b)   ((a) < (b) ? (a) : (b))
#define SPROD(a,b) (((a).x * (b).x) + ((a).y * (b).y) + ((a).z * (b).z))
#define SWAP(A,B,C) (C)=(A);(A)=(B);(B)=(C);

/****************************************************************
 *
 *  include preprocessor flags and certain compile time constants
 *
 ****************************************************************/

#include "defines.h"

/****************************************************************
 *
 *  include type definitions after all preprocessor flags are properly set
 *
 ****************************************************************/

#include "types.h"

/****************************************************************
 *
 *  global variables
 *
 ****************************************************************/

static potfit_calculation       g_calc;
static potfit_configurations    g_config;
static potfit_filenames         g_files;
static potfit_mpi_config        g_mpi;
static potfit_parameters        g_param;
static potfit_potentials        g_pot;
static potfit_memory            g_memory;
static potfit_unknown           g_todo;

/****************************************************************
 *
 *  general functions for setting up and terminating a potfit run
 *
 ****************************************************************/

void  read_input_files(int argc, char**  argv);

void  error(int, const char *, ...);
void  warning(const char *, ...);

/* error reports [errors.c] */
void  write_errors(double *, double);

/* force routines for different potential models [force_xxx.c] */
#ifdef PAIR
static const char interaction_name[5] = "PAIR";
#elif defined EAM && !defined COULOMB
#ifndef TBEAM
static const char interaction_name[4] = "EAM";
#else
static const char interaction_name[6] = "TBEAM";
#endif /* TBEAM */
#elif defined ADP
static const char interaction_name[4] = "ADP";
#elif defined COULOMB && !defined EAM
static const char interaction_name[7] = "ELSTAT";
#elif defined COULOMB && defined EAM
static const char interaction_name[11] = "EAM_ELSTAT";
#elif defined MEAM
static const char interaction_name[5] = "MEAM";
#elif defined STIWEB
static const char interaction_name[7] = "STIWEB";
void  update_stiweb_pointers(double *);
#elif defined TERSOFF
#ifdef TERSOFFMOD
static const char interaction_name[11] = "TERSOFFMOD";
#else
static const char interaction_name[8] = "TERSOFF";
#endif /* TERSOFFMOD */
void  update_tersoff_pointers(double *);
#endif /* interaction type */

/* rescaling functions for EAM [rescale.c] */
#if !defined APOT && ( defined EAM || defined(ADP) || defined MEAM )
double rescale(pot_table_t *, double, int);
void  embed_shift(pot_table_t *);
#endif /* !APOT && (EAM || MEAM) */

#endif /* POTFIT_H */
