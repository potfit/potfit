/****************************************************************
 *
 * params.h: potfit routines for reading the parameter file
 *
 ****************************************************************
 *
 * Copyright 2002-2015
 *      Institute for Theoretical and Applied Physics
 *      University of Stuttgart, D-70550 Stuttgart, Germany
 *      http://potfit.sourceforge.net/
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

#ifndef POTFIT_PARAMS_H
#define POTFIT_PARAMS_H

#include <float.h>
#include <limits.h>

void read_parameters(int argc, char** argv);

int read_parameter_file(char const* param_file);

int get_param_double(
  char const* param_name,
  double* value,
  int line,
  char const* param_file,
  double min,
  double max);

int get_param_int(
  char const* param_name,
  int* value,
  int line,
  char const* param_file,
  int min,
  int max);

int get_param_string(
  char const* param_name,
  char** value,
  int line,
  char const* param_file);

void  check_parameters_complete(char const* paramfile);

#endif // POTFIT_PARAMS_H
