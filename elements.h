/****************************************************************
 *
 * elements.h: header file for periodic table
 *
 ****************************************************************
 *
 * Copyright 2012
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.itap.physik.uni-stuttgart.de/
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
#include "potfit.h"
#endif

typedef struct {
  char  name[20];
  char  short_name[20];
  double  mass;
} element_t;

/* initialize periodic table of elements */
void  init_elements();

/* get mass from number */
double  ele_mass_from_number(int num);

/* get mass from name */
double  ele_mass_from_name(char *name);

/* get number from name */
int   ele_number_from_name(char *name);
