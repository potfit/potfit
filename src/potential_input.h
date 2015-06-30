/****************************************************************
 *
 * potential_input.h: potfit potential functions header file
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

#ifndef POTENTIAL_INPUT_H
#define POTENTIAL_INPUT_H

typedef struct {
  char const* filename;
  int end_header;
  int have_format;
  int have_gradient;
  int num_pots;
} potential_state;

/* reading the potential file */
void read_pot_table(char const* potential_filename);

void read_pot_table0(char const* potential_filename, FILE* pfile);
void read_pot_table3(char const* potential_filename, FILE* pfile,
                     potential_state* pstate);
void read_pot_table4(char const* potential_filename, FILE* pfile,
                     potential_state* pstate);

#ifdef APOT
void update_apot_table(double*);
void update_calc_table(double*, double*, int);
#endif /* APOT */

#endif /* POTENTIAL_INPUT_H */
