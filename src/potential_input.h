/****************************************************************
 *
 * potential_input.h: potfit potential functions header file
 *
 ****************************************************************
 *
 * Copyright 2002-2017 - the potfit development team
 *
 * https://www.potfit.net/
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

#ifndef POTENTIAL_INPUT_H_INCLUDED
#define POTENTIAL_INPUT_H_INCLUDED

typedef struct {
  char const* filename;
  int end_header;
  int have_format;
  int have_gradient;
  int num_pots;
} potential_state;

// reading the potential file
void read_pot_table(char const* potential_filename);

#if defined(APOT)
void update_apot_table(double* xi);
void update_calc_table(double* xi_opt, double* xi_calc, int do_all);
#endif  // APOT

#endif  // POTENTIAL_INPUT_H_INCLUDED
