/****************************************************************
 *
 * potential.h: potfit potential functions header file
 *
 ****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
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

#ifndef POTENTIAL_OUTPUT_H_INCLUDED
#define POTENTIAL_OUTPUT_H_INCLUDED

void write_pot_table_potfit(char const* filename);
void write_pot_table_imd(char const* prefix);
void write_pot_table_lammps();

void write_plotpot_pair(pot_table_t*, const char*);
void write_altplot_pair(pot_table_t*, const char*);
#if defined(PDIST)
void write_pairdist(pot_table_t*, const char*);
#endif  // PDIST

/* functions for electrostatic calculations  */
#if defined(COULOMB)
void write_coulomb_table(void);
#endif  // COULOMB

#if defined(APOT)
void write_imd_data_pair(FILE*, char*, int, int);
#endif  // APOT

#endif  // POTENTIAL_OUTPUT_H_INCLUDED
