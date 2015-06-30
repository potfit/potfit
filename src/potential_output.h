/****************************************************************
 *
 * potential.h: potfit potential functions header file
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

#ifndef POTENTIAL_H
#define POTENTIAL_H

void write_pot_table_potfit(char const *filename);
void write_pot_table_imd(char const *prefix);
void write_pot_table_lammps();

void write_plotpot_pair(pot_table_t *, char *);
void write_altplot_pair(pot_table_t *, char *);
#ifdef PDIST
void write_pairdist(pot_table_t *, char *);
#endif /* PDIST */

/* functions for electrostatic calculations  */
#ifdef COULOMB
void write_coulomb_table(void);
#endif /* COULOMB */

#ifdef APOT
void write_imd_data_pair(FILE *, char *, int, int);
#endif /* APOT */

#endif /* POTENTIAL_H */
