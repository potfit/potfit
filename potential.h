/****************************************************************
 *
 * potential.h: potfit potential functions header file
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

#ifndef POTENTIAL_H
#define POTENTIAL_H

/* reading the potential file */
void  read_pot_table(pot_table_t *, char *);
#ifdef APOT
void  read_pot_table0(pot_table_t *, apot_table_t *, char *, FILE *);
#else
void  read_pot_table3(pot_table_t *, int, char *, FILE *);
void  read_pot_table4(pot_table_t *, int, char *, FILE *);
#endif /* APOT */

/* calculating the potential tables */
void  init_calc_table(pot_table_t *, pot_table_t *);
#ifdef APOT
void  update_apot_table(double *);
void  update_calc_table(double *, double *, int);
#endif /* APOT */

/* writing potentials to files */
#ifdef APOT
void  write_pot_table0(apot_table_t *, char *);
#else
void  write_pot_table3(pot_table_t *, char *);
void  write_pot_table4(pot_table_t *, char *);
#endif /* APOT */
void  write_pot_table_imd(pot_table_t *, char *);
void  write_pot_table_lammps(pot_table_t *);
void  write_plotpot_pair(pot_table_t *, char *);
void  write_altplot_pair(pot_table_t *, char *);
#ifdef PDIST
void  write_pairdist(pot_table_t *, char *);
#endif /* PDIST */

/* functions for electrostatic calculations  */
#ifdef COULOMB
void  init_tails(double);
void  write_coulomb_table(void);
#endif /* COULOMB */

#ifdef APOT
void  write_imd_data_pair(FILE *, char *, int, int);
#endif /* APOT */

#endif /* POTENTIAL_H */
