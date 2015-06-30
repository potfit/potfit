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

#ifndef CONFIG_H
#define CONFIG_H

#ifndef POTFIT_H
#include "potfit.h"
#endif /* POTFIT_H */

typedef struct {
  int atom_count;
  int use_force;
  int line;
  int config;
  int have_elements;
  int have_energy;
  int have_stress;
  int num_fixed_elements;
  vector box_x, box_y, box_z;
  vector tbox_x, tbox_y, tbox_z;
  vector cell_scale;
  int h_box_x, h_box_y, h_box_z;
  sym_tens* stresses;
  int have_contrib_box;
  vector cbox_o;                 /* origin of box of contrib. atoms */
  vector cbox_a, cbox_b, cbox_c; /* box vectors for box of contrib. atoms */
  vector* sphere_centers;        /* centers of the spheres of contrib. atoms */
} config_state;

void read_config(char const*);

void create_memory_for_config(config_state* cstate);
void read_chemical_elements(char* psrc, config_state* cstate);

void init_atom_memory(atom_t* atom);
void init_neigh_memory(neigh_t* neighbor);
#if defined(THREEBODY)
void init_angle_memory(angle_t* angle);
#endif /* THREEBODY */

void init_neighbors(config_state* cstate, double* mindist);
void init_angles(config_state* cstate);

void set_neighbor_slot(neigh_t* neighbor, int col, double r, int neighbor_slot);

double make_box(config_state* cstate);

void write_pair_distribution_file();
void print_minimal_distances_matrix(double const* mindist);

#ifdef CONTRIB
int does_contribute(vector);
#endif /* CONTRIB */

#ifdef APOT
void update_slots(void);
void update_neighbor_slots(neigh_t* neighbor, double r, int neighbor_slot);
#endif /* APOT */

#endif /* CONFIG_H */
