/****************************************************************
 *
 * config.c: Reads atomic configurations and forces.
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

#include <float.h>

#include "potfit.h"

#include "config.h"
#include "memory.h"
#include "utils.h"

typedef struct {
  const char* filename;
  int atom_count;
  int use_force;
  int line;
  int config;
  int have_energy;
  int have_stress;
  int num_fixed_elements;
  vector box_x, box_y, box_z;
  vector tbox_x, tbox_y, tbox_z;
  vector cell_scale;
  int have_box_vector;
  sym_tens* stresses;
#if defined(CONTRIB)
  int have_contrib_box_vector;
  int n_spheres;
  vector cbox_o;                 /* origin of box of contrib. atoms */
  vector cbox_a, cbox_b, cbox_c; /* box vectors for box of contrib. atoms */
  vector* sphere_center;         /* centers of the spheres of contrib. atoms */
  double* sphere_radius;
#endif  // CONTRIB
} config_state;

void reset_cstate(config_state* cstate);
void create_memory_for_config(config_state* cstate);
void init_atom_memory(atom_t* atom);
void read_box_vector(char const* pline, vector* pvect, config_state* cstate);
void read_chemical_elements(char* psrc, config_state* cstate);
void init_box_vectors(config_state* cstate);
void init_neighbors(config_state* cstate, double* mindist);
void set_neighbor_slot(neigh_t* neighbor, int col, double r, int neighbor_slot);
void init_angles(config_state* cstate);

double make_box(config_state* cstate);
void write_pair_distribution_file(void);
void print_minimal_distances_matrix(double const* mindist);

#if defined(CONTRIB)
void read_sphere_center(char const* pline, config_state* cstate);
int does_contribute(vector atom_pos, config_state* cstate);
#endif  // CONTRIB

/****************************************************************
  read_config
****************************************************************/

void read_config(const char* filename)
{
  config_state cstate;

  char buffer[1024];
  char* ptr = NULL;
  char* res = NULL;

  int max_atom_type = -1;
  int with_forces = 0;
  int with_stresses = 0;

  fpos_t filepos;

  memset(&cstate, 0, sizeof(cstate));

  cstate.filename = filename;

  // initialize elements array if not yet done from potential file
  if (g_config.elements == NULL) {
    g_config.elements = (char const**)Malloc(g_param.ntypes * sizeof(char*));
    for (int i = 0; i < g_param.ntypes; i++) {
      g_config.elements[i] = (char*)Malloc(5 * sizeof(char));
      sprintf((char*)g_config.elements[i], "%d", i);
      *((char*)g_config.elements[i] + (i < 10 ? 1 : 2)) = '\0';
    }
  } else
    cstate.num_fixed_elements = g_param.ntypes;

  // initialize minimum distance array
  double* mindist = (double*)Malloc(isquare(g_param.ntypes) * sizeof(double));

  // set maximum cutoff distance as starting value for mindist
  for (int i = 0; i < g_param.ntypes * g_param.ntypes; i++)
    mindist[i] = DBL_MAX;

  for (int i = 0; i < g_param.ntypes; i++) {
    for (int j = 0; j < g_param.ntypes; j++) {
      int k = (i <= j) ? i * g_param.ntypes + j - ((i * (i + 1)) / 2)
                       : j * g_param.ntypes + i - ((j * (j + 1)) / 2);
      mindist[k] = MAX(g_config.rcut[i * g_param.ntypes + j],
                       mindist[i * g_param.ntypes + j]);
    }
  }

  // open file
  FILE* config_file = fopen(filename, "r");

  if (config_file == NULL)
    error(1, "Could not open file %s\n", filename);

  printf(
      "Reading the config file >> %s << and calculating neighbor lists ...\n",
      filename);
  fflush(stdout);

  // read configurations until the end of the file
  do {
    res = fgets(buffer, 1024, config_file);

    ++cstate.line;

    if (res == NULL)
      error(1, "Unexpected end of file in %s\n", filename);

    if (res[0] == '#' && res[1] == 'N') {
      reset_cstate(&cstate);
      ++cstate.config;

      if (sscanf(res + 3, "%d %d", &cstate.atom_count, &cstate.use_force) < 2)
        error(1, "%s: Error in atom number specification on line %d\n",
              filename, cstate.line);
    } else
      continue;

// check if there are enough atoms, 2 for pair and 3 for manybody potentials
#if defined(THREEBODY)
    int min_atom_count = 3;
#else
    int min_atom_count = 2;
#endif  // THREEBODY
    if (cstate.atom_count < min_atom_count)
      error(1,
            "Configuration %d (starting on line %d) has not enough atoms, "
            "please remove it.\n",
            g_config.nconf + 1, cstate.line);

    create_memory_for_config(&cstate);

    for (int i = g_config.natoms; i < g_config.natoms + cstate.atom_count; i++)
      memset(g_config.atoms + i, 0, sizeof(atom_t));

    g_config.inconf[g_config.nconf] = cstate.atom_count;
    g_config.cnfstart[g_config.nconf] = g_config.natoms;
    g_config.useforce[g_config.nconf] = cstate.use_force;
#if defined(STRESS)
    cstate.stresses = g_config.stress + g_config.nconf;
#endif  // STRESS

#if defined(KIM)
    if (strcmp(g_kim.NBC_method, "MI_OPBC_F") == 0 || strcmp(g_kim.NBC_method, "MI_OPBC_H") == 0) {
      /* realloc memory for box_side_len */
      g_kim.box_side_len = (double *) Realloc(g_kim.box_side_len, 3*(g_config.nconf+1)*sizeof(double));
    }
#endif  // KIM

    // read header lines
    do {
      res = fgets(buffer, 1024, config_file);

      if (res == NULL || feof(config_file))
        error(1, "Incomplete header on line %d in configuration file %s\n",
              cstate.line, filename);

      if ((ptr = strchr(res, '\n')) != NULL)
        *ptr = '\0';

      cstate.line++;

      if (res[0] != '#') {
        warning("Ignoring unknown header item in line %d of file >> %s <<:\n",
                cstate.line, filename);
        warning("  \"%s\"\n", res);
      }

      switch (res[1]) {
        /* read the box vectors */
        case 'x':
        case 'X':
          read_box_vector(res + 3, &cstate.box_x, &cstate);
          cstate.have_box_vector |= 1 << 0;
          break;
        case 'y':
        case 'Y':
          read_box_vector(res + 3, &cstate.box_y, &cstate);
          cstate.have_box_vector |= 1 << 1;
          break;
        case 'z':
        case 'Z':
          read_box_vector(res + 3, &cstate.box_z, &cstate);
          cstate.have_box_vector |= 1 << 2;
          break;
#if defined(CONTRIB)
        case 'b':
        case 'B': {
          switch (res[3]) {
            case 'o':
            case 'O':
              if (cstate.have_contrib_box_vector & 1) {
                error(0, "There can only be one box of contributing atoms\n");
                error(1, "  This occured in %s on line %d\n", filename,
                      cstate.line);
              }
              read_box_vector(res + 5, &cstate.cbox_o, &cstate);
              cstate.have_contrib_box_vector |= 1 << 0;
              break;
            case 'a':
            case 'A':
              read_box_vector(res + 5, &cstate.cbox_a, &cstate);
              cstate.have_contrib_box_vector |= 1 << 1;
              break;
            case 'b':
            case 'B':
              read_box_vector(res + 5, &cstate.cbox_b, &cstate);
              cstate.have_contrib_box_vector |= 1 << 2;
              break;
            case 'c':
            case 'C':
              read_box_vector(res + 5, &cstate.cbox_c, &cstate);
              cstate.have_contrib_box_vector |= 1 << 3;
              break;
            case 's':
            case 'S':
              read_sphere_center(res + 5, &cstate);
              cstate.n_spheres++;
              break;
          }
        } break;
#endif  // CONTRIB
        case 'e':
        case 'E':
          if (sscanf(res + 3, "%lf\n", &(g_config.coheng[g_config.nconf])) == 1)
            cstate.have_energy = 1;
          else
            error(1, "%s: Error in energy on line %d\n", filename, cstate.line);
          break;
        case 'w':
        case 'W':
          if (sscanf(res + 3, "%lf\n",
                     &(g_config.conf_weight[g_config.nconf])) != 1)
            error(1, "%s: Error in configuration weight on line %d\n", filename,
                  cstate.line);
          if (g_config.conf_weight[g_config.nconf] < 0.0)
            error(1, "%s: The configuration weight is negative on line %d\n",
                  filename, cstate.line);
          break;
        case 'c':
        case 'C':
          fgetpos(config_file, &filepos);
          read_chemical_elements(res, &cstate);
          fsetpos(config_file, &filepos);
          break;
#if defined(STRESS)
        case 's':
        case 'S':
          if (sscanf(res + 3, "%lf %lf %lf %lf %lf %lf\n",
                     &(cstate.stresses->xx), &(cstate.stresses->yy),
                     &(cstate.stresses->zz), &(cstate.stresses->xy),
                     &(cstate.stresses->yz), &(cstate.stresses->zx)) == 6)
            cstate.have_stress = 1;
          else
            error(1, "Error in stress tensor on line %d\n", cstate.line);
          break;
#endif  // STRESS
        case 'f':
        case 'F':
          break;
        default:
          if (res[1] != '#') {
            warning(
                "Ignoring unknown header item in line %d of file >> %s <<:\n",
                cstate.line, filename);
            warning("  \"%s\"\n", res);
          }
          break;
      }

    } while (res[1] != 'F');

    if (cstate.have_energy == 0)
      error(1, "%s: missing energy in configuration %d!\n", filename,
            g_config.nconf + 1);

    if (cstate.have_box_vector != 7)
      error(1, "Incomplete box vectors for config %d!\n", g_config.nconf + 1);

#if defined(CONTRIB)
    if (cstate.have_contrib_box_vector && cstate.have_contrib_box_vector != 15)
      error(1, "Incomplete box of contributing atoms for config %d!\n",
            g_config.nconf + 1);
#endif  // CONTRIB

    if (cstate.use_force)
      with_forces++;

#if defined(STRESS)
    if (cstate.have_stress == 1) {
      with_stresses++;
      g_config.usestress[g_config.nconf] = 1;
    }
#endif  // STRESS

    g_config.volume[g_config.nconf] = make_box(&cstate);

#if defined(KIM)
    if (strcmp(g_kim.NBC_method, "MI_OPBC_F") == 0 || strcmp(g_kim.NBC_method, "MI_OPBC_H") == 0) {
      double small_value = 1e-8;
      if(cstate.box_x.y > small_value || cstate.box_x.z > small_value
	    || cstate.box_y.z > small_value || cstate.box_y.x > small_value
	    || cstate.box_z.x > small_value || cstate.box_z.y > small_value){
        error(1,"KIM: simulation box of configuration %d is not orthogonal. Try to use 'NEIGH_RVEC' "
	      "instead of 'MI_OPBC'.\n", g_config.nconf);
      } else {
        /* store the box size info in box_side_len */
        g_kim.box_side_len[3*g_config.nconf + 0] = cstate.box_x.x;
        g_kim.box_side_len[3*g_config.nconf + 1] = cstate.box_y.y;
        g_kim.box_side_len[3*g_config.nconf + 2] = cstate.box_z.z;
      }
    }
#endif   // KIM

    // read the atoms
    for (int i = 0; i < cstate.atom_count; i++) {
      atom_t* atom = g_config.atoms + g_config.natoms + i;

      if (7 > fscanf(config_file, "%d %lf %lf %lf %lf %lf %lf\n", &(atom->type),
                     &(atom->pos.x), &(atom->pos.y), &(atom->pos.z),
                     &(atom->force.x), &(atom->force.y), &(atom->force.z)))
        error(1, "Corrupt configuration file on line %d\n", cstate.line + 1);

      cstate.line++;

      if (g_param.global_cell_scale != 1.0) {
        atom->pos.x *= g_param.global_cell_scale;
        atom->pos.y *= g_param.global_cell_scale;
        atom->pos.z *= g_param.global_cell_scale;
      }

      if (atom->type >= g_param.ntypes || atom->type < 0)
        error(
            1,
            "Corrupt configuration file on line %d: Incorrect atom type (%d)\n",
            cstate.line, atom->type);

      atom->absforce = sqrt(dsquare(atom->force.x) + dsquare(atom->force.y) +
                            dsquare(atom->force.z));

      atom->conf = g_config.nconf;

#if defined(CONTRIB)
      if (cstate.have_contrib_box_vector || cstate.n_spheres != 0)
        atom->contrib = does_contribute(atom->pos, &cstate);
      else
        atom->contrib = 1;
#endif  // CONTRIB

      g_config.na_type[g_config.nconf][atom->type] += 1;

      max_atom_type = MAX(max_atom_type, atom->type);
    }

    init_box_vectors(&cstate);

    init_neighbors(&cstate, mindist);

    init_angles(&cstate);

    // increment natoms and configuration number
    g_config.natoms += cstate.atom_count;
    g_config.nconf++;

  } while (!feof(config_file));

  // close config file
  fclose(config_file);

  // the calculation of the neighbor lists is now complete
  printf(
      "Reading the config file >> %s << and calculating neighbor lists ... "
      "done\n",
      filename);

  // calculate the total number of the atom types
  g_config.na_type =
      (int**)Realloc(g_config.na_type, (g_config.nconf + 1) * sizeof(int*));

  g_config.na_type[g_config.nconf] = (int*)Malloc(g_param.ntypes * sizeof(int));

  for (int i = 0; i < g_config.nconf; i++)
    for (int j = 0; j < g_param.ntypes; j++)
      g_config.na_type[g_config.nconf][j] += g_config.na_type[i][j];

  // print diagnostic message
  printf("\nRead %d configurations (%d with forces, %d with stresses)\n",
         g_config.nconf, with_forces, with_stresses);
  printf("with a total of %d atoms (", g_config.natoms);

  for (int i = 0; i < g_param.ntypes; i++) {
    printf("%d %s (%.2f%%)", g_config.na_type[g_config.nconf][i],
           g_config.elements[i],
           100.0 * g_config.na_type[g_config.nconf][i] / g_config.natoms);
    if (i != (g_param.ntypes - 1))
      printf(", ");
  }
  printf(").\n");

  // be pedantic about too large g_param.ntypes
  if ((max_atom_type + 1) < g_param.ntypes) {
    error(0, "There are less than %d atom types in your configurations!\n",
          g_param.ntypes);
    error(1, "Please adjust \"ntypes\" in your parameter file.\n");
  }

  /* mdim is the dimension of the force vector:
     - 3*natoms forces
     - nconf cohesive energies,
     - 6*nconf stress tensor components */
  g_calc.mdim = 3 * g_config.natoms + g_config.nconf;
#if defined(STRESS)
  g_calc.mdim += 6 * g_config.nconf;
#endif  // STRESS

/* mdim has additional components for EAM-like potentials */
#if defined(EAM) || defined(ADP) || defined(MEAM)
  g_calc.mdim += g_config.nconf;     /* nconf limiting constraints */
  g_calc.mdim += 2 * g_param.ntypes; /* g_param.ntypes dummy constraints */
#if defined(TBEAM)
  g_calc.mdim +=
      2 * g_param.ntypes; /* additional dummy constraints for s-band */
#endif                    // TBEAM
#endif                    // EAM || ADP || MEAM

/* mdim has additional components for analytic potentials */
#if defined(APOT)
  /* 1 slot for each analytic parameter -> punishment */
  g_calc.mdim += g_pot.opt_pot.idxlen;
  /* 1 slot for each analytic potential -> punishment */
  g_calc.mdim += g_pot.apot_table.number + 1;
#endif  // APOT

  /* copy forces into single vector */
  g_config.force_0 = (double*)Malloc(g_calc.mdim * sizeof(double));

  int k = 0;

  /* first forces */
  for (int i = 0; i < g_config.natoms; i++) {
    g_config.force_0[k++] = g_config.atoms[i].force.x;
    g_config.force_0[k++] = g_config.atoms[i].force.y;
    g_config.force_0[k++] = g_config.atoms[i].force.z;
  }

  // then cohesive energies
  for (int i = 0; i < g_config.nconf; i++)
    g_config.force_0[k++] = g_config.coheng[i];

#if defined(STRESS)
  // then stresses
  for (int i = 0; i < g_config.nconf; i++) {
    if (g_config.usestress[i]) {
      g_config.force_0[k++] = g_config.stress[i].xx;
      g_config.force_0[k++] = g_config.stress[i].yy;
      g_config.force_0[k++] = g_config.stress[i].zz;
      g_config.force_0[k++] = g_config.stress[i].xy;
      g_config.force_0[k++] = g_config.stress[i].yz;
      g_config.force_0[k++] = g_config.stress[i].zx;
    } else
      k += 6;
  }
#endif  // STRESS

  if (g_param.write_pair == 1)
    write_pair_distribution_file();

/* assign correct distances to different tables */
#if defined(APOT)
  double min = 10.0;

  /* pair potentials */
  for (int i = 0; i < g_param.ntypes; i++) {
    for (int j = 0; j < g_param.ntypes; j++) {
      k = (i <= j) ? i * g_param.ntypes + j - ((i * (i + 1)) / 2)
                   : j * g_param.ntypes + i - ((j * (j + 1)) / 2);
      if (mindist[k] >= 99.9)
        mindist[k] = 2.5;
      g_config.rmin[i * g_param.ntypes + j] = mindist[k];
      g_pot.apot_table.begin[k] = mindist[k] * 0.95;
      g_pot.opt_pot.begin[k] = mindist[k] * 0.95;
      g_pot.calc_pot.begin[k] = mindist[k] * 0.95;
      min = MIN(min, mindist[k]);
    }
  }

/* transfer functions */
#if defined(EAM) || defined(ADP) || defined(MEAM)
  for (int i = g_calc.paircol; i < g_calc.paircol + g_param.ntypes; i++) {
    g_pot.apot_table.begin[i] = min * 0.95;
    g_pot.opt_pot.begin[i] = min * 0.95;
    g_pot.calc_pot.begin[i] = min * 0.95;
  }
#if defined(TBEAM)
  for (int i = g_calc.paircol + 2 * g_param.ntypes;
       i < g_calc.paircol + 3 * g_param.ntypes; i++) {
    g_pot.apot_table.begin[i] = min * 0.95;
    g_pot.opt_pot.begin[i] = min * 0.95;
    g_pot.calc_pot.begin[i] = min * 0.95;
  }
#endif  // TBEAM
#endif  // EAM || ADP || MEAM

/* dipole and quadrupole functions */
#if defined(ADP)
  for (int i = 0; i < g_calc.paircol; i++) {
    int j = g_calc.paircol + 2 * g_param.ntypes + i;
    g_pot.apot_table.begin[j] = min * 0.95;
    g_pot.opt_pot.begin[j] = min * 0.95;
    g_pot.calc_pot.begin[j] = min * 0.95;
    j = 2 * g_calc.paircol + 2 * g_param.ntypes + i;
    g_pot.apot_table.begin[j] = min * 0.95;
    g_pot.opt_pot.begin[j] = min * 0.95;
    g_pot.calc_pot.begin[j] = min * 0.95;
  }
#endif  // ADP

#if defined(MEAM)
  /* f_ij */
  for (int i = 0; i < g_calc.paircol; i++) {
    int j = g_calc.paircol + 2 * g_param.ntypes + i;
    g_pot.apot_table.begin[j] = min * 0.95;
    g_pot.opt_pot.begin[j] = min * 0.95;
    g_pot.calc_pot.begin[j] = min * 0.95;
  }
  /* g_i */
  /* g_i takes cos(theta) as an argument, so we need to tabulate it only
     in the range of [-1:1]. Actually we use [-1.1:1.1] to be safe. */
  for (int i = 0; i < g_param.ntypes; i++) {
    int j = 2 * g_calc.paircol + 2 * g_param.ntypes + i;
    g_pot.apot_table.begin[j] = -1.1;
    g_pot.opt_pot.begin[j] = -1.1;
    g_pot.calc_pot.begin[j] = -1.1;
    g_pot.apot_table.end[j] = 1.1;
    g_pot.opt_pot.end[j] = 1.1;
    g_pot.calc_pot.end[j] = 1.1;
  }
#endif  // MEAM

  /* recalculate step, invstep and xcoord for new tables */
  for (int i = 0; i < g_pot.calc_pot.ncols; i++) {
    g_pot.calc_pot.step[i] =
        (g_pot.calc_pot.end[i] - g_pot.calc_pot.begin[i]) / (APOT_STEPS - 1);
    g_pot.calc_pot.invstep[i] = 1.0 / g_pot.calc_pot.step[i];
    for (int j = 0; j < APOT_STEPS; j++) {
      int index = i * APOT_STEPS + (i + 1) * 2 + j;
      g_pot.calc_pot.xcoord[index] =
          g_pot.calc_pot.begin[i] + j * g_pot.calc_pot.step[i];
    }
  }

  update_slots();
#endif  // APOT

  print_minimal_distances_matrix(mindist);
}

#if defined(APOT)

/****************************************************************
  update_slots:
    recalculate the slots of the atoms for analytic potential
****************************************************************/

void update_slots(void)
{
  for (int i = 0; i < g_config.natoms; i++) {
    for (int j = 0; j < g_config.atoms[i].num_neigh; j++) {
      double r = g_config.atoms[i].neigh[j].r;

      // update slots for pair potential part, slot 0
      update_neighbor_slots(g_config.atoms[i].neigh + j, r, 0);

#if defined EAM || defined ADP || defined MEAM
      // update slots for eam transfer functions, slot 1
      update_neighbor_slots(g_config.atoms[i].neigh + j, r, 1);
#if defined(TBEAM)
      // update slots for tbeam transfer functions, s-band, slot 2
      update_neighbor_slots(g_config.atoms[i].neigh + j, r, 2);
#endif  // TBEAM
#endif  // EAM || ADP || MEAM

#if defined(MEAM)
      // update slots for MEAM f functions, slot 2
      update_neighbor_slots(g_config.atoms[i].neigh + j, r, 2);
#endif  // MEAM

#if defined(ADP)
      // update slots for adp dipole functions, slot 2
      update_neighbor_slots(g_config.atoms[i].neigh + j, r, 2);
      // update slots for adp quadrupole functions, slot 3
      update_neighbor_slots(g_config.atoms[i].neigh + j, r, 3);
#endif  // ADP

    }  // end loop over all neighbors
  }    // end loop over all atoms

#if defined(THREEBODY) && defined(MEAM)
  // update angular slots
  for (int i = 0; i < g_config.natoms; i++) {
    for (int j = 0; j < g_config.atoms[i].num_angles; j++) {
      double rr = g_config.atoms[i].angle_part[j].cos + 1.1;
      int col =
          2 * g_calc.paircol + 2 * g_param.ntypes + g_config.atoms[i].type;
      g_config.atoms[i].angle_part[j].slot =
          (int)(rr * g_pot.calc_pot.invstep[col]);
      g_config.atoms[i].angle_part[j].step = g_pot.calc_pot.step[col];
      g_config.atoms[i].angle_part[j].shift =
          (rr -
           g_config.atoms[i].angle_part[j].slot * g_pot.calc_pot.step[col]) *
          g_pot.calc_pot.invstep[col];
      // move slot to the correct potential
      g_config.atoms[i].angle_part[j].slot += g_pot.calc_pot.first[col];
    }
  }
#endif  // THREEBODY && MEAM

#if defined(STIWEB)
  g_pot.apot_table.sw.init = 0;
#endif  // STIWEB

#if defined(TERSOFF)
  g_pot.apot_table.tersoff.init = 0;
#endif  // TERSOFF
}

/****************************************************************
  update_neighbor_slots
    recalculate the slots of the neighbors for analytic potential
****************************************************************/

void update_neighbor_slots(neigh_t* neighbor, double r, int neighbor_slot)
{
  int col = neighbor->col[neighbor_slot];

  if (r < g_pot.calc_pot.end[col]) {
    double rr = r - g_pot.calc_pot.begin[col];
    neighbor->slot[neighbor_slot] = (int)(rr * g_pot.calc_pot.invstep[col]);
    neighbor->step[neighbor_slot] = g_pot.calc_pot.step[col];
    neighbor->shift[neighbor_slot] =
        (rr - neighbor->slot[neighbor_slot] * g_pot.calc_pot.step[col]) *
        g_pot.calc_pot.invstep[col];
    // move slot to the correct potential
    neighbor->slot[neighbor_slot] += g_pot.calc_pot.first[col];
  }
}

#endif  // APOT

/****************************************************************
  reset_cstate
****************************************************************/

void reset_cstate(config_state* cstate)
{
  cstate->atom_count = 0;
  cstate->use_force = 0;
  cstate->have_energy = 0;
  cstate->have_stress = 0;
  cstate->have_stress = 0;
  cstate->box_x.x = 0.0;
  cstate->box_x.y = 0.0;
  cstate->box_x.z = 0.0;
  cstate->box_y.x = 0.0;
  cstate->box_y.y = 0.0;
  cstate->box_y.z = 0.0;
  cstate->box_z.x = 0.0;
  cstate->box_z.y = 0.0;
  cstate->box_z.z = 0.0;
  cstate->tbox_x.x = 0.0;
  cstate->tbox_x.y = 0.0;
  cstate->tbox_x.z = 0.0;
  cstate->tbox_y.x = 0.0;
  cstate->tbox_y.y = 0.0;
  cstate->tbox_y.z = 0.0;
  cstate->tbox_z.x = 0.0;
  cstate->tbox_z.y = 0.0;
  cstate->tbox_z.z = 0.0;
  cstate->cell_scale.x = 0.0;
  cstate->cell_scale.y = 0.0;
  cstate->cell_scale.z = 0.0;
  cstate->have_box_vector = 0;
  cstate->stresses = NULL;
#if defined(CONTRIB)
  cstate->have_contrib_box_vector = 0;
  cstate->n_spheres = 0;
  cstate->cbox_o.x = 0.0;
  cstate->cbox_o.y = 0.0;
  cstate->cbox_o.z = 0.0;
  cstate->cbox_a.x = 0.0;
  cstate->cbox_a.y = 0.0;
  cstate->cbox_a.z = 0.0;
  cstate->cbox_b.x = 0.0;
  cstate->cbox_b.y = 0.0;
  cstate->cbox_b.z = 0.0;
  cstate->cbox_c.x = 0.0;
  cstate->cbox_c.y = 0.0;
  cstate->cbox_c.z = 0.0;
  cstate->sphere_center = NULL;
  cstate->sphere_radius = NULL;
#endif  // CONTRIB
}

/****************************************************************
  create_memory_for_config
****************************************************************/

void create_memory_for_config(config_state* cstate)
{
  int nconf = g_config.nconf + 1;

  // increase memory for this many additional atoms
  g_config.atoms = (atom_t*)Realloc(
      g_config.atoms, (g_config.natoms + cstate->atom_count) * sizeof(atom_t));

  for (int i = 0; i < cstate->atom_count; i++)
    g_config.atoms[g_config.natoms + i].neigh =
        (neigh_t*)Malloc(sizeof(neigh_t));

  g_config.coheng = (double*)Realloc(g_config.coheng, nconf * sizeof(double));
  g_config.coheng[g_config.nconf] = 0.0;

  g_config.conf_weight =
      (double*)Realloc(g_config.conf_weight, nconf * sizeof(double));
  g_config.conf_weight[g_config.nconf] = 1.0;

  g_config.volume = (double*)Realloc(g_config.volume, nconf * sizeof(double));
  g_config.volume[g_config.nconf] = 0.0;

#if defined(STRESS)
  g_config.stress =
      (sym_tens*)Realloc(g_config.stress, nconf * sizeof(sym_tens));
  memset(&g_config.stress[g_config.nconf], 0,
         sizeof(g_config.stress[g_config.nconf]));
  g_config.usestress = (int*)Realloc(g_config.usestress, nconf * sizeof(int));
  g_config.usestress[g_config.nconf] = 0;
#endif  // STRESS

  g_config.inconf = (int*)Realloc(g_config.inconf, nconf * sizeof(int));
  g_config.inconf[g_config.nconf] = 0;

  g_config.cnfstart = (int*)Realloc(g_config.cnfstart, nconf * sizeof(int));
  g_config.cnfstart[g_config.nconf] = 0;

  g_config.useforce = (int*)Realloc(g_config.useforce, nconf * sizeof(int));
  g_config.useforce[g_config.nconf] = 0;

  g_config.na_type =
      (int**)Realloc(g_config.na_type, (nconf + 1) * sizeof(int*));
  g_config.na_type[g_config.nconf] = (int*)Malloc(g_param.ntypes * sizeof(int));

  for (int i = 0; i < g_param.ntypes; i++)
    g_config.na_type[g_config.nconf][i] = 0;
}

/****************************************************************
  read_box_vector
****************************************************************/

void read_box_vector(char const* pline, vector* pvect, config_state* cstate)
{
  if (sscanf(pline, "%lf %lf %lf\n", &pvect->x, &pvect->y, &pvect->z) == 3) {
    if (g_param.global_cell_scale != 1.0) {
      pvect->x *= g_param.global_cell_scale;
      pvect->y *= g_param.global_cell_scale;
      pvect->z *= g_param.global_cell_scale;
    }
  } else
    error(1, "%s:%d Error reading box vector\n", cstate->filename,
          cstate->line);
}

#if defined(CONTRIB)

/****************************************************************
  read_sphere_center
****************************************************************/

void read_sphere_center(char const* pline, config_state* cstate)
{
  int n = cstate->n_spheres;

  cstate->sphere_center =
      (vector*)Realloc(cstate->sphere_center, (n + 1) * sizeof(vector));
  cstate->sphere_radius =
      (double*)Realloc(cstate->sphere_radius, (n + 1) * sizeof(double));

  if (sscanf(pline, "%lf %lf %lf %lf\n", &cstate->sphere_center[n].x,
             &cstate->sphere_center[n].y, &cstate->sphere_center[n].z,
             &cstate->sphere_radius[n]) == 4) {
    if (g_param.global_cell_scale != 1.0) {
      cstate->sphere_center[n].x *= g_param.global_cell_scale;
      cstate->sphere_center[n].y *= g_param.global_cell_scale;
      cstate->sphere_center[n].z *= g_param.global_cell_scale;
      cstate->sphere_radius[n] *= g_param.global_cell_scale;
    }
  } else
    error(1, "%s:%d Error reading sphere of contributing atoms\n",
          cstate->filename, cstate->line);
}

#endif  // CONTRIB

/****************************************************************
  read_chemical_elements
****************************************************************/

void read_chemical_elements(char* psrc, config_state* cstate)
{
  int i = 0;
  char const* pchar = strtok(psrc + 3, " \t");

  if (!cstate->num_fixed_elements) {
    while (pchar != NULL && i < g_param.ntypes) {
      int len = max(strlen(pchar), 4);
      strncpy((char*)g_config.elements[i], pchar, len);
      *((char*)g_config.elements[i] + len) = '\0';
      pchar = strtok(NULL, " \t");
      i++;
      cstate->num_fixed_elements++;
    }
  } else {
    while (pchar != NULL && i < g_param.ntypes) {
      if (strcmp(pchar, g_config.elements[i]) != 0) {
        if (atoi(g_config.elements[i]) == i && i > cstate->num_fixed_elements) {
          int len = max(strlen(pchar), 4);
          strncpy((char*)g_config.elements[i], pchar, len);
          *((char*)g_config.elements[i] + len) = '\0';
          cstate->num_fixed_elements++;
        } else {
          error(0, "Mismatch found in configuration %d, line %d.\n",
                cstate->config, cstate->line);
          error(0, "Expected element >> %s << but found element >> %s <<.\n",
                g_config.elements[i], pchar);
          error(0, "You can use list_config to identify that configuration.\n");
          error(1, "Please check your configuration files!\n");
        }
      }
      pchar = strtok(NULL, " \t");
      i++;
    }
  }
}

/****************************************************************
  init_box_vectors
****************************************************************/

void init_box_vectors(config_state* cstate)
{
  // check cell size

  // inverse height in each direction
  vector iheight;

  iheight.x = sqrt(SPROD(cstate->tbox_x, cstate->tbox_x));
  iheight.y = sqrt(SPROD(cstate->tbox_y, cstate->tbox_y));
  iheight.z = sqrt(SPROD(cstate->tbox_z, cstate->tbox_z));

  if ((ceil(g_config.rcutmax * iheight.x) > 30000) ||
      (ceil(g_config.rcutmax * iheight.y) > 30000) ||
      (ceil(g_config.rcutmax * iheight.z) > 30000))
    error(1, "Very bizarre small cell size - aborting\n");

  cstate->cell_scale.x = (int)ceil(g_config.rcutmax * iheight.x);
  cstate->cell_scale.y = (int)ceil(g_config.rcutmax * iheight.y);
  cstate->cell_scale.z = (int)ceil(g_config.rcutmax * iheight.z);

  if (cstate->cell_scale.x > 1 || cstate->cell_scale.y > 1 ||
      cstate->cell_scale.z > 1) {
    warning(
        "The box size of configuration %d is smaller than the cutoff "
        "distance.\n",
        cstate->config);
    warning(
        "  Using additional periodic images for energy and force "
        "calculations.\n");
  }

#if defined(DEBUG)
  fprintf(stderr, "\nChecking cell size for configuration %d:\n",
          g_config.nconf);
  fprintf(stderr, "Box dimensions:\n");
  fprintf(stderr, "     %10.6f %10.6f %10.6f\n", cstate->box_x.x,
          cstate->box_x.y, cstate->box_x.z);
  fprintf(stderr, "     %10.6f %10.6f %10.6f\n", cstate->box_y.x,
          cstate->box_y.y, cstate->box_y.z);
  fprintf(stderr, "     %10.6f %10.6f %10.6f\n", cstate->box_z.x,
          cstate->box_z.y, cstate->box_z.z);
  fprintf(stderr, "Box normals:\n");
  fprintf(stderr, "     %10.6f %10.6f %10.6f\n", cstate->tbox_x.x,
          cstate->tbox_x.y, cstate->tbox_x.z);
  fprintf(stderr, "     %10.6f %10.6f %10.6f\n", cstate->tbox_y.x,
          cstate->tbox_y.y, cstate->tbox_y.z);
  fprintf(stderr, "     %10.6f %10.6f %10.6f\n", cstate->tbox_z.x,
          cstate->tbox_z.y, cstate->tbox_z.z);
  fprintf(stderr, "Box heights:\n");
  fprintf(stderr, "     %10.6f %10.6f %10.6f\n", 1.0 / iheight.x,
          1.0 / iheight.y, 1.0 / iheight.z);
  fprintf(stderr, "Potential range:  %f\n", g_config.rcutmax);
  fprintf(stderr, "Periodic images needed: %d %d %d\n\n",
          (int)(2 * cstate->cell_scale.x + 1),
          (int)(2 * cstate->cell_scale.y + 1),
          (int)(2 * cstate->cell_scale.z + 1));
#endif  // DEBUG
}

/****************************************************************
  init_neighbors
****************************************************************/

void init_neighbors(config_state* cstate, double* mindist)
{
  vector d;
  vector dd;

  // compute the neighbor table
  for (int i = g_config.natoms; i < g_config.natoms + cstate->atom_count; i++) {
/* loop over all atoms for threebody interactions */
#if defined(THREEBODY)
    for (int j = g_config.natoms; j < g_config.natoms + cstate->atom_count; j++)
#else
    int j_start = i;
#if defined(KIM)
    if (g_kim.is_half_neighbors != 1)
      j_start = g_config.natoms;
#endif
    for (int j = j_start; j < g_config.natoms + cstate->atom_count; j++)
#endif  // THREEBODY
    {
      d.x = g_config.atoms[j].pos.x - g_config.atoms[i].pos.x;
      d.y = g_config.atoms[j].pos.y - g_config.atoms[i].pos.y;
      d.z = g_config.atoms[j].pos.z - g_config.atoms[i].pos.z;

      for (int ix = -cstate->cell_scale.x; ix <= cstate->cell_scale.x; ix++) {
        for (int iy = -cstate->cell_scale.y; iy <= cstate->cell_scale.y; iy++) {
          for (int iz = -cstate->cell_scale.z; iz <= cstate->cell_scale.z;
               iz++) {
            if ((i == j) && (ix == 0) && (iy == 0) && (iz == 0))
              continue;
            dd.x = d.x + ix * cstate->box_x.x + iy * cstate->box_y.x +
                   iz * cstate->box_z.x;
            dd.y = d.y + ix * cstate->box_x.y + iy * cstate->box_y.y +
                   iz * cstate->box_z.y;
            dd.z = d.z + ix * cstate->box_x.z + iy * cstate->box_y.z +
                   iz * cstate->box_z.z;
            double r = sqrt(SPROD(dd, dd));
            int type1 = g_config.atoms[i].type;
            int type2 = g_config.atoms[j].type;

            if (r <= g_config.rcut[type1 * g_param.ntypes + type2]) {
              int short_distance = 0;

              if (r <= g_config.rmin[type1 * g_param.ntypes + type2]) {
                warning("Configuration %i: Distance %i\n", cstate->config, r);
                warning(" atom %d (type %d) at pos: %f %f %f\n",
                        i - g_config.natoms, type1, g_config.atoms[i].pos.x,
                        g_config.atoms[i].pos.y, g_config.atoms[i].pos.z);
                warning(" atom %d (type %d) at pos: %f %f %f\n",
                        j - g_config.natoms, type2, dd.x, dd.y, dd.z);
#if !defined(KIM)
                short_distance = 1;
#endif // !KIM
              }
              g_config.atoms[i].neigh = (neigh_t*)Realloc(
                  g_config.atoms[i].neigh,
                  (g_config.atoms[i].num_neigh + 1) * sizeof(neigh_t));
              dd.x /= r;
              dd.y /= r;
              dd.z /= r;

              int k = g_config.atoms[i].num_neigh++;

              neigh_t* n = g_config.atoms[i].neigh + k;

              memset(n, 0, sizeof(neigh_t));

              n->type = type2;
              n->nr = j;
              n->r = r;
              n->r2 = r * r;
              n->inv_r = 1.0 / r;
              n->dist_r = dd;
              n->dist.x = dd.x * r;
              n->dist.y = dd.y * r;
              n->dist.z = dd.z * r;

#if defined(ADP)
              n->sqrdist.xx = dd.x * dd.x * r * r;
              n->sqrdist.yy = dd.y * dd.y * r * r;
              n->sqrdist.zz = dd.z * dd.z * r * r;
              n->sqrdist.yz = dd.y * dd.z * r * r;
              n->sqrdist.zx = dd.z * dd.x * r * r;
              n->sqrdist.xy = dd.x * dd.y * r * r;
#endif  // ADP

              /* pre-compute index and shift into potential table */

              if (!short_distance) {
                /* pair potential */
                int col = (type1 <= type2)
                              ? type1 * g_param.ntypes + type2 -
                                    ((type1 * (type1 + 1)) / 2)
                              : type2 * g_param.ntypes + type1 -
                                    ((type2 * (type2 + 1)) / 2);
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 0);

                mindist[col] = MIN(mindist[col], r);

#if defined(EAM) || defined(ADP) || defined(MEAM)
                /* transfer function */
                col = g_calc.paircol + type2;
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 1);
#if defined(TBEAM)
                /* transfer function - d band */
                col = g_calc.paircol + 2 * g_param.ntypes + type2;
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 2);
#endif  // TBEAM
#endif  // EAM || ADP || MEAM

#if defined(MEAM)
                /* Store slots and stuff for f(r_ij) */
                col = g_calc.paircol + 2 * g_param.ntypes +
                      g_config.atoms[i].neigh[k].col[0];
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 2);
#endif  // MEAM

#if defined(ADP)
                /* dipole part */
                col = g_calc.paircol + 2 * g_param.ntypes +
                      g_config.atoms[i].neigh[k].col[0];
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 2);

                /* quadrupole part */
                col = 2 * g_calc.paircol + 2 * g_param.ntypes +
                      g_config.atoms[i].neigh[k].col[0];
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 3);
#endif  // ADP

#if defined(STIWEB)
                /* Store slots and stuff for exp. function */
                col = g_calc.paircol + g_config.atoms[i].neigh[k].col[0];
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 1);
#endif  // STIWEB

              } /* !sh_dist */
            }   /* r < r_cut */
          }     /* loop over images in z direction */
        }       /* loop over images in y direction */
      }         /* loop over images in x direction */
    }           /* second loop over atoms (neighbors) */
  }             /* first loop over atoms */
}

/****************************************************************
  set_neighbor_slot
    compute box transformation matrix
****************************************************************/

void set_neighbor_slot(neigh_t* neighbor, int col, double r, int store_slot)
{
  int slot = 0;
  double step = 0.0;
  double istep = 0.0;
  double shift = 0.0;

  switch (g_pot.format_type) {
    case POTENTIAL_FORMAT_ANALYTIC:
    case POTENTIAL_FORMAT_TABULATED_EQ_DIST: {
      double rr = r - g_pot.calc_pot.begin[col];

      if (rr < 0) {
        // TODO: rephrase
        error(0, "The distance %f is smaller than the beginning\n", r);
        error(0, "of the potential #%d (r_begin=%f).\n", col,
              g_pot.calc_pot.begin[col]);
        fflush(stdout);
        error(1, "Short distance!\n");
      }

      istep = g_pot.calc_pot.invstep[col];
      slot = (int)(rr * istep);
      shift = (rr - slot * g_pot.calc_pot.step[col]) * istep;
      slot += g_pot.calc_pot.first[col];
      step = g_pot.calc_pot.step[col];
      break;
    }
    case POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST: {
      int klo = g_pot.calc_pot.first[col];
      int khi = g_pot.calc_pot.last[col];

      /* bisection */
      while (khi - klo > 1) {
        slot = (khi + klo) >> 1;
        if (g_pot.calc_pot.xcoord[slot] > r)
          khi = slot;
        else
          klo = slot;
      }

      slot = klo;
      step = g_pot.calc_pot.xcoord[khi] - g_pot.calc_pot.xcoord[klo];
      shift = (r - g_pot.calc_pot.xcoord[klo]) / step;
      break;
    }
    case POTENTIAL_FORMAT_KIM:
      return;
    case POTENTIAL_FORMAT_UNKNOWN:
      error(1, "Unknown potential format detected.\n");
  }

  // independent of format - we should be left of last index
  if (slot >= g_pot.calc_pot.last[col]) {
    slot--;
    shift += 1.0;
  }

  neighbor->shift[store_slot] = shift;
  neighbor->slot[store_slot] = slot;
  neighbor->step[store_slot] = step;
  neighbor->col[store_slot] = col;
}

/****************************************************************
  init_angles
****************************************************************/

void init_angles(config_state* cstate)
{
// compute the angular part
// for TERSOFF we create a full neighbor list,
// for all other potentials only a half list
#if defined(THREEBODY)
  for (int i = g_config.natoms; i < g_config.natoms + cstate->atom_count; i++) {
    int nnn = g_config.atoms[i].num_neigh;
    int ijk = 0;
    g_config.atoms[i].angle_part = (angle_t*)Malloc(sizeof(angle_t));
#if defined(TERSOFF)
    for (int j = 0; j < nnn; j++) {
#else
    for (int j = 0; j < nnn - 1; j++) {
#endif  // TERSOFF
      g_config.atoms[i].neigh[j].ijk_start = ijk;
#if defined(TERSOFF)
      for (int k = 0; k < nnn; k++) {
        if (j == k)
          continue;
#else
      for (int k = j + 1; k < nnn; k++) {
#endif  // TERSOFF

        g_config.atoms[i].angle_part = (angle_t*)Realloc(
            g_config.atoms[i].angle_part, (ijk + 1) * sizeof(angle_t));

        memset(g_config.atoms[i].angle_part + ijk, 0, sizeof(angle_t));

        double ccos = g_config.atoms[i].neigh[j].dist_r.x *
                          g_config.atoms[i].neigh[k].dist_r.x +
                      g_config.atoms[i].neigh[j].dist_r.y *
                          g_config.atoms[i].neigh[k].dist_r.y +
                      g_config.atoms[i].neigh[j].dist_r.z *
                          g_config.atoms[i].neigh[k].dist_r.z;

        g_config.atoms[i].angle_part[ijk].cos = ccos;

        int col =
            2 * g_calc.paircol + 2 * g_param.ntypes + g_config.atoms[i].type;

        if (g_pot.format_type == POTENTIAL_FORMAT_ANALYTIC ||
            g_pot.format_type == POTENTIAL_FORMAT_TABULATED_EQ_DIST) {
          if ((fabs(ccos) - 1.0) > 1e-10) {
            int type1 = g_config.atoms[i].type;
            int type2 = g_config.atoms[i].neigh[j].type;
            printf("%.20f %f %d %d %d\n", ccos, g_pot.calc_pot.begin[col], col,
                   type1, type2);
            fflush(stdout);
            error(1, "cos out of range, it is strange!\n");
          }
#if defined(MEAM)
          double istep = g_pot.calc_pot.invstep[col];
          int slot = (int)((ccos + 1) * istep);
          double shift = ((ccos + 1) - slot * g_pot.calc_pot.step[col]) * istep;
          slot += g_pot.calc_pot.first[col];
          //             double step = g_pot.calc_pot.step[col];

          /* Don't want lower bound spline knot to be final knot or upper
             bound knot will cause trouble since it goes beyond the array */
          if (slot >= g_pot.calc_pot.last[col]) {
            slot--;
            shift += 1.0;
          }
#endif  // MEAM
        }
#if defined(MEAM)
// TODO: how did this ever work ???
//         g_config.atoms[i].angle_part[ijk].shift = shift;
//         g_config.atoms[i].angle_part[ijk].slot = slot;
//         g_config.atoms[i].angle_part[ijk].step = step;
#endif  // MEAM
        ijk++;
      } /* third loop over atoms */
    }   /* second loop over atoms */
    g_config.atoms[i].num_angles = ijk;
  }     /* first loop over atoms */
#endif  // THREEBODY
}

/****************************************************************
  make_box
****************************************************************/

double make_box(config_state* cstate)
{
  double volume = 0;

  // compute tbox_j such that SPROD(box_i,tbox_j) == delta_ij

  // first unnormalized
  cstate->tbox_x = vec_prod(cstate->box_y, cstate->box_z);
  cstate->tbox_y = vec_prod(cstate->box_z, cstate->box_x);
  cstate->tbox_z = vec_prod(cstate->box_x, cstate->box_y);

  // volume
  volume = SPROD(cstate->box_x, cstate->tbox_x);

  if (0.0 == volume)
    error(1, "Box edges are parallel\n");

  // normalization
  cstate->tbox_x.x /= volume;
  cstate->tbox_x.y /= volume;
  cstate->tbox_x.z /= volume;
  cstate->tbox_y.x /= volume;
  cstate->tbox_y.y /= volume;
  cstate->tbox_y.z /= volume;
  cstate->tbox_z.x /= volume;
  cstate->tbox_z.y /= volume;
  cstate->tbox_z.z /= volume;

  return volume;
}

/****************************************************************
  write_pair_distribution_file
****************************************************************/

void write_pair_distribution_file()
{
  char pairname[255];
#if defined(APOT)
  int pair_steps = APOT_STEPS / 2;
#else
  int pair_steps = 500;
#endif  // APOT
  if (g_calc.paircol == 0) {
    warning("write_pair_distribution_file failed, paircol = 0\n");
    return;
  }
  double* pair_table =
      (double*)Malloc(g_calc.paircol * pair_steps * sizeof(double));
  double* pair_dist = (double*)Malloc(g_calc.paircol * sizeof(double));
  int pos = 0;
  int max_count = 0;

  sprintf(pairname, "%s.pair", g_files.config);
  FILE* pairfile = fopen(pairname, "w");

  if (pairfile == NULL) {
    warning("Error opening pair distribution file %s: %s\n", pairname,
            strerror(errno));
    return;
  }

  fprintf(pairfile, "# radial distribution file for %d potential(s)\n",
          g_calc.paircol);

  for (int i = 0; i < g_param.ntypes; i++) {
    for (int k = 0; k < g_param.ntypes; k++) {
      int index = (i <= k) ? i * g_param.ntypes + k - (i * (i + 1) / 2)
                           : k * g_param.ntypes + i - (k * (k + 1) / 2);
      pair_dist[index] = g_config.rcut[i * g_param.ntypes + k] / pair_steps;
    }
  }

  for (int k = 0; k < g_calc.paircol; k++) {
    for (int i = 0; i < g_config.natoms; i++) {
      for (int j = 0; j < g_config.atoms[i].num_neigh; j++) {
        int col = g_config.atoms[i].neigh[j].col[0];

        if (col == k) {
          pos = (int)(g_config.atoms[i].neigh[j].r / pair_dist[k]);
#if defined(DEBUG)
          if (g_config.atoms[i].neigh[j].r <= 1) {
            warning("Short distance (%f) found.\n",
                    g_config.atoms[i].neigh[j].r);
            warning("\tatom=%d neighbor=%d\n", i, j);
          }
#endif  // DEBUG
          pair_table[k * pair_steps + pos]++;
          if ((int)pair_table[k * pair_steps + pos] > max_count)
            max_count = (int)pair_table[k * pair_steps + pos];
        }
      }
    }
  }

  for (int k = 0; k < g_calc.paircol; k++) {
    for (int i = 0; i < pair_steps; i++) {
      pair_table[k * pair_steps + i] /= max_count;
      fprintf(pairfile, "%f %f\n", i * pair_dist[k],
              pair_table[k * pair_steps + i]);
    }
    if (k != (g_calc.paircol - 1))
      fprintf(pairfile, "\n\n");
  }
  fclose(pairfile);
}

/****************************************************************
  print_minimal_distances_matrix
****************************************************************/

void print_minimal_distances_matrix(const double* mindist)
{
  int i = 0;
  int j = 0;
  int k = 0;

  printf("\nMinimal Distances Matrix:\n");
  printf("Atom\t");
  for (i = 0; i < g_param.ntypes; i++)
    printf("%8s\t", g_config.elements[i]);
  printf("with\n");
  for (i = 0; i < g_param.ntypes; i++) {
    printf("%s\t", g_config.elements[i]);
    for (j = 0; j < g_param.ntypes; j++) {
      k = (i <= j) ? i * g_param.ntypes + j - ((i * (i + 1)) / 2)
                   : j * g_param.ntypes + i - ((j * (j + 1)) / 2);
      printf("%f\t", mindist[k]);
    }
    printf("\n");
  }
  printf("\n");
}

#if defined(CONTRIB)

/****************************************************************
  does_contribute
    check if the atom does contribute to the error sum
****************************************************************/

int does_contribute(vector pos, config_state* cstate)
{
  vector dist;

  if (cstate->have_contrib_box_vector) {
    dist.x = pos.x - cstate->cbox_o.x;
    dist.y = pos.y - cstate->cbox_o.y;
    dist.z = pos.z - cstate->cbox_o.z;
    int n_a =
        SPROD(dist, cstate->cbox_a) / SPROD(cstate->cbox_a, cstate->cbox_a);
    int n_b =
        SPROD(dist, cstate->cbox_b) / SPROD(cstate->cbox_b, cstate->cbox_b);
    int n_c =
        SPROD(dist, cstate->cbox_c) / SPROD(cstate->cbox_c, cstate->cbox_c);
    if (n_a >= 0 && n_a <= 1)
      if (n_b >= 0 && n_b <= 1)
        if (n_c >= 0 && n_c <= 1)
          return 1;
  }

  for (int i = 0; i < cstate->n_spheres; i++) {
    dist.x = (pos.x - cstate->sphere_center[i].x);
    dist.y = (pos.y - cstate->sphere_center[i].y);
    dist.z = (pos.z - cstate->sphere_center[i].z);
    double r = SPROD(dist, dist);
    r = sqrt(r);
    if (r < cstate->sphere_radius[i])
      return 1;
  }

  return 0;
}

#endif  // CONTRIB
