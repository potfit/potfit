/****************************************************************
 *
 * config.c: Reads atomic configurations and forces.
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
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

#include "potfit.h"

#include "config.h"
#include "utils.h"

/****************************************************************
 *
 *  read the configurations
 *
 ****************************************************************/

void read_config(char const* filename)
{
  config_state cstate;

  memset(&cstate, 0, sizeof(cstate));

  char buffer[1024];
  char* ptr = NULL;
  char* res = NULL;

  int i = 0;
  int j = 0;
  int k = 0;
  int line = 0;
  int max_atom_type = -1;
  int with_forces = 0;
  int with_stresses = 0;

  double* mindist = NULL;

  FILE* config_file = NULL;

  fpos_t filepos;

  /* initialize elements array */
  g_config.elements = (char **)malloc(g_param.ntypes * sizeof(char *));
  if (g_config.elements == NULL)
    error(1, "Cannot allocate memory for element names.");
  reg_for_free(g_config.elements, "elements");
  for (i = 0; i < g_param.ntypes; i++) {
    g_config.elements[i] = (char *)malloc(3 * sizeof(char));
    if (g_config.elements[i] == NULL)
      error(1, "Cannot allocate memory for %d. element name.\n", i + 1);
    reg_for_free(g_config.elements[i], "elements[%d]", i);
    snprintf(g_config.elements[i], 3, "%d", i);
  }

  /* initialize minimum distance array */
  mindist = (double *)malloc(g_param.ntypes * g_param.ntypes * sizeof(double));
  if (mindist == NULL)
    error(1, "Cannot allocate memory for minimal distance.");

  /* set maximum cutoff distance as starting value for mindist */
  for (i = 0; i < g_param.ntypes * g_param.ntypes; i++)
    mindist[i] = 99.9;
  for (i = 0; i < g_param.ntypes; i++) {
    for (j = 0; j < g_param.ntypes; j++) {
      k = (i <= j) ? i * g_param.ntypes + j - ((i * (i + 1)) / 2) : j * g_param.ntypes + i - ((j * (j + 1)) / 2);
      mindist[k] = MAX(g_config.rcut[i * g_param.ntypes + j], mindist[i * g_param.ntypes + j]);
    }
  }

  g_config.nconf = 0;

  /* open file */
  config_file = fopen(filename, "r");
  if (config_file == NULL)
    error(1, "Could not open file %s\n", filename);

  printf("Reading the config file >> %s << and calculating neighbor lists ... ", filename);
  fflush(stdout);

  /* read configurations until the end of the file */
  do
  {

    res = fgets(buffer, 1024, config_file);
    line++;

    if (NULL == res)
      error(1, "Unexpected end of file in %s", filename);

    if (res[0] == '#' && res[1] == 'N')
    {
      memset(&cstate, 0, sizeof(cstate));

      if (sscanf(res + 3, "%d %d", &cstate.atom_count, &cstate.use_force) < 2)
        error(1, "%s: Error in atom number specification on line %d\n", filename, line);
    }
    else
      continue;

    /* check if there are enough atoms, 2 for pair and 3 for manybody potentials */
#if !defined(THREEBODY)
    if (cstate.atom_count < 2)
      error(1, "The configuration %d (starting on line %d) has not enough atoms. Please remove it.", g_config.nconf + 1, line);
#else
    if (cstate.atom_count < 3)
      error(1, "The configuration %d (starting on line %d) has not enough atoms. Please remove it.", g_config.nconf + 1, line);
#endif /* !THREEBODY */

    create_memory_for_config(&cstate);

    for (i = g_config.natoms; i < g_config.natoms + cstate.atom_count; i++)
      init_atom_memory(g_config.atoms + i);

    for (i = 0; i < g_param.ntypes; i++)
      g_config.na_type[g_config.nconf][i] = 0;

    g_config.inconf[g_config.nconf] = cstate.atom_count;
    g_config.cnfstart[g_config.nconf] = g_config.natoms;
    g_config.useforce[g_config.nconf] = cstate.use_force;
#if defined(STRESS)
    cstate.stresses = g_config.stress + g_config.nconf;
#endif /* STRESS */
#if defined(CONTRIB)
    g_config.have_contrib = 0;
    g_config.have_contrib_box = 0;
#endif /* CONTRIB */

    do
    {
      res = fgets(buffer, 1024, config_file);

      if ((ptr = strchr(res, '\n')) != NULL)
        *ptr = '\0';

      line++;

      /* read the box vectors: x */
      if (res[1] == 'X') {
        if (sscanf(res + 3, "%lf %lf %lf\n", &cstate.box_x.x, &cstate.box_x.y, &cstate.box_x.z) == 3) {
          cstate.h_box_x = 1;
        if (g_param.global_cell_scale != 1.0) {
            cstate.box_x.x *= g_param.global_cell_scale;
            cstate.box_x.y *= g_param.global_cell_scale;
            cstate.box_x.z *= g_param.global_cell_scale;
          }
        } else
          error(1, "%s: Error in box vector x, line %d\n", filename, line);
      }

      /* read the box vectors: y */
      else if (res[1] == 'Y') {
        if (sscanf(res + 3, "%lf %lf %lf\n", &cstate.box_y.x, &cstate.box_y.y, &cstate.box_y.z) == 3) {
          cstate.h_box_y = 1;
        if (g_param.global_cell_scale != 1.0) {
          cstate.box_y.x *= g_param.global_cell_scale;
          cstate.box_y.y *= g_param.global_cell_scale;
          cstate.box_y.z *= g_param.global_cell_scale;
          }
        } else
          error(1, "%s: Error in box vector y, line %d\n", filename, line);
      }

      /* read the box vectors: z */
      else if (res[1] == 'Z') {
        if (sscanf(res + 3, "%lf %lf %lf\n", &cstate.box_z.x, &cstate.box_z.y, &cstate.box_z.z) == 3) {
          cstate.h_box_z = 1;
        if (g_param.global_cell_scale != 1.0) {
          cstate.box_z.x *= g_param.global_cell_scale;
          cstate.box_z.y *= g_param.global_cell_scale;
          cstate.box_z.z *= g_param.global_cell_scale;
          }
        } else
          error(1, "%s: Error in box vector z, line %d\n", filename, line);
      }
#if defined(CONTRIB)
      /* box of contributing particles: origin */
      else if (strncmp(res + 1, "B_O", 3) == 0) {
        if (cstate.have_contrib_box == 1) {
          error(0, "There can only be one box of contributing atoms\n");
          error(1, "This occured in %s on line %d", filename, line);
        }
        if (sscanf(res + 5, "%lf %lf %lf\n", &cstate.cbox_o.x, &cstate.cbox_o.y, &cstate.cbox_o.z) == 3) {
          cstate.have_contrib_box = 1;
          if (g_param.global_cell_scale != 1.0) {
            cstate.cbox_o.x *= g_param.global_cell_scale;
            cstate.cbox_o.y *= g_param.global_cell_scale;
            cstate.cbox_o.z *= g_param.global_cell_scale;
          }
        } else
          error(1, "%s: Error in box of contributing atoms, line %d\n",
                filename, line);
      }

      /* box of contributing particles: a */
      else if (strncmp(res + 1, "B_A", 3) == 0) {
        if (sscanf(res + 5, "%lf %lf %lf\n", &cbox_a.x, &cbox_a.y,
                    &cbox_a.z) == 3) {
          have_contrib++;
          if (global_cell_scale != 1.0) {
            cbox_a.x *= global_cell_scale;
            cbox_a.y *= global_cell_scale;
            cbox_a.z *= global_cell_scale;
          }
        } else
          error(1, "%s: Error in box of contributing atoms, line %d\n",
                filename, line);
      }

      /* box of contributing particles: b */
      else if (strncmp(res + 1, "B_B", 3) == 0) {
        if (sscanf(res + 5, "%lf %lf %lf\n", &cbox_b.x, &cbox_b.y,
                    &cbox_b.z) == 3) {
          have_contrib++;
          if (global_cell_scale != 1.0) {
            cbox_b.x *= global_cell_scale;
            cbox_b.y *= global_cell_scale;
            cbox_b.z *= global_cell_scale;
          }
        } else
          error(1, "%s: Error in box of contributing atoms, line %d\n",
                filename, line);
      }

      /* box of contributing particles: c */
      else if (strncmp(res + 1, "B_C", 3) == 0) {
        if (sscanf(res + 5, "%lf %lf %lf\n", &cbox_c.x, &cbox_c.y,
                    &cbox_c.z) == 3) {
          have_contrib++;
          if (global_cell_scale != 1.0) {
            cbox_c.x *= global_cell_scale;
            cbox_c.y *= global_cell_scale;
            cbox_c.z *= global_cell_scale;
          }
        } else
          error(1, "%s: Error in box of contributing atoms, line %d\n",
                filename, line);
      }

      /* sphere of contributing particles */
      else if (strncmp(res + 1, "B_S", 3) == 0) {
        sphere_centers = (vector *)realloc(sphere_centers,
                                            (n_spheres + 1) * sizeof(vector));
        r_spheres =
            (double *)realloc(r_spheres, (n_spheres + 1) * sizeof(double));
        if (sscanf(res + 5, "%lf %lf %lf %lf\n", &sphere_centers[n_spheres].x,
                    &sphere_centers[n_spheres].y, &sphere_centers[n_spheres].z,
                    &r_spheres[n_spheres]) == 4) {
          n_spheres++;
          if (global_cell_scale != 1.0) {
            sphere_centers[n_spheres].x *= global_cell_scale;
            sphere_centers[n_spheres].y *= global_cell_scale;
            sphere_centers[n_spheres].z *= global_cell_scale;
            r_spheres[n_spheres] *= global_cell_scale;
          }
        } else
          error(1, "%s: Error in sphere of contributing atoms, line %d\n",
                filename, line);
      }
#endif /* CONTRIB */

      /* cohesive energy */
      else if (res[1] == 'E') {
        if (sscanf(res + 3, "%lf\n", &(g_config.coheng[g_config.nconf])) == 1)
          cstate.have_energy = 1;
        else
          error(1, "%s: Error in energy on line %d\n", filename, line);
      }

      /* configuration weight */
      else if (res[1] == 'W') {
        if (sscanf(res + 3, "%lf\n", &(g_config.conf_weight[g_config.nconf])) != 1)
          error(1, "%s: Error in configuration weight on line %d\n", filename,
                line);
          if (g_config.conf_weight[g_config.nconf] < 0.0)
          error(1, "%s: The configuration weight is negative on line %d\n",
                filename, line);
      }

      /* chemical elements */
      else if (res[1] == 'C') {
        fgetpos(config_file, &filepos);
        read_chemical_elements(res, &cstate);
        fsetpos(config_file, &filepos);
      }

#if defined(STRESS)
      /* read stress */
      else if (res[1] == 'S') {
        if (sscanf(res + 3, "%lf %lf %lf %lf %lf %lf\n", &(cstate.stresses->xx),
          &(cstate.stresses->yy), &(cstate.stresses->zz), &(cstate.stresses->xy),
                   &(cstate.stresses->yz), &(cstate.stresses->zx)) == 6)
          cstate.have_stress = 1;
        else
          error(1, "Error in stress tensor on line %d\n", line);
      }
#endif /* STRESS */

      else if (res[1] != '#' && res[1] != 'F') {
        warning("Unknown header line in %s detected:\n", filename);
        warning("Line %d : %s\n", line, res);
      }

    } while (res[1] != 'F');

    if (0 == cstate.have_energy)
      error(1, "%s: missing energy in configuration %d!", filename, g_config.nconf);

    if (!(cstate.h_box_x && cstate.h_box_y && cstate.h_box_z))
      error(1, "Incomplete box vectors for config %d!", g_config.nconf);

#if defined(CONTRIB)
    if (have_contrib_box && have_contrib != 4)
      error(1, "Incomplete box of contributing atoms for config %d!", g_config.nconf);
#endif /* CONTRIB */

#if defined(STRESS)
    g_config.usestress[g_config.nconf] = cstate.have_stress; /* no stress tensor available */
#endif /* STRESS */

#if defined(STRESS)
    if (cstate.have_stress)
      with_stresses++;
#endif /* STRESS */

    if (cstate.use_force)
      with_forces++;

    g_config.volume[g_config.nconf] = make_box(&cstate);

    /* read the atoms */
    for (i = 0; i < cstate.atom_count; i++)
    {
      atom_t* atom = g_config.atoms + g_config.natoms + i;

      if (7 > fscanf(config_file, "%d %lf %lf %lf %lf %lf %lf\n", &(atom->type),
                     &(atom->pos.x), &(atom->pos.y), &(atom->pos.z),
                     &(atom->force.x), &(atom->force.y), &(atom->force.z)))
        error(1, "Corrupt configuration file on line %d\n", line + 1);

      line++;

      if (g_param.global_cell_scale != 1.0) {
        atom->pos.x *= g_param.global_cell_scale;
        atom->pos.y *= g_param.global_cell_scale;
        atom->pos.z *= g_param.global_cell_scale;
      }

      if (atom->type >= g_param.ntypes || atom->type < 0)
        error( 1, "Corrupt configuration file on line %d: Incorrect atom type (%d)\n", line, atom->type);

      atom->absforce = sqrt(dsquare(atom->force.x) + dsquare(atom->force.y) + dsquare(atom->force.z));

      atom->conf = g_config.nconf;

#if defined(CONTRIB)
      if (have_contrib_box || n_spheres != 0)
        atom->contrib = does_contribute(atom->pos);
      else
        atom->contrib = 1;
#endif // CONTRIB

      g_config.na_type[g_config.nconf][atom->type] += 1;
      max_atom_type = MAX(max_atom_type, atom->type);
    }

    /* check cell size */
    /* inverse height in direction */
    vector iheight;

    iheight.x = sqrt(SPROD(cstate.tbox_x, cstate.tbox_x));
    iheight.y = sqrt(SPROD(cstate.tbox_y, cstate.tbox_y));
    iheight.z = sqrt(SPROD(cstate.tbox_z, cstate.tbox_z));

    if ((ceil(g_config.rcutmax * iheight.x) > 30000) ||
        (ceil(g_config.rcutmax * iheight.y) > 30000) ||
        (ceil(g_config.rcutmax * iheight.z) > 30000))
      error(1, "Very bizarre small cell size - aborting");

    cstate.cell_scale.x = (int)ceil(g_config.rcutmax * iheight.x);
    cstate.cell_scale.y = (int)ceil(g_config.rcutmax * iheight.y);
    cstate.cell_scale.z = (int)ceil(g_config.rcutmax * iheight.z);

    if (cstate.cell_scale.x > 1 || cstate.cell_scale.y > 1 || cstate.cell_scale.z > 1)
    {
      // TODO
      warning("The box size of at least one configuration is smaller than the "
        "cutoff distance.\n");
      warning("Using additional periodic images for energy and force "
        "calculations.\n");
    }

#ifdef DEBUG
//     fprintf(stderr, "\nChecking cell size for configuration %d:\n", nconf);
//     fprintf(stderr, "Box dimensions:\n");
//     fprintf(stderr, "     %10.6f %10.6f %10.6f\n", box_x.x, box_x.y, box_x.z);
//     fprintf(stderr, "     %10.6f %10.6f %10.6f\n", box_y.x, box_y.y, box_y.z);
//     fprintf(stderr, "     %10.6f %10.6f %10.6f\n", box_z.x, box_z.y, box_z.z);
//     fprintf(stderr, "Box normals:\n");
//     fprintf(stderr, "     %10.6f %10.6f %10.6f\n", tbox_x.x, tbox_x.y,
//             tbox_x.z);
//     fprintf(stderr, "     %10.6f %10.6f %10.6f\n", tbox_y.x, tbox_y.y,
//             tbox_y.z);
//     fprintf(stderr, "     %10.6f %10.6f %10.6f\n", tbox_z.x, tbox_z.y,
//             tbox_z.z);
//     fprintf(stderr, "Box heights:\n");
//     fprintf(stderr, "     %10.6f %10.6f %10.6f\n", 1.0 / iheight.x,
//             1.0 / iheight.y, 1.0 / iheight.z);
//     fprintf(stderr, "Potential range:  %f\n", rcutmax);
//     fprintf(stderr, "Periodic images needed: %d %d %d\n\n",
//             2 * cell_scale[0] + 1, 2 * cell_scale[1] + 1,
//             2 * cell_scale[2] + 1);
#endif /* DEBUG */

    init_neighbors(&cstate, mindist);

    init_angles(&cstate);

    // Addition by AI 17/7/2015
    #ifdef LMP
       g_config.lattice = (lattice_t *)realloc(g_config.lattice,(g_config.nconf+1)*sizeof(lattice_t));
       if(g_config.lattice==NULL){error(1, "Cannot allocate memory for lattice vector");}
       g_config.lattice[g_config.nconf].xx = g_config.box_x.x;
       g_config.lattice[g_config.nconf].xy = g_config.box_x.y;
       g_config.lattice[g_config.nconf].xz = g_config.box_x.z;
       g_config.lattice[g_config.nconf].yx = g_config.box_y.x;
       g_config.lattice[g_config.nconf].yy = g_config.box_y.y;
       g_config.lattice[g_config.nconf].yz = g_config.box_y.z;
       g_config.lattice[g_config.nconf].zx = g_config.box_z.x;
       g_config.lattice[g_config.nconf].zy = g_config.box_z.y;
       g_config.lattice[g_config.nconf].zz = g_config.box_z.z;
    #endif
    // End addition

    /* increment natoms and configuration number */
    g_config.natoms += cstate.atom_count;
    g_config.nconf++;

  } while (!feof(config_file));

  /* close config file */
  fclose(config_file);

  /* the calculation of the neighbor lists is now complete */
  printf("done\n");

  /* calculate the total number of the atom types */
  g_config.na_type = (int **)realloc(g_config.na_type, (g_config.nconf + 1) * sizeof(int *));

  reg_for_free(g_config.na_type, "na_type");

  if (NULL == g_config.na_type)
    error(1, "Cannot allocate memory for na_type");

  g_config.na_type[g_config.nconf] = (int *)malloc(g_param.ntypes * sizeof(int));

  reg_for_free(g_config.na_type[g_config.nconf], "na_type[%d]", g_config.nconf);

  for (i = 0; i < g_param.ntypes; i++)
    g_config.na_type[g_config.nconf][i] = 0;

  for (i = 0; i < g_config.nconf; i++)
    for (j = 0; j < g_param.ntypes; j++)
      g_config.na_type[g_config.nconf][j] += g_config.na_type[i][j];

  /* print diagnostic message */
  printf("\nRead %d configurations (%d with forces, %d with stresses)\n", g_config.nconf, with_forces, with_stresses);
  printf("with a total of %d atoms (", g_config.natoms);
  for (i = 0; i < g_param.ntypes; i++) {
    if (cstate.have_elements)
      printf("%d %s (%.2f%%)", g_config.na_type[g_config.nconf][i], g_config.elements[i],
             100.0 * g_config.na_type[g_config.nconf][i] / g_config.natoms);
    else
      printf("%d type %d (%.2f%%)", g_config.na_type[g_config.nconf][i], i,
             100.0 * g_config.na_type[g_config.nconf][i] / g_config.natoms);
      if (i != (g_param.ntypes - 1))
      printf(", ");
  }
  printf(").\n");

  /* be pedantic about too large g_param.ntypes */
  if ((max_atom_type + 1) < g_param.ntypes) {
    error(0, "There are less than %d atom types in your configurations!\n", g_param.ntypes);
    error(1, "Please adjust \"ntypes\" in your parameter file.");
  }

  reg_for_free(g_config.atoms, "atoms");
  reg_for_free(g_config.coheng, "coheng");
  reg_for_free(g_config.conf_weight, "conf_weight");
  reg_for_free(g_config.volume, "volume");
#ifdef STRESS
  reg_for_free(g_config.stress, "stress");
#endif /* STRESS */
  reg_for_free(g_config.inconf, "inconf");
  reg_for_free(g_config.cnfstart, "cnfstart");
  reg_for_free(g_config.useforce, "useforce");
#ifdef STRESS
  reg_for_free(g_config.usestress, "usestress");
#endif /* STRESS */
#ifdef CONTRIB
  if (n_spheres > 0) {
    reg_for_free(r_spheres, "sphere radii");
    reg_for_free(sphere_centers, "sphere centers");
  }
#endif /* CONTRIB */

  /* mdim is the dimension of the force vector:
     - 3*natoms forces
     - nconf cohesive energies,
     - 6*nconf stress tensor components */
  g_calc.mdim = 3 * g_config.natoms + g_config.nconf;
#ifdef STRESS
  g_calc.mdim += 6 * g_config.nconf;
#endif /* STRESS */

/* mdim has additional components for EAM-like potentials */
#if defined EAM || defined ADP || defined MEAM
  g_calc.mdim += g_config.nconf;      /* nconf limiting constraints */
  g_calc.mdim += 2 * g_param.ntypes; /* g_param.ntypes dummy constraints */
#ifdef TBEAM
  g_calc.mdim += 2 * g_param.ntypes; /* additional dummy constraints for s-band */
#endif                /* TBEAM */
#endif                /* EAM || ADP || MEAM */

/* mdim has additional components for analytic potentials */
#ifdef APOT
  /* 1 slot for each analytic parameter -> punishment */
  g_calc.mdim += g_pot.opt_pot.idxlen;
  /* 1 slot for each analytic potential -> punishment */
  g_calc.mdim += g_pot.apot_table.number + 1;
#endif /* APOT */

  /* copy forces into single vector */
  if (NULL == (g_config.force_0 = (double *)malloc(g_calc.mdim * sizeof(double))))
    error(1, "Cannot allocate force vector");
  reg_for_free(g_config.force_0, "force_0");

  k = 0;
  for (i = 0; i < g_config.natoms; i++) { /* first forces */
    g_config.force_0[k++] = g_config.atoms[i].force.x;
    g_config.force_0[k++] = g_config.atoms[i].force.y;
    g_config.force_0[k++] = g_config.atoms[i].force.z;
  }
  for (i = 0; i < g_config.nconf; i++) { /* then cohesive energies */
    g_config.force_0[k++] = g_config.coheng[i];
  }

#ifdef STRESS
  for (i = 0; i < g_config.nconf; i++) { /* then stresses */
    if (g_config.usestress[i]) {
      g_config.force_0[k++] = g_config.stress[i].xx;
      g_config.force_0[k++] = g_config.stress[i].yy;
      g_config.force_0[k++] = g_config.stress[i].zz;
      g_config.force_0[k++] = g_config.stress[i].xy;
      g_config.force_0[k++] = g_config.stress[i].yz;
      g_config.force_0[k++] = g_config.stress[i].zx;
    } else {
      for (j = 0; j < 6; j++)
        g_config.force_0[k++] = 0.0;
    }
  }
#endif /* STRESS */

#if defined EAM || defined ADP || defined MEAM
  for (i = 0; i < g_config.nconf; i++)
    g_config.force_0[k++] = 0.0; /* punishment rho out of bounds */
  for (i = 0; i < 2 * g_param.ntypes; i++)
    g_config.force_0[k++] = 0.0; /* constraint on U(n=0):=0 */
#ifdef TBEAM
  for (i = 0; i < 2 * g_param.ntypes; i++)
    g_config.force_0[k++] = 0.0; /* constraint on U(n=0):=0 for s-band */
#endif                  /* TBEAM */
#endif                  /* EAM || ADP || MEAM */

#ifdef APOT
  for (i = 0; i < g_pot.opt_pot.idxlen; i++)
    g_config.force_0[k++] = 0.0; /* punishment for individual parameters */
  for (i = 0; i <= g_pot.apot_table.number; i++)
    g_config.force_0[k++] = 0.0; /* punishment for potential functions */
#endif                  /* APOT */

  write_pair_distribution_file();

/* assign correct distances to different tables */
#ifdef APOT
  double min = 10.0;

  /* pair potentials */
  for (i = 0; i < g_param.ntypes; i++) {
    for (j = 0; j < g_param.ntypes; j++) {
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
#if defined EAM || defined ADP || defined MEAM
for (i = g_calc.paircol; i < g_calc.paircol + g_param.ntypes; i++) {
    g_pot.apot_table.begin[i] = min * 0.95;
    g_pot.opt_pot.begin[i] = min * 0.95;
    g_pot.calc_pot.begin[i] = min * 0.95;
  }
#ifdef TBEAM
for (i = g_calc.paircol + 2 * g_param.ntypes; i < g_calc.paircol + 3 * g_param.ntypes; i++) {
    g_pot.apot_table.begin[i] = min * 0.95;
    g_pot.opt_pot.begin[i] = min * 0.95;
    g_pot.calc_pot.begin[i] = min * 0.95;
  }
#endif /* TBEAM */
#endif /* EAM || ADP || MEAM */

/* dipole and quadrupole functions */
#ifdef ADP
  for (i = 0; i < g_calc.paircol; i++) {
    j = g_calc.paircol + 2 * g_param.ntypes + i;
    g_pot.apot_table.begin[j] = min * 0.95;
    g_pot.opt_pot.begin[j] = min * 0.95;
    g_pot.calc_pot.begin[j] = min * 0.95;
    j = 2 * g_calc.paircol + 2 * g_param.ntypes + i;
    g_pot.apot_table.begin[j] = min * 0.95;
    g_pot.opt_pot.begin[j] = min * 0.95;
    g_pot.calc_pot.begin[j] = min * 0.95;
  }
#endif /* ADP */

#ifdef MEAM
  /* f_ij */
  for (i = 0; i < g_calc.paircol; i++) {
    j = g_calc.paircol + 2 * g_param.ntypes + i;
    g_pot.apot_table.begin[j] = min * 0.95;
    g_pot.opt_pot.begin[j] = min * 0.95;
    g_pot.calc_pot.begin[j] = min * 0.95;
  }
  /* g_i */
  /* g_i takes cos(theta) as an argument, so we need to tabulate it only
     in the range of [-1:1]. Actually we use [-1.1:1.1] to be safe. */
  for (i = 0; i < g_param.ntypes; i++) {
    j = 2 * g_calc.paircol + 2 * g_param.ntypes + i;
    g_pot.apot_table.begin[j] = -1.1;
    g_pot.opt_pot.begin[j] = -1.1;
    g_pot.calc_pot.begin[j] = -1.1;
    g_pot.apot_table.end[j] = 1.1;
    g_pot.opt_pot.end[j] = 1.1;
    g_pot.calc_pot.end[j] = 1.1;
  }
#endif /* MEAM */

  /* recalculate step, invstep and xcoord for new tables */
  for (i = 0; i < g_pot.calc_pot.ncols; i++) {
    g_pot.calc_pot.step[i] = (g_pot.calc_pot.end[i] - g_pot.calc_pot.begin[i]) / (APOT_STEPS - 1);
    g_pot.calc_pot.invstep[i] = 1.0 / g_pot.calc_pot.step[i];
    for (j = 0; j < APOT_STEPS; j++) {
      int index = i * APOT_STEPS + (i + 1) * 2 + j;
      g_pot.calc_pot.xcoord[index] = g_pot.calc_pot.begin[i] + j * g_pot.calc_pot.step[i];
    }
  }

  update_slots();
#endif /* APOT */

  print_minimal_distances_matrix(mindist);

  free(mindist);

  return;
}

/****************************************************************
 *
 *  create_memory_for_config
 *      bla bla
 *
 ****************************************************************/

void create_memory_for_config(config_state* cstate)
{
  int i = 0;

  /* increase memory for this many additional atoms */
  g_config.atoms = (atom_t *)realloc(g_config.atoms, (g_config.natoms + cstate->atom_count) * sizeof(atom_t));

  if (g_config.atoms == NULL)
    error(1, "Cannot allocate memory for atoms");

  for (i = 0; i < cstate->atom_count; i++) {
    g_config.atoms[g_config.natoms + i].neigh = (neigh_t *)malloc(sizeof(neigh_t));
    if (g_config.atoms[g_config.natoms + i].neigh == NULL)
      error(1, "Cannot allocate memory for neighbors");
    reg_for_free(g_config.atoms[g_config.natoms + i].neigh, "test neigh");
  }

  g_config.coheng = (double *)realloc(g_config.coheng, (g_config.nconf + 1) * sizeof(double));

  if (g_config.coheng == NULL)
    error(1, "Cannot allocate memory for cohesive energy");

  g_config.conf_weight = (double *)realloc(g_config.conf_weight, (g_config.nconf + 1) * sizeof(double));

  if (g_config.conf_weight == NULL)
    error(1, "Cannot allocate memory for configuration weights");
  else
    g_config.conf_weight[g_config.nconf] = 1.0;

  g_config.volume = (double *)realloc(g_config.volume, (g_config.nconf + 1) * sizeof(double));

  if (g_config.volume == NULL)
    error(1, "Cannot allocate memory for volume");

#if defined(STRESS)
  g_config.stress = (sym_tens *)realloc(g_config.stress, (g_config.nconf + 1) * sizeof(sym_tens));

  if (g_config.stress == NULL)
    error(1, "Cannot allocate memory for stress");
#endif /* STRESS */

  g_config.inconf = (int *)realloc(g_config.inconf, (g_config.nconf + 1) * sizeof(int));

  if (g_config.inconf == NULL)
    error(1, "Cannot allocate memory for atoms in conf");

  g_config.cnfstart = (int *)realloc(g_config.cnfstart, (g_config.nconf + 1) * sizeof(int));

  if (g_config.cnfstart == NULL)
    error(1, "Cannot allocate memory for start of conf");

  g_config.useforce = (int *)realloc(g_config.useforce, (g_config.nconf + 1) * sizeof(int));

  if (g_config.useforce == NULL)
    error(1, "Cannot allocate memory for useforce");

#if defined(STRESS)
  g_config.usestress = (int *)realloc(g_config.usestress, (g_config.nconf + 1) * sizeof(int));

  if (g_config.usestress == NULL)
    error(1, "Cannot allocate memory for usestress");
#endif /* STRESS */

  g_config.na_type = (int **)realloc(g_config.na_type, (g_config.nconf + 2) * sizeof(int *));

  if (g_config.na_type == NULL)
    error(1, "Cannot allocate memory for na_type");

  g_config.na_type[g_config.nconf] = (int *)malloc(g_param.ntypes * sizeof(int));

  reg_for_free(g_config.na_type[g_config.nconf], "na_type[%d]", g_config.nconf);

  if (g_config.na_type[g_config.nconf] == NULL)
    error(1, "Cannot allocate memory for na_type");
}

/****************************************************************
 *
 *  read chemical elements from config file
 *
 ****************************************************************/

void read_chemical_elements(char* psrc, config_state* cstate)
{
  int i = 0;
  char const* pchar = strtok(psrc + 3, " \t");

  if (!cstate->have_elements)
  {
    while (pchar != NULL && i < g_param.ntypes)
    {
      strcpy(g_config.elements[i], pchar);
      pchar = strtok(NULL, " \t");
      i++;
      cstate->num_fixed_elements++;
    }
    cstate->have_elements = 1;
  }
  else
  {
    while (pchar != NULL && i < g_param.ntypes)
    {
      if (strcmp(pchar, g_config.elements[i]) != 0)
      {
        if (atoi(g_config.elements[i]) == i && i > cstate->num_fixed_elements)
        {
          strcpy(g_config.elements[i], pchar);
          cstate->num_fixed_elements++;
        } else {
          error(0, "Mismatch found in configuration %d, line %d.\n", cstate->config, cstate->line);
          error(0, "Expected element >> %s << but found element >> %s <<.\n", g_config.elements[i], pchar);
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
 *
 *  init_atom
 *
 ****************************************************************/

void init_atom_memory(atom_t *atom)
{
  atom->type = 0;
  atom->num_neigh = 0;
  atom->pos.x = 0.0;
  atom->pos.y = 0.0;
  atom->pos.z = 0.0;
  atom->force.x = 0.0;
  atom->force.y = 0.0;
  atom->force.z = 0.0;
  atom->absforce = 0.0;
  atom->conf = 0;

  #ifdef CONTRIB
  atom->contrib = 0;
  #endif /* CONTRIB */

  #if defined EAM || defined ADP || defined MEAM
  atom->rho = 0.0;
  atom->gradF = 0.0;
  #ifdef TBEAM
  atom->rho_s = 0.0;
  atom->gradF_s = 0.0;
  #endif /* TBEAM */
  #endif /* EAM || ADP || MEAM */

  #ifdef ADP
  atom->mu.x = 0.0;
  atom->mu.y = 0.0;
  atom->mu.z = 0.0;
  atom->lambda.xx = 0.0;
  atom->lambda.yy = 0.0;
  atom->lambda.zz = 0.0;
  atom->lambda.xy = 0.0;
  atom->lambda.yz = 0.0;
  atom->lambda.zx = 0.0;
  atom->nu = 0.0;
  #endif /* ADP */

  #ifdef DIPOLE
  atom->E_stat.x = 0.0;
  atom->E_stat.y = 0.0;
  atom->E_stat.z = 0.0;
  atom->p_sr.x = 0.0;
  atom->p_sr.y = 0.0;
  atom->p_sr.z = 0.0;
  atom->E_ind.x = 0.0;
  atom->E_ind.y = 0.0;
  atom->E_ind.z = 0.0;
  atom->p_ind.x = 0.0;
  atom->p_ind.y = 0.0;
  atom->p_ind.z = 0.0;
  atom->E_old.x = 0.0;
  atom->E_old.y = 0.0;
  atom->E_old.z = 0.0;
  atom->E_tot.x = 0.0;
  atom->E_tot.y = 0.0;
  atom->E_tot.z = 0.0;
  #endif /* DIPOLE */

  #ifdef THREEBODY
  atom->num_angles = 0;
  #ifdef MEAM
  atom->rho_eam = 0.0;
  #endif /* MEAM */
  #endif /* MANYBODY */

  atom->neigh = NULL;
  #ifdef THREEBODY
  atom->angle_part = NULL;
  #endif /* THREEBODY */
}

/****************************************************************
 *
 *  init_neigh
 *
 ****************************************************************/

void init_neigh_memory(neigh_t *neigh)
{
  int i = 0;

  neigh->type = 0;
  neigh->nr = 0;
  neigh->r = 0.0;
  neigh->r2 = 0.0;
  neigh->inv_r = 0.0;
  neigh->dist.x = 0.0;
  neigh->dist.y = 0.0;
  neigh->dist.z = 0.0;
  neigh->dist_r.x = 0.0;
  neigh->dist_r.y = 0.0;
  neigh->dist_r.z = 0.0;
  for (i = 0; i < SLOTS; i++) {
    neigh->slot[i] = 0;
    neigh->shift[i] = 0.0;
    neigh->step[i] = 0.0;
    neigh->col[i] = 0;
  }

  #ifdef ADP
  neigh->sqrdist.xx = 0.0;
  neigh->sqrdist.yy = 0.0;
  neigh->sqrdist.zz = 0.0;
  neigh->sqrdist.xy = 0.0;
  neigh->sqrdist.yz = 0.0;
  neigh->sqrdist.zx = 0.0;
  neigh->u_val = 0.0;
  neigh->u_grad = 0.0;
  neigh->w_val = 0.0;
  neigh->w_grad = 0.0;
  #endif /* APOT */

  #ifdef COULOMB
  neigh->fnval_el = 0.0;
  neigh->grad_el = 0.0;
  neigh->ggrad_el = 0.0;
  #endif /* COULOMB */

  #ifdef THREEBODY
  neigh->f = 0.0;
  neigh->df = 0.0;
  neigh->ijk_start = 0;
  #endif /* THREEBODY */

  #ifdef MEAM
  neigh->drho = 0.0;
  #endif /* MEAM */

  #ifdef TERSOFF
  neigh->dzeta.x = 0.0;
  neigh->dzeta.y = 0.0;
  neigh->dzeta.z = 0.0;
  #endif /* TERSOFF */
}

#ifdef THREEBODY
void init_angle(angle_t * angle)
{
  angle->cos = 0.0;

  #ifdef MEAM
  angle->slot = 0;
  angle->shift = 0.0;
  angle->step = 0.0;
  angle->g = 0.0;
  angle->dg = 0.0;
  #endif /* MEAM */
}
#endif /* THREEBODY */

/****************************************************************
 *
 *  init_neighbors
 *      compute box transformation matrix
 *
 ****************************************************************/

void init_neighbors(config_state* cstate, double* mindist)
{
  int i = 0;
  int j = 0;
  int ix = 0;
  int iy = 0;
  int iz = 0;

  vector d;
  vector dd;

  /* compute the neighbor table */
  for (i = g_config.natoms; i < g_config.natoms + cstate->atom_count; i++)
  {
    g_config.atoms[i].num_neigh = 0;
    /* loop over all atoms for threebody interactions */
#if defined(THREEBODY)
    for (j = g_config.natoms; j < g_config.natoms + cstate->atom_count; j++)
#else
    for (j = i; j < g_config.natoms + cstate->atom_count; j++)
#endif /* THREEBODY */
    {
      d.x = g_config.atoms[j].pos.x - g_config.atoms[i].pos.x;
      d.y = g_config.atoms[j].pos.y - g_config.atoms[i].pos.y;
      d.z = g_config.atoms[j].pos.z - g_config.atoms[i].pos.z;
      for (ix = -cstate->cell_scale.x; ix <= cstate->cell_scale.x; ix++)
      {
        for (iy = -cstate->cell_scale.y; iy <= cstate->cell_scale.y; iy++)
        {
          for (iz = -cstate->cell_scale.z; iz <= cstate->cell_scale.z; iz++)
          {
            if ((i == j) && (ix == 0) && (iy == 0) && (iz == 0))
              continue;
            dd.x = d.x + ix * cstate->box_x.x + iy * cstate->box_y.x + iz * cstate->box_z.x;
            dd.y = d.y + ix * cstate->box_x.y + iy * cstate->box_y.y + iz * cstate->box_z.y;
            dd.z = d.z + ix * cstate->box_x.z + iy * cstate->box_y.z + iz * cstate->box_z.z;
            double r = sqrt(SPROD(dd, dd));
            int type1 = g_config.atoms[i].type;
            int type2 = g_config.atoms[j].type;

            if (r <= g_config.rcut[type1 * g_param.ntypes + type2])
            {
              int short_distance = 0;

              if (r <= g_config.rmin[type1 * g_param.ntypes + type2])
              {
                warning("Configuration %i: Distance %i\n", cstate->config, r);
                warning(" atom %d (type %d) at pos: %f %f %f\n", i - g_config.natoms, type1, g_config.atoms[i].pos.x, g_config.atoms[i].pos.y, g_config.atoms[i].pos.z);
                warning(" atom %d (type %d) at pos: %f %f %f\n", j - g_config.natoms, type2, dd.x, dd.y, dd.z);
                short_distance = 1;
              }
              g_config.atoms[i].neigh = (neigh_t *)realloc(g_config.atoms[i].neigh,(g_config.atoms[i].num_neigh + 1) * sizeof(neigh_t));
              dd.x /= r;
              dd.y /= r;
              dd.z /= r;
              int k = g_config.atoms[i].num_neigh++;
              init_neigh_memory(g_config.atoms[i].neigh + k);
              g_config.atoms[i].neigh[k].type = type2;
              g_config.atoms[i].neigh[k].nr = j;
              g_config.atoms[i].neigh[k].r = r;
              g_config.atoms[i].neigh[k].r2 = r * r;
              g_config.atoms[i].neigh[k].inv_r = 1.0 / r;
              g_config.atoms[i].neigh[k].dist_r = dd;
              g_config.atoms[i].neigh[k].dist.x = dd.x * r;
              g_config.atoms[i].neigh[k].dist.y = dd.y * r;
              g_config.atoms[i].neigh[k].dist.z = dd.z * r;
              #ifdef ADP
              g_config.atoms[i].neigh[k].sqrdist.xx = dd.x * dd.x * r * r;
              g_config.atoms[i].neigh[k].sqrdist.yy = dd.y * dd.y * r * r;
              g_config.atoms[i].neigh[k].sqrdist.zz = dd.z * dd.z * r * r;
              g_config.atoms[i].neigh[k].sqrdist.yz = dd.y * dd.z * r * r;
              g_config.atoms[i].neigh[k].sqrdist.zx = dd.z * dd.x * r * r;
              g_config.atoms[i].neigh[k].sqrdist.xy = dd.x * dd.y * r * r;
              #endif /* ADP */

              // TODO
//               mindist[col] = MIN(mindist[col], r);

              /* pre-compute index and shift into potential table */

              /* pair potential */
              if (!short_distance)
              {
                int col = (type1 <= type2)
                  ? type1 * g_param.ntypes + type2 - ((type1 * (type1 + 1)) / 2)
                  : type2 * g_param.ntypes + type1 - ((type2 * (type2 + 1)) / 2);
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 0);

                #if defined EAM || defined ADP || defined MEAM
                /* transfer function */
                col = g_calc.paircol + type2;
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 1);

                #ifdef TBEAM
                /* transfer function - d band */
                col = g_calc.paircol + 2 * g_param.ntypes + type2;
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 2);
                #endif /* TBEAM */

                #endif /* EAM || ADP || MEAM */

                #ifdef MEAM
                /* Store slots and stuff for f(r_ij) */
                col = g_calc.paircol + 2 * g_param.ntypes + g_config.atoms[i].neigh[k].col[0];
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 2);
                #endif /* MEAM */

                #ifdef ADP
                /* dipole part */
                col = g_calc.paircol + 2 * g_param.ntypes + g_config.atoms[i].neigh[k].col[0];
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 2);

                /* quadrupole part */
                col = 2 * g_calc.paircol + 2 * g_param.ntypes + g_config.atoms[i].neigh[k].col[0];
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 3);
                #endif /* ADP */

                #ifdef STIWEB
                /* Store slots and stuff for exp. function */
                col = g_calc.paircol + atoms[i].neigh[k].col[0];
                set_neighbor_slot(g_config.atoms[i].neigh + k, col, r, 1);
                #endif /* STIWEB */

              } /* !sh_dist */
            }   /* r < r_cut */
          }     /* loop over images in z direction */
        }       /* loop over images in y direction */
      }         /* loop over images in x direction */
    }           /* second loop over atoms (neighbors) */

    reg_for_free(g_config.atoms[i].neigh, "neighbor table atom %d", i);
  } /* first loop over atoms */
}

/****************************************************************
 *
 *  set_neighbor_slot
 *      compute box transformation matrix
 *
 ****************************************************************/

void set_neighbor_slot(neigh_t* neighbor, int col, double r, int store_slot)
{
  int slot = 0;
  double step = 0.0;
  double istep = 0.0;
  double shift = 0.0;

  if (g_pot.format == 0 || g_pot.format == 3)
  {
    double rr = r - g_pot.calc_pot.begin[col];

    if (rr < 0)
    {
      // TODO: rephrase
      error(0, "The distance %f is smaller than the beginning\n", r);
      error(0, "of the potential #%d (r_begin=%f).\n", col, g_pot.calc_pot.begin[col]);
      fflush(stdout);
      error(1, "Short distance!");
    }

    istep = g_pot.calc_pot.invstep[col];
    slot = (int)(rr * istep);
    shift = (rr - slot * g_pot.calc_pot.step[col]) * istep;
    slot += g_pot.calc_pot.first[col];
    step = g_pot.calc_pot.step[col];
  }
  else /* format == 4 ! */
  {
    int klo = g_pot.calc_pot.first[col];
    int khi = g_pot.calc_pot.last[col];

    /* bisection */
    while (khi - klo > 1)
    {
      slot = (khi + klo) >> 1;
      if (g_pot.calc_pot.xcoord[slot] > r)
        khi = slot;
      else
        klo = slot;
    }

    slot = klo;
    step = g_pot.calc_pot.xcoord[khi] - g_pot.calc_pot.xcoord[klo];
    shift = (r - g_pot.calc_pot.xcoord[klo]) / step;
  }

  /* independent of format - we should be left of last index */
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
 *
 *  init_angles
 *      compute box transformation matrix
 *
 ****************************************************************/

void init_angles(config_state* cstate)
{
  /* compute the angular part */
/* For TERSOFF we create a full neighbor list, for all other potentials only a
 * half list */
#ifdef THREEBODY
    for (i = natoms; i < natoms + count; i++) {
      nnn = atoms[i].num_neigh;
      ijk = 0;
      atoms[i].angle_part = (angle_t *)malloc(sizeof(angle_t));
#ifdef TERSOFF
      for (j = 0; j < nnn; j++) {
#else
      for (j = 0; j < nnn - 1; j++) {
#endif /* TERSOFF */
        atoms[i].neigh[j].ijk_start = ijk;
#ifdef TERSOFF
        for (k = 0; k < nnn; k++) {
          if (j == k)
            continue;
#else
        for (k = j + 1; k < nnn; k++) {
#endif /* TERSOFF */
          atoms[i].angle_part = (angle_t *)realloc(atoms[i].angle_part,
                                                   (ijk + 1) * sizeof(angle_t));
          init_angle(atoms[i].angle_part + ijk);
          ccos = atoms[i].neigh[j].dist_r.x * atoms[i].neigh[k].dist_r.x +
                 atoms[i].neigh[j].dist_r.y * atoms[i].neigh[k].dist_r.y +
                 atoms[i].neigh[j].dist_r.z * atoms[i].neigh[k].dist_r.z;

          atoms[i].angle_part[ijk].cos = ccos;

          col = 2 * g_calc.paircol + 2 * g_param.ntypes + atoms[i].type;
          if (0 == format || 3 == format) {
            if ((fabs(ccos) - 1.0) > 1e-10) {
              printf("%.20f %f %d %d %d\n", ccos, g_pot.calc_pot.begin[col], col,
                     type1, type2);
              fflush(stdout);
              error(1, "cos out of range, it is strange!");
            }
#ifdef MEAM
            istep = g_pot.calc_pot.invstep[col];
            slot = (int)((ccos + 1) * istep);
            shift = ((ccos + 1) - slot * g_pot.calc_pot.step[col]) * istep;
            slot += g_pot.calc_pot.first[col];
            step = g_pot.calc_pot.step[col];

            /* Don't want lower bound spline knot to be final knot or upper
               bound knot will cause trouble since it goes beyond the array */
            if (slot >= g_pot.calc_pot.last[col]) {
              slot--;
              shift += 1.0;
            }
#endif /* !MEAM */
          }
#ifdef MEAM
          atoms[i].angle_part[ijk].shift = shift;
          atoms[i].angle_part[ijk].slot = slot;
          atoms[i].angle_part[ijk].step = step;
#endif /* MEAM */
          ijk++;
        } /* third loop over atoms */
      }   /* second loop over atoms */
      atoms[i].num_angles = ijk;
      reg_for_free(atoms[i].angle_part, "angular part atom %d", i);
    }  /* first loop over atoms */
#endif /* THREEBODY */
}

/****************************************************************
 *
 *  compute box transformation matrix
 *
 ****************************************************************/

double make_box(config_state* cstate)
{
  double volume = 0;

  /* compute tbox_j such that SPROD(box_i,tbox_j) == delta_ij */
  /* first unnormalized */
  cstate->tbox_x = vec_prod(cstate->box_y, cstate->box_z);
  cstate->tbox_y = vec_prod(cstate->box_z, cstate->box_x);
  cstate->tbox_z = vec_prod(cstate->box_x, cstate->box_y);

  /* volume */
  volume = SPROD(cstate->box_x, cstate->tbox_x);

  if (0.0 == volume)
    error(1, "Box edges are parallel\n");

  /* normalization */
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
 *
 *  compute box transformation matrix
 *
 ****************************************************************/

void write_pair_distribution_file ()
{
  int i = 0;
  int j = 0;
  int k = 0;

  if (1 == g_param.write_pair)
  {
    char pairname[255];
    FILE *pairfile;
    #ifdef APOT
    int pair_steps = APOT_STEPS / 2;
    #else
    int pair_steps = 500;
    #endif /* APOT */
    double pair_table[g_calc.paircol * pair_steps];
    double pair_dist[g_calc.paircol];
    int pos, max_count = 0;

    strcpy(pairname, g_files.config);
    strcat(pairname, ".pair");
    pairfile = fopen(pairname, "w");
    fprintf(pairfile, "# radial distribution file for %d potential(s)\n", g_calc.paircol);

    for (i = 0; i < g_calc.paircol * pair_steps; i++)
      pair_table[i] = 0.0;

    for (i = 0; i < g_param.ntypes; i++)
      for (k = 0; k < g_param.ntypes; k++)
        pair_dist[(i <= k) ? i * g_param.ntypes + k - (i * (i + 1) / 2)
        : k * g_param.ntypes + i - (k * (k + 1) / 2)] =
        g_config.rcut[i * g_param.ntypes + k] / pair_steps;

      for (k = 0; k < g_calc.paircol; k++)
      {
        for (i = 0; i < g_config.natoms; i++)
        {
          for (j = 0; j < g_config.atoms[i].num_neigh; j++)
          {
            int col = g_config.atoms[i].neigh[j].col[0];

            if (col == k)
            {
              pos = (int)(g_config.atoms[i].neigh[j].r / pair_dist[k]);
              #ifdef DEBUG
              if (g_config.atoms[i].neigh[j].r <= 1) {
                warning("Short distance (%f) found.\n", g_config.atoms[i].neigh[j].r);
                warning("\tatom=%d neighbor=%d\n", i, j);
              }
              #endif /* DEBUG */
              pair_table[k * pair_steps + pos]++;
              if ((int)pair_table[k * pair_steps + pos] > max_count)
                max_count = (int)pair_table[k * pair_steps + pos];
            }
          }
        }
      }

      for (k = 0; k < g_calc.paircol; k++) {
        for (i = 0; i < pair_steps; i++) {
          pair_table[k * pair_steps + i] /= max_count;
          fprintf(pairfile, "%f %f\n", i * pair_dist[k],
                  pair_table[k * pair_steps + i]);
        }
        if (k != (g_calc.paircol - 1))
          fprintf(pairfile, "\n\n");
      }
      fclose(pairfile);
  }
}

/****************************************************************
 *
 *  print_minimal_distances_matrix
 *      bla
 *      bla
 *
 ****************************************************************/

void print_minimal_distances_matrix(const double *mindist)
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


#ifdef CONTRIB

/****************************************************************
 *
 *  check if the atom does contribute to the error sum
 *
 ****************************************************************/

int does_contribute(vector pos) {
  int i;
  double n_a, n_b, n_c, r;
  vector dist;

  if (have_contrib_box) {
    dist.x = pos.x - cbox_o.x;
    dist.y = pos.y - cbox_o.y;
    dist.z = pos.z - cbox_o.z;
    n_a = SPROD(dist, cbox_a) / SPROD(cbox_a, cbox_a);
    n_b = SPROD(dist, cbox_b) / SPROD(cbox_b, cbox_b);
    n_c = SPROD(dist, cbox_c) / SPROD(cbox_c, cbox_c);
    if (n_a >= 0 && n_a <= 1)
      if (n_b >= 0 && n_b <= 1)
        if (n_c >= 0 && n_c <= 1)
          return 1;
  }

  for (i = 0; i < n_spheres; i++) {
    dist.x = (pos.x - sphere_centers[i].x);
    dist.y = (pos.y - sphere_centers[i].y);
    dist.z = (pos.z - sphere_centers[i].z);
    r = SPROD(dist, dist);
    r = sqrt(r);
    if (r < r_spheres[i])
      return 1;
  }

  return 0;
}

#endif /* CONTRIB */

#if defined(APOT)

/****************************************************************
 *
 * recalculate the slots of the atoms for analytic potential
 *
 ****************************************************************/

void update_slots(void)
{
  int i = 0;
  int j = 0;
  double r = 0.0;

  for (i = 0; i < g_config.natoms; i++)
  {
    for (j = 0; j < g_config.atoms[i].num_neigh; j++)
    {
      r = g_config.atoms[i].neigh[j].r;

      /* update slots for pair potential part, slot 0 */
      update_neighbor_slots(g_config.atoms[i].neigh + j, r, 0);

#if defined EAM || defined ADP || defined MEAM
      /* update slots for eam transfer functions, slot 1 */
      update_neighbor_slots(g_config.atoms[i].neigh + j, r, 1);
#ifdef TBEAM
      /* update slots for tbeam transfer functions, s-band, slot 2 */
      update_neighbor_slots(g_config.atoms[i].neigh + j, r, 2);
#endif /* TBEAM */
#endif /* EAM || ADP || MEAM */

#ifdef MEAM
      /* update slots for MEAM f functions, slot 2 */
      update_neighbor_slots(g_config.atoms[i].neigh + j, r, 2);
#endif /* MEAM */

#ifdef ADP
      /* update slots for adp dipole functions, slot 2 */
      update_neighbor_slots(g_config.atoms[i].neigh + j, r, 2);

      /* update slots for adp quadrupole functions, slot 3 */
      update_neighbor_slots(g_config.atoms[i].neigh + j, r, 3);
#endif /* ADP */

    } /* end loop over all neighbors */
  }   /* end loop over all atoms */

#ifdef THREEBODY
  /* update angular slots */
  for (i = 0; i < natoms; i++) {
    for (j = 0; j < atoms[i].num_angles; j++) {
      rr = atoms[i].angle_part[j].cos + 1.1;
#ifdef MEAM
      col = 2 * g_calc.paircol + 2 * g_param.ntypes + atoms[i].type;
      atoms[i].angle_part[j].slot = (int)(rr * g_pot.calc_pot.invstep[col]);
      atoms[i].angle_part[j].step = g_pot.calc_pot.step[col];
      atoms[i].angle_part[j].shift =
          (rr - atoms[i].angle_part[j].slot * g_pot.calc_pot.step[col]) *
          g_pot.calc_pot.invstep[col];
      /* move slot to the right potential */
      atoms[i].angle_part[j].slot += g_pot.calc_pot.first[col];
#endif /* MEAM */
    }
  }
#endif /* THREEBODY */

#ifdef STIWEB
  g_pot.apot_table.sw.init = 0;
#endif /* STIWEB */

#ifdef TERSOFF
  g_pot.apot_table.tersoff.init = 0;
#endif /* TERSOFF */
}

/****************************************************************
 *
 *  update_neighbor_slots
 *      recalculate the slots of the atoms for analytic potential
 *
 ****************************************************************/

void  update_neighbor_slots(neigh_t* neighbor, double r, int neighbor_slot)
{
  int col = neighbor->col[neighbor_slot];

  if (r < g_pot.calc_pot.end[col])
  {
    double rr = r - g_pot.calc_pot.begin[col];
    neighbor->slot[neighbor_slot] = (int)(rr * g_pot.calc_pot.invstep[col]);
    neighbor->step[neighbor_slot] = g_pot.calc_pot.step[col];
    neighbor->shift[neighbor_slot] =
      (rr - neighbor->slot[neighbor_slot] * g_pot.calc_pot.step[col]) * g_pot.calc_pot.invstep[col];
    /* move slot to the right potential */
    neighbor->slot[neighbor_slot] += g_pot.calc_pot.first[col];
  }
}

#endif /* APOT */
