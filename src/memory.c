/****************************************************************
 *
 * memory.c: potfit memory management header file
 *
 ****************************************************************
 *
 * Copyright 2002-2018 - the potfit development team
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
 *****************************************************************/

#include "potfit.h"

#include "memory.h"

#if defined(__SANITIZE_ADDRESS__)
  #define ENABLE_ASAN 1
#endif

#if defined(__has_feature)
  #if __has_feature(address_sanitizer)
    #define ENABLE_ASAN 1
  #endif
#endif

// forward declarations of helper functions

void init_interaction_name(const char* name);
void free_atom_table(atom_t* patom, int natoms);
void free_neigh_table(neigh_t* pneigh, int nneigh);
#if defined(THREEBODY)
void free_angle_table(angle_t* pangle, int nangle);
#endif  // THREEBODY
void free_pot_table(pot_table_t* ppot);

typedef struct {
  void** pointers;
  int num_pointers;
} potfit_memory;

static potfit_memory g_memory;

/****************************************************************
 *
 *  Malloc:
 *    allocate memory, initialize it and register it for freeing
 *
 ****************************************************************/

void* Malloc(size_t size)
{
  if (size == 0) {
#if defined(ENABLE_ASAN)
    // trigger asan on purpose
    int* p = NULL;
    *p = 1;
#else
    error(1, "Allocating memory with size 0!\n");
#endif
  }

  void* p = malloc(size);

  if (p == NULL)
    error(1, "Error allocating resources\n");

  memset(p, 0, size);

  g_memory.pointers = (void**)realloc(
      g_memory.pointers, sizeof(void*) * (g_memory.num_pointers + 1));

  if (g_memory.pointers == NULL)
    error(1, "Error allocating resources\n");

  g_memory.pointers[g_memory.num_pointers] = p;
  g_memory.num_pointers++;

  return p;
}

/****************************************************************
 *
 *  Realloc:
 *    reallocate memory and update the freeing array accordingly
 *    initialization is difficult because we don't know to former size
 *
 ****************************************************************/

void* Realloc(void* pvoid, size_t size)
{
  if (size == 0)
    error(1, "Reallocating memory with size 0!\n");

  if (pvoid == NULL)
    return Malloc(size);

  void* temp = realloc(pvoid, size);

  if (temp == NULL)
    error(1, "Error reallocating resources\n");

  for (int i = 0; i < g_memory.num_pointers; i++) {
    if (pvoid == g_memory.pointers[i])
      g_memory.pointers[i] = temp;
  }

  return temp;
}

/****************************************************************
 *
 *  initialize_global_variables
 *
 ****************************************************************/

void initialize_global_variables()
{
  memset(&g_calc, 0, sizeof(g_calc));
  memset(&g_config, 0, sizeof(g_config));
  memset(&g_files, 0, sizeof(g_files));

  g_config.rcutmin = 999.9;

#if defined(KIM)
  memset(&g_kim, 0, sizeof(g_kim));
#endif // KIM

  g_mpi.init_done = 0;
  g_mpi.myid = 0;
  g_mpi.num_cpus = 1;
  g_mpi.firstatom = 0;
  g_mpi.firstconf = 0;
  g_mpi.myatoms = 0;
  g_mpi.myconf = 0;
#if defined(MPI)
  g_mpi.atom_dist = NULL;
  g_mpi.atom_len = NULL;
  g_mpi.conf_dist = NULL;
  g_mpi.conf_len = NULL;
#endif  // MPI

  memset(&g_param, 0, sizeof(g_param));
  g_param.imdpotsteps = 500;
  g_param.lammpspotsteps = 500;
  g_param.sweight = -1.0;
  g_param.global_cell_scale = 1.0;
#if defined(EVO)
  g_param.evo_threshold = 1.0e-6;
#endif  // EVO

  g_pot.interaction_name = NULL;
  g_pot.gradient = NULL;
  g_pot.invar_pot = NULL;
  g_pot.format_type = POTENTIAL_FORMAT_UNKNOWN;
  g_pot.have_invar = 0;
  memset(&g_pot.calc_pot, 0, sizeof(g_pot.calc_pot));
  memset(&g_pot.opt_pot, 0, sizeof(g_pot.opt_pot));
#if defined(APOT)
  g_pot.smooth_pot = NULL;
  g_pot.cp_start = 0;
  g_pot.calc_list = NULL;
  g_pot.compnodelist = NULL;
  memset(&g_pot.apot_table, 0, sizeof(g_pot.apot_table));
#endif  // APOT
#if defined(APOT) || defined(KIM)
  g_pot.global_idx = 0;
  g_pot.global_pot = 0;
  g_pot.have_globals = 0;
#endif // APOT || KIM

  g_memory.pointers = NULL;
  g_memory.num_pointers = 0;

#if defined(PAIR)
  #if !defined(ANG)
    init_interaction_name("PAIR");
  #endif // !ANG
#elif defined(ANG)
  #if !defined(COULOMB)
    init_interaction_name("ANG");
  #elif defined(COULOMB)
    init_interaction_name("ANG_ELSTAT");
  #endif // !COULOMB
#elif defined(EAM) && !defined(COULOMB)
  #if !defined(TBEAM)
    init_interaction_name("EAM");
  #else
    init_interaction_name("TBEAM");
  #endif  // TBEAM
#elif defined(ADP)
  init_interaction_name("ADP");
#elif defined(COULOMB) && !defined(EAM)
  init_interaction_name("ELSTAT");
#elif defined(COULOMB) && defined(EAM)
  init_interaction_name("EAM_ELSTAT");
#elif defined(MEAM)
  init_interaction_name("MEAM");
#elif defined(STIWEB)
  init_interaction_name("STIWEB");
#elif defined(TERSOFF)
  #if defined(TERSOFFMOD)
    init_interaction_name("TERSOFFMOD");
  #else
    init_interaction_name("TERSOFF");
  #endif  // TERSOFFMOD
#elif defined(KIM)
  init_interaction_name("KIM");
#endif  // interaction type
}

/****************************************************************
    init_interaction_name
****************************************************************/

void init_interaction_name(const char* name)
{
  int len = strlen(name);
  char* temp = (char*)Malloc((len + 1) * sizeof(char));
  sprintf(temp, "%s", name);
  g_pot.interaction_name = temp;
}

/****************************************************************
 *
 *  free_allocated_memory -- de-allocate memory of global variables
 *
 ****************************************************************/

void free_allocated_memory()
{
  for (int i = 0; i < g_memory.num_pointers; i++)
    free(g_memory.pointers[i]);

  free(g_memory.pointers);
}

/****************************************************************
 *
 *  free_char_pointer2
 *
 ****************************************************************/

void free_char_pointer2(char** pchar, int len)
{
  if (pchar != NULL) {
    for (int i = 0; i < len; i++)
      free(pchar[i]);
    free(pchar);
  }
}

/****************************************************************
 *
 *  free_double_pointer2
 *
 ****************************************************************/

void free_double_pointer2(double** pdouble, int len)
{
  if (pdouble != NULL) {
    for (int i = 0; i < len; i++)
      free(pdouble[i]);
    free(pdouble);
  }
}

/****************************************************************
 *
 *  free_int_pointer2
 *
 ****************************************************************/

void free_int_pointer2(int** pint, int len)
{
  if (pint != NULL) {
    for (int i = 0; i < len; i++)
      free(pint[i]);
    free(pint);
  }
}

/****************************************************************
 *
 *  free_void_pointer
 *
 ****************************************************************/

void free_void_pointer(void* pvoid)
{
  if (pvoid != NULL)
    free(pvoid);
}

/****************************************************************
 *
 *  free_atom_table
 *
 ****************************************************************/

void free_atom_table(atom_t* patom, int natoms)
{
  if (patom != NULL && natoms > 0) {
    for (int i = 0; i < natoms; i++) {
      free_neigh_table(patom[i].neigh, patom[i].num_neigh);
      free(patom[i].neigh);
#if defined(THREEBODY)
      free_angle_table(patom[i].angle_part, patom[i].num_angles);
      free(patom[i].angle_part);
#endif  // THREEBODY
    }
    free(patom);
  }
}

/****************************************************************
 *  free_neigh_table
 ****************************************************************/

void free_neigh_table(neigh_t* pneigh, int nneigh)
{
  // no dynamically allocated memory
}

#if defined(THREEBODY)

/****************************************************************
 *  free_angle_table
 ****************************************************************/

void free_angle_table(angle_t* pangle, int nangle)
{
  // no dynamically allocated memory
}

#endif  // THREEBODY

/****************************************************************
 *
 *  free_pot_table
 *
 ****************************************************************/

void free_pot_table(pot_table_t* ppot)
{
  free_void_pointer(ppot->begin);
  free_void_pointer(ppot->end);
  free_void_pointer(ppot->step);
  free_void_pointer(ppot->invstep);
  free_void_pointer(ppot->first);
  free_void_pointer(ppot->last);
  free_void_pointer(ppot->xcoord);
  free_void_pointer(ppot->table);
  free_void_pointer(ppot->d2tab);
  free_void_pointer(ppot->idx);
}
