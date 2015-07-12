/****************************************************************
 *
 * memory.c: potfit memory management header file
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
 *****************************************************************/

#include "potfit.h"

#include "memory.h"

// forward declarations of helper functions

void init_interaction_name(const char* name);
void free_char_pointer2(char** pint, int len);
void free_double_pointer2(double** pdouble, int len);
void free_int_pointer2(int** pint, int len);
void free_void_pointer(void* pvoid);
void free_atom_table(atom_t* patom, int natoms);
void free_neigh_table(neigh_t* pneigh, int nneigh);
#if defined(THREEBODY)
void free_angle_table(angle_t* pangle, int nangle);
#endif  // THREEBODY
void free_pot_table(pot_table_t* ppot);
#if defined(APOT)
void free_apot_table(apot_table_t* papot);
#endif

typedef struct
{
  void** pointers;
  int num_pointers;
} potfit_memory;

potfit_memory g_memory;

/****************************************************************
 *
 *  Malloc
 *
 ****************************************************************/

void* Malloc(size_t size)
{
  void* p = malloc(size);

  if (p == NULL)
    error(1, "Error allocating resources");

  memset(p, 0, size);

  g_memory.pointers = (void**)realloc(g_memory.pointers, g_memory.num_pointers + 1);
  g_memory.pointers[g_memory.num_pointers] = p;
  g_memory.num_pointers++;

  return p;
}

/****************************************************************
 *
 *  initialize_global_variables
 *
 ****************************************************************/

void* Realloc(void* pvoid, size_t size)
{
  void* temp = realloc(pvoid, size);

  if (temp == NULL)
    error(1, "Error allocating resources");

  if (pvoid == NULL)
  {
    g_memory.pointers = (void**)realloc(g_memory.pointers, g_memory.num_pointers + 1);
    g_memory.pointers[g_memory.num_pointers] = temp;
    g_memory.num_pointers++;
  }
  else if (temp != pvoid)
  {
    for (int i = 0; i < g_memory.num_pointers; i++)
    {
      if (pvoid == g_memory.pointers[i])
      {
        g_memory.pointers[i] = temp;
      }
    }
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
#endif

  memset(&g_param, 0, sizeof(g_param));

  g_param.global_cell_scale = 1.0;

  g_pot.gradient = NULL;
  g_pot.invar_pot = NULL;
  g_pot.format = -1;
  g_pot.have_invar = 0;
  memset(&g_pot.calc_pot, 0, sizeof(g_pot.calc_pot));
  memset(&g_pot.opt_pot, 0, sizeof(g_pot.opt_pot));
#if defined(APOT)
  g_pot.smooth_pot = NULL;
  g_pot.cp_start = 0;
  g_pot.global_idx = 0;
  g_pot.global_pot = 0;
  g_pot.have_globals = 0;
  g_pot.calc_list = NULL;
  g_pot.compnodelist = NULL;
  memset(&g_pot.apot_table, 0, sizeof(g_pot.apot_table));
#endif  // APOT

#if defined(PAIR)
  init_interaction_name("PAIR");
#elif defined(EAM) && !defined(COULOMB)
#if !defined(TBEAM)
  init_interaction_name("EAM");
#else
  init_interaction_name("TBEAM");
#endif /* TBEAM */
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
#endif /* TERSOFFMOD */
#endif /* interaction type */

  g_memory.pointers = NULL;
  g_memory.num_pointers = 0;
}

/****************************************************************
 *
 *  init_interaction_name
 *
 ****************************************************************/

void init_interaction_name(const char* name)
{
  int len = strlen(name);
  g_todo.interaction_name = (char*)Malloc((len + 1) * sizeof(char));
  strncpy(g_todo.interaction_name, name, len);
  g_todo.interaction_name[len] = '\0';
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
  if (pchar != NULL)
  {
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
  if (pdouble != NULL)
  {
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
  if (pint != NULL)
  {
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
  if (patom != NULL && natoms > 0)
  {
    for (int i = 0; i < natoms; i++)
    {
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

#if defined(APOT)

/****************************************************************
 *
 *  free_apot_table
 *
 ****************************************************************/

void free_apot_table(apot_table_t* papot)
{
  free_void_pointer(papot->idxpot);
  free_char_pointer2(papot->names, papot->number);
  free_void_pointer(papot->begin);
  free_void_pointer(papot->end);
  free_void_pointer(papot->n_par);
  free_void_pointer(papot->idxparam);
  free_int_pointer2(papot->invar_par, papot->number);
  //   free_char_pointer3(papot->param_name, papot->number);
  free_double_pointer2(papot->pmin, papot->number);
  free_double_pointer2(papot->values, papot->number);
  free_double_pointer2(papot->pmax, papot->number);
  free_void_pointer(papot->n_glob);
//   free_int_pointer3(papot->global_idx, ???);
#if defined(PAIR)
//   free_double_pointer2(papot->chempot, g_param.ntypes);
#endif  // PAIR
}

#endif // APOT
