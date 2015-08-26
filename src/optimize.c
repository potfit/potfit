/****************************************************************
 *
 * optimize.c: Contains main potfit program
 *
 ****************************************************************
 *
 * Copyright 2002-2014
 *      Institute for Theoretical and Applied Physics
 *      University of Stuttgart, D-70550 Stuttgart, Germany
 *      http://potfit.sourceforge.net/
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

#include "potfit.h"

#include "optimize.h"

// individual optimization algorithms
void run_simulated_annealing(double* const xi);
void run_differential_evolution(double* const xi);
void run_powell_lsq(double* const xi);

void run_optimization()
{
  printf("\nStarting optimization with %d parameters.\n", g_calc.ndim);
  fflush(stdout);

  double* xi = g_pot.opt_pot.table;

#if !defined(EVO)
  run_simulated_annealing(xi);
#else  /* EVO */
  run_differential_evolution(xi);
#endif /* EVO */

  printf("\nStarting powell minimization ...\n");

  run_powell_lsq(xi);

  printf("\nFinished powell minimization, calculating errors ...\n");
}
