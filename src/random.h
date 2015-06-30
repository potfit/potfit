/****************************************************************
 *
 * random.h: Contains the dSFMT, a prng
 *      adapted from Mutsuo Saito and Makoto Matsumoto
 *
 *****************************************************************
 *
 * Copyright 2002-2014
 *      Institute for Theoretical and Applied Physics
 *      University of Stuttgart, D-70550 Stuttgart, Germany
 *      http://potfit.sourceforge.net/
 *
 *****************************************************************
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

#ifndef POTFIT_RANDOM_H
#define POTFIT_RANDOM_H

void init_rng(int seed);
double eqdist();
double normdist();

#endif  // POTFIT_RANDOM_H
