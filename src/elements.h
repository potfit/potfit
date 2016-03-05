/****************************************************************
 *
 * elements.h: header file for periodic table
 *
 ****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
 *
 * http://potfit.sourceforge.net/
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

#ifndef ELEMENTS_H_INCLUDED
#define ELEMENTS_H_INCLUDED

// get mass from number
double ele_mass_from_number(int num);

// get mass from name
double ele_mass_from_name(const char* name);

// get number from name
int ele_number_from_name(const char* name);

#endif  // ELEMENTS_H_INCLUDED
