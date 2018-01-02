/****************************************************************
 *
 * config.h: header file for reading atomic configurations
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
 ****************************************************************/

#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

void read_config(const char* filename);

#if defined(APOT)
void update_slots(void);
void update_neighbor_slots(neigh_t* neighbor, double r, int neighbor_slot);
#endif  // APOT

#endif  // CONFIG_H_INCLUDED
