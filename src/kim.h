/****************************************************************
 *
 * kim.h: header file for KIM interface
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

#ifndef KIM_H_INCLUDED
#define KIM_H_INCLUDED

#if defined(KIM)

#include "KIM_API_C.h"
#include "KIM_API_status.h"

#include "kim.h"

#if !defined(DIM)
#define DIM 3
#else
STATIC_ASSERT(DIM==3, potfit_kim_support_requires_DIM_eq_3);
#endif

// called from potfit.c
void initialize_KIM();
void shutdown_KIM();

#endif // KIM

#endif // KIM_H_INCLUDED
