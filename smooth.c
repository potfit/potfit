/****************************************************************
*
*  smooth.c: Calculates a cutoff radius for a smooth cutoff
*
*****************************************************************/
/*
*   Copyright 2008-2009 Daniel Schopf
*             Institute for Theoretical and Applied Physics
*             University of Stuttgart, D-70550 Stuttgart, Germany
*             http://www.itap.physik.uni-stuttgart.de/
*
*****************************************************************/
/*
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
*   along with potfit; if not, write to the Free Software
*   Foundation, Inc., 51 Franklin St, Fifth Floor,
*   Boston, MA  02110-1301  USA
*/
/****************************************************************
* $Revision: 1.9 $
* $Date: 2010/01/11 09:03:08 $
*****************************************************************/

#ifdef APOT

#include "potfit.h"

/******************************************************************************
*
* function for smooth cutoff radius
*
******************************************************************************/

real cutoff(real r, real r0, real h)
{
  real  val = 0;

  if ((r - r0) > 0)
    return 0;
  val = (r - r0) / h;
  val *= val;
  val *= val;

  return val / (1 + val);
}

#endif /* APOT */
