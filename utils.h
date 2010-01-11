/****************************************************************
*
*  utils.h: potfit utilities header file
*
*****************************************************************/
/*
*   Copyright 2002-2007 Peter Brommer
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
* $Revision: 1.6 $
* $Date: 2010/01/11 09:03:08 $
*****************************************************************/

#define ANSI

real  sqrreal(real x);
int  *vect_int(long dim);
real *vect_real(long dim);
real **mat_real(long rowdim, long coldim);
void  free_vect_real(real *vect);
void  free_vect_int(int *vect);
void  free_mat_real(real **matrix);
void  reg_for_free(void *p, char *name);
void  free_all_pointers();
