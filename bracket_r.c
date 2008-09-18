/****************************************************************
* 
*  bracket_r.c: Brackets a minimum of a function.
* 
*****************************************************************/
/* Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi 
*                 (gsl/min/bracketing.c)
*            2005-2007 Peter Brommer
*            Institute for Theoretical and Applied Physics
*            University of Stuttgart, D-70550 Stuttgart, Germany
*            http://www.itap.physik.uni-stuttgart.de/
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
* $Revision: 1.4 $
* $Date: 2008/09/18 14:34:09 $
*****************************************************************/

#include <math.h>
#include "potfit.h"
#include "utils.h"

#define CGOLD 0.3819660
#define MAX_IT 100

#define P_SWAP(A,B,C) (C)=(A);(A)=(B);(B)=(C);


extern real *xicom, *delcom;

void bracket_r(real *x_lower, real *x_minimum, real *x_upper,
	       real *f_lower, real *f_minimum, real *f_upper,
	       real *f_vec1, real *f_vec2)
{
  /* The three following variables must be declared volatile to avoid storage
     in extended precision registers available on some architecture. The code
     relies on the ability to compare double values. As the values will be
     store in regular memory, the extended precision will then be lost and
     values that are different in extended precision might have equal
     representation in double precision. This behavior might break the
     algorithm. 
   */
  volatile double f_left = *f_lower;
  volatile double f_right = *f_upper;
  volatile double f_center;
  double x_left = *x_lower;
  double x_right = *x_upper;
  double x_center;
  static char errmsg[255];
  static real *vecu = NULL;	/* Vector of location u */
  static real *f_vec3 = NULL;	/* 3rd target vector */
  static real *p_left, *p_right, *p_center, *p_temp;
  int   j;
  int   last = 0;		/* indicates whether upwards is left or right */
  long  nb_eval = 0;
  if (vecu == NULL)
    vecu = vect_real(ndimtot);
  if (f_vec3 == NULL)
    f_vec3 = vect_real(mdim);

  p_left = f_vec1;
  p_right = f_vec2;
  p_center = f_vec3;

  if (f_right >= f_left) {
    x_center = x_left;
    f_center = f_left;
    P_SWAP(p_center, p_left, p_temp);
    x_left = -(x_right - x_center) / CGOLD + x_right;
    nb_eval++;
    for (j = 0; j < ndimtot; j++)
      vecu[j] = xicom[j] + x_left * delcom[j];	/*set vecu */
    f_left = (*calc_forces) (vecu, p_left, 0);
  } else {
    x_center = x_right;
    f_center = f_right;
    P_SWAP(p_center, p_right, p_temp);
    x_right = (x_center - x_left) / CGOLD + x_left;
    nb_eval++;
    for (j = 0; j < ndimtot; j++)
      vecu[j] = xicom[j] + x_right * delcom[j];	/*set vecu */
    f_right = (*calc_forces) (vecu, p_right, 0);
  }

  do {
    if (f_center < f_left) {
      if (f_center < f_right) {
	/* SUCCESS: Minimum in bracket! */
	*x_lower = x_left;
	*x_upper = x_right;
	*x_minimum = x_center;
	*f_lower = f_left;
	*f_upper = f_right;
	*f_minimum = f_center;
	for (j = 0; j < mdim; j++)
	  f_vec1[j] = p_center[j];
	return;
      } else if (f_center > f_right) {
	/* OK, go right! */
	x_left = x_center;
	f_left = f_center;
	P_SWAP(p_left, p_center, p_temp);
	x_center = x_right;
	f_center = f_right;
	P_SWAP(p_center, p_right, p_temp);
	x_right = (x_center - x_left) / CGOLD + x_left;
	nb_eval++;
	for (j = 0; j < ndimtot; j++)
	  vecu[j] = xicom[j] + x_right * delcom[j];
	f_right = (*calc_forces) (vecu, p_right, 0);
      } else {			/* f_center == f_right */

	/* Pathological: Search between center and right */
	/* This means a change from original algorithm */
#ifdef DEBUG
	sprintf(errmsg, "Pathological  @%li %f %f %f! center-right!\n",
		nb_eval, x_left, x_center, x_right);
	warning(errmsg);
#endif /* DEBUG */
	x_left = x_center;
	f_left = f_center;
	P_SWAP(p_left, p_center, p_temp);
	x_center = -(x_right - x_left) * CGOLD + x_right;
	nb_eval++;
	for (j = 0; j < ndimtot; j++)
	  vecu[j] = xicom[j] + x_center * delcom[j];
	f_center = (*calc_forces) (vecu, p_center, 0);
	last = 1;
      }
    } else if (f_center > f_left)
      /* Search to the left */
    {

      x_right = x_center;
      f_right = f_center;
      P_SWAP(p_right, p_center, p_temp);
      x_center = x_left;
      f_center = f_left;
      P_SWAP(p_center, p_left, p_temp);
      x_left = -(x_right - x_center) / CGOLD + x_right;
      nb_eval++;
      for (j = 0; j < ndimtot; j++)
	vecu[j] = xicom[j] + x_left * delcom[j];
      f_left = (*calc_forces) (vecu, p_left, 0);
    } else {			/* f_center == f_left */

      if (f_center < f_right) {
	/* between center and left */
#ifdef DEBUG
	sprintf(errmsg, "Pathological  @%li %f %f %f! center-left!\n",
		nb_eval, x_left, x_center, x_right);
	warning(errmsg);
#endif /* DEBUG */
	x_right = x_center;
	f_right = f_center;
	P_SWAP(p_right, p_center, p_temp);
	x_center = (x_right - x_left) * CGOLD + x_left;
	nb_eval++;
	for (j = 0; j < ndimtot; j++)
	  vecu[j] = xicom[j] + x_center * delcom[j];
	f_center = (*calc_forces) (vecu, p_center, 0);
	last = 2;
      } else if (f_center > f_right) {
	/* Search to the right */
	x_left = x_center;
	f_left = f_center;
	P_SWAP(p_left, p_center, p_temp);
	x_center = x_right;
	f_center = f_right;
	P_SWAP(p_center, p_right, p_temp);
	x_right = (x_center - x_left) / CGOLD + x_left;
	nb_eval++;
	for (j = 0; j < ndimtot; j++)
	  vecu[j] = xicom[j] + x_right * delcom[j];
	f_right = (*calc_forces) (vecu, p_right, 0);
      } else {			/* f_center==f_left==f_right */

	/* Kind of pathological case: go left/right in turns */
	if (last == 2) {
	  /* go further to left, it goes up towards the right */
#ifdef DEBUG

	  sprintf(errmsg, "Pathological  @%li %f %f %f! Go left!\n",
		  nb_eval, x_left, x_center, x_right);
	  warning(errmsg);
#endif /* DEBUG */
	  x_center = x_left;
	  f_center = f_left;
	  P_SWAP(p_center, p_left, p_temp);
	  x_left = -(x_right - x_center) / CGOLD + x_right;
	  nb_eval++;
	  for (j = 0; j < ndimtot; j++)
	    vecu[j] = xicom[j] + x_left * delcom[j];
	  f_left = (*calc_forces) (vecu, p_left, 0);
	  last = 1;
	} else {		/* go further to the right, to left it went up */

#ifdef DEBUG
	  sprintf(errmsg, "Pathological @%li %f %f %f! go right!\n",
		  nb_eval, x_left, x_center, x_right);
	  warning(errmsg);
#endif /* DEBUG */
	  x_center = x_right;
	  f_center = f_right;
	  P_SWAP(p_center, p_right, p_temp);
	  x_right = (x_center - x_left) / CGOLD + x_left;
	  nb_eval++;
	  for (j = 0; j < ndimtot; j++)
	    vecu[j] = xicom[j] + x_right * delcom[j];
	  f_right = (*calc_forces) (vecu, p_right, 0);
	  last = 2;
	}
      }
    }
  }
  while (nb_eval < MAX_IT);
#ifdef DEBUG
  sprintf(errmsg,
	  "Problems with bracketing minimum in %li tries: F(%.16g)=%.16g, F(%.16g)=%.16g, F(%.16g)=%.16g.",
	  nb_eval, x_left, f_left, x_center, f_center, x_right, f_right);
  error(errmsg);
#else /* DEBUG */
  error("Problems with bracketing of minimum, aborting");
#endif /* DEBUG */
  return;
}
