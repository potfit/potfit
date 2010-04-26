/****************************************************************
*
*  powell_lsq.h: Header file used by powell routines.
*
*****************************************************************/
/*
*   Copyright 2002-2005 Peter Brommer
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
*
*****************************************************************/

#ifdef POWELL
#define EXPOW			/* define Variables in powell */
#define INPOW(data) =data	/* initialize data only in powell */
#else
#define EXPOW extern		/* declare them extern otherwise */
#define INPOW(data)		/* skip initialization otherwise */
#endif

void  bracket_r(real *x_lower, real *x_minimum, real *x_upper, real *f_lower,
		real *f_minimum, real *f_upper, real *f_vec1, real *f_vec2);

real  brent_r(real ax, real bx, real cx, real fbx, real tol,
	      real *xmin, real *xmin2, real *fxmin, real *fxmin2);

void  copy_matrix(real **a, real **b, int n, int m);

void  copy_vector(real *a, real *b, int n);

int   gamma_init(real **gamma, real **d, real *xi, real *force_xi);

int   gamma_update(real **gamma, real a, real b, real *fa, real *fb,
		   real *delta, int j, int m, int n, real fmin);

void  lineqsys_init(real **gamma, real **lineqsys, real *deltaforce,
		    real *p, int n, int m);

void  lineqsys_update(real **gamma, real **lineqsys, real *force_xi,
		      real *p, int i, int n, int m);

real  linmin_r(real p[], real xi[], real fxi1, real *x1, real *x2,
	       real *fret1, real *fret2);

void  matdotvec(real **a, real *x, real *y, int n, int m);

real  normalize_vector(real *v, int n);
