/****************************************************************
 *
 * bracket.h: header file for bracketing and brent minimization
 *
 ****************************************************************
 *
 * Copyright 2002-2013
 * 	Institute for Theoretical and Applied Physics
 * 	University of Stuttgart, D-70550 Stuttgart, Germany
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
 ****************************************************************/

#ifndef BRACKET_H
#define BRACKET_H

/* for bracket.c */
#define CGOLD 0.3819660
#define MAX_IT 100

/* for brent.c */
#define ITMAX 100
#define ZEPS 1.0e-9
#define SHIFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

extern double *xicom, *delcom;

void  bracket(double *, double *, double *, double *, double *, double *, double *, double *);
double brent(double, double, double, double, double, double *, double *, double *, double *);

double linmin(double *, double *, double, double *, double *, double *, double *);

#endif /* BRACKET_H */
