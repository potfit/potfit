/****************************************************************
*
*  erfc.h: approximation of error function for wolf summation
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
* $Revision: 1.1 $
* $Date: 2010/04/14 15:01:22 $
*****************************************************************/


real erfc(real x)
{

#define A_1 0.254829592
#define A_2 -0.284496736
#define A_3 1.421413741
#define A_4 -1.453152027
#define A_5 1.061405429
#define P   0.3275911

  real t, xsq, tp;

  t = 1.0 / ( 1.0 + P * x );

  xsq = x * x;

  tp = t * ( A_1 + t * ( A_2 + t * ( A_3 + t * ( A_4 + t * A_5 ) ) ) );
  
  return ( tp * exp( -xsq ) );

} 
