/****************************************************************
 *
 * defines.h: potfit header file
 *
 ****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
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

#ifndef DEFINES_H_INCLUDED
#define DEFINES_H_INCLUDED

#if !defined(M_PI)
#define M_PI 3.14159265358979323846f
#endif  // M_PI

// general flag for threebody potentials (MEAM, Tersoff, SW, ...)
#if defined(MEAM) || defined(STIWEB) || defined(TERSOFF)
#define THREEBODY
#endif  // MEAM || TERSOFF || STIWEB

#if defined(KIM)
#define APOT
#define PAIR
#endif // KIM

#if defined(APOT)
#define APOT_STEPS 500    // number of sampling points for analytic pot
#define APOT_PUNISH 10e6  // general value for apot punishments
#endif                    // APOT

#if defined(EAM) || defined(ADP) || defined(MEAM)
#define DUMMY_WEIGHT 100.0
#endif  // EAM || ADP || MEAM

#define FORCE_EPS 0.1

#if defined(COULOMB) || defined(DIPOLE)
#define DP_EPS 14.40  // this is e^2/(4*pi*epsilon_0) in eV A
#endif                // COULOMB || DIPOLE

/****************************************************************
 *
 *  SLOTS: number of different distance tables used in the force calculations
 *
 *  In potfit all potentials are calculated via spline interpolations of pre-
 *  calculated potential tables. To speed up the spline calculation,
 *  the exact position of a neighbor distance has to be known with respect to
 *  the tabulated values. For each potential function in a force routine there
 *  should be a different slot with the corresponding potential table
 *information.
 *
 *  SLOTS = 1 for the following interactions:
 *      PAIR, COULOMB, DIPOLE, TERSOFF
 *      0 ... pair distance
 *
 *  EAM:        SLOTS = 2
 *      0 ... pair distance
 *      1 ... transfer function
 *
 *  STIWEB:     SLOTS = 2
 *      0 ... pair distance
 *      1 ... exponential functions
 *
 *  TBEAM:      SLOTS = 3
 *      0 ... pair distance
 *      1 ... transfer function
 *      2 ... 2nd transfer function
 *
 *  MEAM:       SLOTS = 3
 *      0 ... pair distance
 *      1 ... transfer function
 *      2 ... f(r_ij)
 *
 *  ADP:        SLOTS = 4
 *      0 ... pair distance
 *      1 ... transfer function
 *      2 ... dipole term
 *      3 ... quadrupole term
 *
 ****************************************************************/

#define SLOTS 1

#if defined(EAM) || defined(STIWEB)
#undef SLOTS
#define SLOTS 2
#endif  // EAM || STIWEB

#if defined(TBEAM) || defined(MEAM)
#undef SLOTS
#define SLOTS 3
#endif  // TBEAM || MEAM

#if defined(ADP)
#undef SLOTS
#define SLOTS 4
#endif  // ADP

#endif  // DEFINES_H_INCLUDED
