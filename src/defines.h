/****************************************************************
 *
 * defines.h: potfit header file
 *
 ****************************************************************
 *
 * Copyright 2002-2015
 *      Institute for Theoretical and Applied Physics
 *      University of Stuttgart, D-70550 Stuttgart, Germany
 *      http://potfit.sourceforge.net/
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

#ifndef POTFIT_DEFINES_H
#define POTFIT_DEFINES_H

/* general flag for threebody potentials (MEAM, Tersoff, SW, ...) */
#if defined(MEAM) || defined(STIWEB) || defined(TERSOFF)
#define THREEBODY
#endif /* MEAM || TERSOFF || STIWEB */

/* define EAM if TBEAM is defined */
#if defined(TBEAM) && !defined(EAM)
#define EAM
#endif /* TBEAM && !EAM */

#if defined(APOT)
#define APOT_STEPS 500          /* number of sampling points for analytic pot */
#define APOT_PUNISH 10e6        /* general value for apot punishments */
#endif /* APOT */

#if defined(EAM) || defined(ADP) || defined(MEAM)
#define DUMMY_WEIGHT 100.0
#endif /* EAM || ADP || MEAM */

#define FORCE_EPS 0.1

/****************************************************************
 *
 *  SLOTS: number of different distance tables used in the force calculations
 *
 *  In potfit all potentials are calculated via spline interpolations of pre-
 *  calculated potential tables. To speed up the spline calculation,
 *  the exact position of a neighbor distance has to be known with respect to the
 *  tabulated values. For each potential function in a force routine there should
 *  be a different slot with the corresponding potential table information.
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
#endif /* EAM || STIWEB */

#if defined(TBEAM) || defined(MEAM)
#undef SLOTS
#define SLOTS 3
#endif /* TBEAM || MEAM */

#if defined(ADP)
#undef SLOTS
#define SLOTS 4
#endif /* ADP */

#endif // POTFIT_DEFINES_H
