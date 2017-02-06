############################################################################
#
# potfit - open-source force-matching
#
# Copyright 2002-2017 - the potfit development team
#
# https://www.potfit.net/
#
############################################################################
#
#   This file is part of potfit.
#
#   potfit is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   potfit is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with potfit; if not, see <http://www.gnu.org/licenses/>.
#
############################################################################

include ../Makefile.inc

###########################################################################
#
#  Defaults for some variables
#
###########################################################################

STRIP 		= $(shell which strip 2> /dev/null)
LIBS		+= -lm
MPI_FLAGS	+= -DMPI
DEBUG_FLAGS	+= -DDEBUG
ACML4PATH 	= ${ACML4DIR}/lib
ACML5PATH 	= ${ACML5DIR}/lib
BIN_DIR 	= bin/

###########################################################################
#
#  Parallelization method
#
###########################################################################

# default is serial
PARALLEL = SERIAL
# MPI
ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
PARALLEL = MPI
endif

###########################################################################
#
#  Compiler, flags, libraries
#
###########################################################################

# compiler; if empty, we issue an error later
CC = ${CC_${PARALLEL}}

# optimization flags
OPT_FLAGS   += ${${PARALLEL}_FLAGS} ${OPT_${PARALLEL}_FLAGS} -DNDEBUG
DEBUG_FLAGS += ${${PARALLEL}_FLAGS} ${DEBUG_${PARALLEL}_FLAGS}

# libraries
LIBS += ${${PARALLEL}_LIBS}

# optimization or debug
CFLAGS = ${FLAGS}

ifneq (,$(findstring asan,${MAKETARGET}))
  CFLAGS += ${ASAN_FLAGS} ${DEBUG_FLAGS}
  LFLAGS_SERIAL = ${ASAN_LFLAGS}
  LFLAGS_MPI = ${ASAN_LFLAGS}
else ifneq (,$(findstring debug,${MAKETARGET}))
  CFLAGS += ${DEBUG_FLAGS}
else
  CFLAGS += ${OPT_FLAGS}
endif

# profiling support
ifneq (,$(findstring prof,${MAKETARGET}))
CFLAGS += ${PROF_FLAGS}
LIBS   += ${PROF_LIBS}
endif

###########################################################################
#
# potfit sources
#
###########################################################################

# header files needed for all options

POTFITHDR	+= bracket.h
POTFITHDR	+= config.h
POTFITHDR	+= defines.h
POTFITHDR	+= elements.h
POTFITHDR	+= errors.h
POTFITHDR	+= force.h
POTFITHDR	+= memory.h
POTFITHDR	+= mpi_utils.h
POTFITHDR	+= optimize.h
POTFITHDR	+= params.h
POTFITHDR	+= potfit.h
POTFITHDR	+= potential_input.h
POTFITHDR	+= potential_output.h
POTFITHDR	+= random.h
POTFITHDR	+= splines.h
POTFITHDR	+= types.h
POTFITHDR	+= utils.h

# source files needed for all options

POTFITSRC	+= bracket.c
POTFITSRC	+= brent.c
POTFITSRC	+= config.c
POTFITSRC	+= elements.c
POTFITSRC	+= errors.c
POTFITSRC	+= force_common.c
POTFITSRC	+= linmin.c
POTFITSRC	+= memory.c
POTFITSRC	+= mpi_utils.c
POTFITSRC	+= optimize.c
POTFITSRC	+= params.c
POTFITSRC	+= potential_input.c
POTFITSRC	+= potential_output.c
POTFITSRC	+= potential_output_imd.c
POTFITSRC	+= potential_output_lammps.c
POTFITSRC	+= potfit.c
POTFITSRC	+= powell_lsq.c
POTFITSRC	+= random.c
POTFITSRC	+= random_dsfmt.c
POTFITSRC	+= simann.c
POTFITSRC	+= splines.c
POTFITSRC	+= utils.c

# potential representations

POTFITSRC	+= potential_input_f0.c
POTFITSRC	+= potential_input_f3.c
POTFITSRC	+= potential_input_f4.c
POTFITSRC	+= potential_input_f5.c

# force routines

ifneq (,$(strip $(findstring pair,${MAKETARGET})))
  POTFITSRC	+= force_pair.c
endif

ifneq (,$(strip $(findstring eam,${MAKETARGET})))
  ifneq (,$(strip $(findstring meam,${MAKETARGET})))
    POTFITSRC	+= force_meam.c
  else ifneq (,$(strip $(findstring coulomb,${MAKETARGET})))
    POTFITSRC	+= force_eam_elstat.c
  else ifneq (,$(strip $(findstring dipole,${MAKETARGET})))
    POTFITSRC	+= force_eam_elstat.c
  else
    POTFITSRC	+= force_eam.c
  endif
endif

ifneq (,$(strip $(findstring coulomb,${MAKETARGET})))
  ifeq (,$(strip $(findstring eam,${MAKETARGET})))
    POTFITSRC	+= force_elstat.c
  endif
endif

ifneq (,$(strip $(findstring dipole,${MAKETARGET})))
  ifeq (,$(strip $(findstring eam,${MAKETARGET})))
    POTFITSRC	+= force_elstat.c
  endif
endif

ifneq (,$(strip $(findstring adp,${MAKETARGET})))
  POTFITSRC	+= force_adp.c
endif

ifneq (,$(strip $(findstring stiweb,${MAKETARGET})))
  POTFITSRC	+= force_stiweb.c
endif

ifneq (,$(strip $(findstring tersoff,${MAKETARGET})))
  POTFITSRC	+= force_tersoff.c
endif

# additional files for certain options

ifneq (,$(strip $(findstring apot,${MAKETARGET})))
  POTFITHDR	+= functions.h
  POTFITSRC	+= functions.c functions_impl.c
  ifneq (,$(strip $(findstring pair,${MAKETARGET})))
    POTFITHDR	+= chempot.h
    POTFITSRC	+= chempot.c
  endif
else ifneq (,$(strip $(findstring kim,${MAKETARGET})))
  POTFITHDR	+= chempot.h
  POTFITHDR	+= functions.h
  POTFITHDR	+= kim.h

  POTFITSRC	+= chempot.c
  POTFITSRC	+= force_kim.c
  POTFITSRC	+= functions.c
  POTFITSRC	+= functions_impl.c
  POTFITSRC	+= kim.c
else
  ifneq (,$(strip $(findstring meam,${MAKETARGET})))
    POTFITSRC	+= rescale_meam.c
  else
    POTFITSRC	+= rescale.c
  endif
endif

ifneq (,$(strip $(findstring evo,${MAKETARGET})))
  POTFITSRC	+= diff_evo.c
endif

ifneq (,$(strip $(findstring parab,${MAKETARGET})))
  POTFITSRC	+= parabola.c
endif

#########################################################
#
# potfit Configuration rules
#
#########################################################

HEADERS := ${POTFITHDR}
SOURCES	:= ${POTFITSRC}

###  INTERACTIONS  #######################################

INTERACTION = 0

# pair potentials
ifneq (,$(findstring pair,${MAKETARGET}))
  CFLAGS += -DPAIR
  INTERACTION = 1
endif

# embedded atom method (EAM) potentials
ifneq (,$(strip $(findstring eam,${MAKETARGET})))
  ifneq (,$(findstring 1,${INTERACTION}))
    ERROR += "More than one potential model specified!\n"
  endif
  ifneq (,$(strip $(findstring meam,${MAKETARGET})))
    CFLAGS += -DMEAM
  else ifneq (,$(strip $(findstring tbeam,${MAKETARGET})))
    CFLAGS += -DEAM -DTBEAM
  else
    CFLAGS += -DEAM
  endif
  INTERACTION = 1
endif

# COULOMB
ifneq (,$(strip $(findstring coulomb,${MAKETARGET})))
  ifeq (,$(strip $(findstring eam,${MAKETARGET})))
    ifneq (,$(findstring 1,${INTERACTION}))
      ERROR += "More than one potential model specified!\n"
    endif
  endif
  ifeq (,$(strip $(findstring apot,${MAKETARGET})))
    ERROR += "COULOMB does not support tabulated potentials!\n"
  endif
  CFLAGS += -DCOULOMB
  ifneq (,$(strip $(findstring dsf,${MAKETARGET})))
    CFLAGS += -DDSF
  endif
  INTERACTION = 1
endif

# DIPOLE
ifneq (,$(strip $(findstring dipole,${MAKETARGET})))
  ifeq (,$(strip $(findstring eam,${MAKETARGET})))
    ifneq (,$(findstring 1,${INTERACTION}))
      ERROR += "More than one potential model specified!\n"
    endif
  endif
  ifneq (,$(strip $(findstring dsf,${MAKETARGET})))
      ERROR += "DIPOLE and DSF are not currently compatible!\\n"
  endif
  ifeq (,$(strip $(findstring apot,${MAKETARGET})))
    ERROR += "DIPOLE does not support tabulated potentials!\n"
  endif
  CFLAGS += -DCOULOMB -DDIPOLE
  INTERACTION = 1
endif

# angular dependent potentials (ADP)
ifneq (,$(strip $(findstring adp,${MAKETARGET})))
  ifneq (,$(findstring 1,${INTERACTION}))
    ERROR += "More than one potential model specified!\n"
  endif
  CFLAGS += -DADP
  INTERACTION = 1
endif

# Stillinger-Weber potentials (STIWEB)
ifneq (,$(strip $(findstring stiweb,${MAKETARGET})))
  ifneq (,$(findstring 1,${INTERACTION}))
    ERROR += "More than one potential model specified!\n"
  endif
  ifeq (,$(findstring apot,${MAKETARGET}))
    ERROR += "STIWEB does not work without the apot flag!\n"
  endif
  CFLAGS += -DSTIWEB
  INTERACTION = 1
endif

# Tersoff potentials (TERSOFF and TERSOFFMOD)
ifneq (,$(strip $(findstring tersoff,${MAKETARGET})))
  ifneq (,$(findstring 1,${INTERACTION}))
    ERROR += "More than one potential model specified!\n"
  endif
  ifeq (,$(findstring apot,${MAKETARGET}))
    ERROR += "TERSOFF does not work without the apot flag!\n"
  endif
  ifneq (,$(findstring tersoffmod,${MAKETARGET}))
    CFLAGS += -DTERSOFFMOD
  endif
  CFLAGS += -DTERSOFF
  INTERACTION = 1
endif

ifneq (,$(findstring kim,${MAKETARGET}))
  ifneq (,$(findstring 1,${INTERACTION}))
    ERROR += "More than one potential model specified.\n"
  endif
  INTERACTION = 1
endif

ifneq (,$(findstring 0,${INTERACTION}))
  ERROR += "No interaction model specified!\n"
endif

# EVO - for differential evolution
ifneq (,$(findstring evo,${MAKETARGET}))
  CFLAGS += -DEVO
endif

# APOT - for analytic potentials
ifneq (,$(findstring apot,${MAKETARGET}))
  ifneq (,$(findstring resc,${MAKETARGET}))
    ERROR += "Analytic potentials are incompatible with the rescale option!\n"
  endif
  CFLAGS += -DAPOT
endif

# Stress
ifneq (,$(findstring stress,${MAKETARGET}))
  CFLAGS += -DSTRESS
endif

# Disable gauge punishments for EAM/ADP
ifneq (,$(findstring nopunish,${MAKETARGET}))
  CFLAGS += -DNOPUNISH
endif

ifneq (,$(findstring limit,${MAKETARGET}))
  WARNING += "limit is now mandatory -- "
endif

ifneq (,$(findstring parab,${MAKETARGET}))
  ERROR += "parab is no longer supported, please remove it from your target!\n"
endif

ifneq (,$(findstring wzero,${MAKETARGET}))
  ERROR += "wzero is no longer supported, please remove it from your target!\n"
endif

ifneq (,$(findstring dist,${MAKETARGET}))
  ifeq (,$(findstring MPI,${PARALLEL}))
    CFLAGS += -DPDIST
  else
    ERROR += "dist is not mpi parallelized!\n"
  endif
endif

ifneq (,$(findstring newscale,${MAKETARGET}))
  ERROR += "newscale is no longer supported, please remove it from your target!\n"
endif

ifneq (,$(findstring fweight,${MAKETARGET}))
  CFLAGS += -DFWEIGHT
endif

ifneq (,$(findstring contrib,${MAKETARGET}))
  CFLAGS += -DCONTRIB
endif

ifneq (,$(findstring mkl,${MAKETARGET}))
  ifneq (,$(findstring acml,${MAKETARGET}))
    ERROR += "ACML cannot be used together with MKL.\n"
  endif
endif

# force acml4 or acml5 over acml
ifneq (,$(findstring acml,${MAKETARGET}))
  ifneq (,$(findstring mkl,${MAKETARGET}))
    ERROR += "ACML cannot be used together with MKL.\n"
  endif
  ifeq (,$(findstring acml4,${MAKETARGET}))
    ifeq (,$(findstring acml5,${MAKETARGET}))
      ERROR += "The acml target is obsolete. Please use acml4 or acml5.\n"
    endif
  endif
  CFLAGS += -DACML
endif

ifneq (,$(findstring resc,${MAKETARGET}))
  CFLAGS += -DRESCALE
endif

ifneq (,$(findstring kim,${MAKETARGET}))
  CFLAGS += -DKIM
endif

# compability with KIM
ifneq ($(findstring kim,${MAKETARGET}),)
  ifeq ($(findstring kim-toolchain,${SYSTEM}),)
    ERROR += "KIM support requires the SYSTEM variable setting 'kim-toolchain'.\n"
  endif
endif
ifneq ($(findstring kim-toolchain,${SYSTEM}),)
  ifeq ($(findstring kim,${MAKETARGET}),)
    ERROR += "SYSTEM=kim-toolchain is not supported for non-KIM builds.\n"
  endif
endif

ifneq (,$(findstring kim,${MAKETARGET}))
  ifeq ($(shell kim-api-build-config --version 2> /dev/null),)
    ERROR += "kim-api-build-config could not be found. Check your PATH variable.\n"
  endif
  ifneq (,$(findstring apot,${MAKETARGET}))
    ERROR += "KIM option uses the OpenKIM internal data format, please remove apot from target!\n"
  endif
  ifneq (,$(findstring dist,${MAKETARGET}))
    ERROR += "KIM potential does not support dist option!\n"
  endif
  ifneq (,$(findstring resc,${MAKETARGET}))
    ERROR += "KIM potential does not support resc option!\n"
  endif
endif

CFLAGS += -D${TARGET} -D${MATH_LIB}

# Substitute .o for .c to get the names of the object files
OBJECTS := $(subst .c,.o,${SOURCES})

###########################################################################
#
# targets for building potfit
#
###########################################################################

# all objects depend on headers
${OBJECTS}: ${HEADERS}

# How to compile *.c files
%.o: %.c
	@echo " [CC] ${@:.o=.c}"
	@${CC} ${CFLAGS} ${CINCLUDE} -c $< || { \
		echo -e "The following command failed with the above error:\n"; \
		echo -e ${CC} ${CFLAGS} -c $<"\n"; \
		exit -1; \
		}

# check for errors here
.PHONY: errors
errors:
ifneq (,${ERROR})
	@echo -e "\nThe following error(s) occurred:"
	@echo -e ${ERROR}
	@exit -1
endif

# How to link
${MAKETARGET}: errors ${OBJECTS}
	@echo " [LD] $@"
	@${CC} ${LFLAGS_${PARALLEL}} -o $@ ${OBJECTS} ${LIBS}
ifneq (,${STRIP})
  ifeq (,$(findstring prof,${MAKETARGET}))
    ifeq (,$(findstring debug,${MAKETARGET}))
      ifeq (,$(findstring asan,${MAKETARGET}))
	@${STRIP} --strip-unneeded -R .comment $@
      endif
    endif
  endif
endif
ifneq (,${BIN_DIR})
  ifneq (,$(wildcard ../${BIN_DIR}))
	@mv $@ ../${BIN_DIR}
  endif
endif
	@echo -e "\nCompiled with options"
	@echo -e "CC=${CC} \nCFLAGS=${CFLAGS} \nLIBS=${LIBS}\n"
	@echo -e "Building $@ was sucessfull."

###########################################################################
#
# misc. targets
#
###########################################################################

clean:
	@rm -f *.o *.u *~ \#* *.V *.T *.O *.il

distclean:
	@rm -f *.o *.u *~ \#* *.V *.T *.O *.il
