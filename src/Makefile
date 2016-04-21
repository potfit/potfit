############################################################################
#
# potfit - open-source force-matching
#
# Copyright 2002-2016 - the potfit development team
#
# http://potfit.sourceforge.net/
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
BIN_DIR 	= ../bin/

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
CFLAGS := ${FLAGS}
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${DEBUG_FLAGS}
else
ifneq (,$(findstring asan,${MAKETARGET}))
CFLAGS += ${ASAN_FLAGS}
LFLAGS_SERIAL = ${ASAN_LFLAGS}
LFLAGS_MPI = ${ASAN_LFLAGS}
else
CFLAGS += ${OPT_FLAGS}
endif
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

POTFITHDR   	+= bracket.h
POTFITHDR 	+= config.h
POTFITHDR 	+= defines.h
POTFITHDR 	+= elements.h
POTFITHDR 	+= errors.h
POTFITHDR 	+= force.h
POTFITHDR 	+= memory.h
POTFITHDR 	+= mpi_utils.h
POTFITHDR 	+= optimize.h
POTFITHDR 	+= params.h
POTFITHDR 	+= potfit.h
POTFITHDR 	+= potential_input.h
POTFITHDR 	+= potential_output.h
POTFITHDR 	+= random.h
POTFITHDR 	+= splines.h
POTFITHDR 	+= types.h
POTFITHDR 	+= utils.h

POTFITSRC 	+= bracket.c
POTFITSRC 	+= brent.c
POTFITSRC 	+= config.c
POTFITSRC 	+= elements.c
POTFITSRC 	+= errors.c
POTFITSRC 	+= force_common.c
POTFITSRC 	+= linmin.c
POTFITSRC 	+= memory.c
POTFITSRC 	+= mpi_utils.c
POTFITSRC 	+= optimize.c
POTFITSRC	+= params.c
POTFITSRC 	+= potential_input.c
POTFITSRC 	+= potential_output.c
POTFITSRC 	+= potential_output_imd.c
POTFITSRC 	+= potential_output_lammps.c
POTFITSRC 	+= potfit.c
POTFITSRC 	+= powell_lsq.c
POTFITSRC 	+= random.c
POTFITSRC 	+= random_dsfmt.c
POTFITSRC 	+= simann.c
POTFITSRC 	+= splines.c
POTFITSRC 	+= utils.c

ifneq (,$(strip $(findstring kim,${MAKETARGET})))
  POTFITSRC      += kim.c
  POTFITSRC      += force_kim.c
  POTFITSRC      += potential_input_f5.c
  POTFITHDR      += kim.h
  ifeq (,$(strip $(findstring nolimits,${MAKETARGET})))
    POTFITHDR      += functions.h
    POTFITSRC      += functions.c
    POTFITSRC      += chempot.c
  endif
else # kim

  ifneq (,$(strip $(findstring apot,${MAKETARGET})))
    POTFITSRC 	+= potential_input_f0.c
  else
    POTFITSRC 	+= potential_input_f3.c
    POTFITSRC 	+= potential_input_f4.c
  endif

  ifneq (,$(strip $(findstring pair,${MAKETARGET})))
    POTFITSRC      += force_pair.c 
  endif

  ifneq (,$(strip $(findstring eam,${MAKETARGET})))
    ifneq (,$(strip $(findstring meam,${MAKETARGET})))
      POTFITSRC      += force_meam.c
    else ifneq (,$(strip $(findstring coulomb,${MAKETARGET})))
      POTFITSRC      += force_eam_elstat.c
    else ifneq (,$(strip $(findstring dipole,${MAKETARGET})))
      POTFITSRC      += force_eam_elstat.c
    else
      POTFITSRC      += force_eam.c 
    endif
  endif

  ifneq (,$(strip $(findstring coulomb,${MAKETARGET})))
    ifeq (,$(strip $(findstring eam,${MAKETARGET})))
      POTFITSRC      += force_elstat.c
    endif
  endif

  ifneq (,$(strip $(findstring dipole,${MAKETARGET})))
    ifeq (,$(strip $(findstring eam,${MAKETARGET})))
      POTFITSRC      += force_elstat.c
    endif
  endif

  ifneq (,$(strip $(findstring adp,${MAKETARGET})))
    POTFITSRC      += force_adp.c
  endif

  ifneq (,$(strip $(findstring stiweb,${MAKETARGET})))
    POTFITSRC      += force_stiweb.c
  endif

  ifneq (,$(strip $(findstring tersoff,${MAKETARGET})))
    POTFITSRC      += force_tersoff.c
  endif


  ifneq (,$(strip $(findstring apot,${MAKETARGET})))
    POTFITHDR      += functions.h
    POTFITSRC      += functions.c functions_impl.c
    ifneq (,$(strip $(findstring pair,${MAKETARGET})))
      POTFITHDR      += chempot.h
      POTFITSRC      += chempot.c
    endif
  else
    ifneq (,$(strip $(findstring meam,${MAKETARGET})))
      POTFITSRC      += rescale_meam.c
    else
      POTFITSRC      += rescale.c
    endif
  endif

endif #kim 
# added ends

ifneq (,$(strip $(findstring evo,${MAKETARGET})))
  POTFITSRC      += diff_evo.c
endif

ifneq (,$(strip $(findstring parab,${MAKETARGET})))
  POTFITSRC      += parabola.c
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
    ERROR += More than one potential model specified
  endif
  ifneq (,$(strip $(findstring meam,${MAKETARGET})))
    CFLAGS  += -DMEAM
  else ifneq (,$(strip $(findstring tbeam,${MAKETARGET})))
    CFLAGS  += -DEAM -DTBEAM
  else
    CFLAGS  += -DEAM
  endif
  INTERACTION = 1
endif

# COULOMB
ifneq (,$(strip $(findstring coulomb,${MAKETARGET})))
  ifeq (,$(strip $(findstring eam,${MAKETARGET})))
    ifneq (,$(findstring 1,${INTERACTION}))
      ERROR += More than one potential model specified
    endif
  endif
  ifeq (,$(strip $(findstring apot,${MAKETARGET})))
    ERROR += COULOMB does not support tabulated potentials
  endif
  CFLAGS  += -DCOULOMB
  INTERACTION = 1
endif

# DIPOLE
ifneq (,$(strip $(findstring dipole,${MAKETARGET})))
  ifeq (,$(strip $(findstring eam,${MAKETARGET})))
    ifneq (,$(findstring 1,${INTERACTION}))
      ERROR += More than one potential model specified
    endif
  endif
  ifeq (,$(strip $(findstring apot,${MAKETARGET})))
    ERROR += DIPOLE does not support tabulated potentials
  endif
  CFLAGS  += -DCOULOMB -DDIPOLE
  INTERACTION = 1
endif

# angular dependent potentials (ADP)
ifneq (,$(strip $(findstring adp,${MAKETARGET})))
  ifneq (,$(findstring 1,${INTERACTION}))
    ERROR += More than one potential model specified
  endif
  CFLAGS  += -DADP
  INTERACTION = 1
endif

# Stillinger-Weber potentials (STIWEB)
ifneq (,$(strip $(findstring stiweb,${MAKETARGET})))
  ifneq (,$(findstring 1,${INTERACTION}))
    ERROR += More than one potential model specified
  endif
  ifeq (,$(findstring apot,${MAKETARGET}))
    ERROR += STIWEB does not work without the apot flag
  endif
  CFLAGS  += -DSTIWEB
  INTERACTION = 1
endif

# Tersoff potentials (TERSOFF and TERSOFFMOD)
ifneq (,$(strip $(findstring tersoff,${MAKETARGET})))
  ifneq (,$(findstring 1,${INTERACTION}))
    ERROR += "More than one potential model specified."
  endif
  ifeq (,$(findstring apot,${MAKETARGET}))
    ERROR += "TERSOFF does not work without the apot flag!"
  endif
  ifneq (,$(findstring tersoffmod,${MAKETARGET}))
    CFLAGS += -DTERSOFFMOD
  endif
  CFLAGS  += -DTERSOFF
  INTERACTION = 1
endif

# added
ifneq (,$(findstring kim,${MAKETARGET}))
  INTERACTION = 1
endif 
# added ends

ifneq (,$(findstring 0,${INTERACTION}))
  ERROR += "No interaction model specified."
endif

# EVO - for differential evolution
ifneq (,$(findstring evo,${MAKETARGET}))
CFLAGS += -DEVO
endif

# APOT - for analytic potentials
ifneq (,$(findstring apot,${MAKETARGET}))
  ifneq (,$(findstring resc,${MAKETARGET}))
    ERROR += "Analytic potentials are incompatible with the rescale option!"
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
ERROR += "parab is no longer supported, please remove it from your target -- "
endif

ifneq (,$(findstring wzero,${MAKETARGET}))
ERROR += "wzero is no longer supported, please remove it from your target -- "
endif

ifneq (,$(findstring dist,${MAKETARGET}))
ifeq (,$(findstring MPI,${PARALLEL}))
CFLAGS += -DPDIST
else
ERROR += "dist is not mpi parallelized -- "
endif
endif

ifneq (,$(findstring newscale,${MAKETARGET}))
ERROR += "newscale is no longer supported, please remove it from your target -- "
endif

ifneq (,$(findstring fweight,${MAKETARGET}))
CFLAGS += -DFWEIGHT
endif

ifneq (,$(findstring contrib,${MAKETARGET}))
CFLAGS += -DCONTRIB
endif

# force acml4 or acml5 over acml
ifneq (,$(findstring acml,${MAKETARGET}))
ifeq (,$(findstring acml4,${MAKETARGET}))
ifeq (,$(findstring acml5,${MAKETARGET}))
ERROR += The acml target is obsolete. Please use acml4 or acml5.
endif
endif
endif

ifneq (,$(findstring acml4,${MAKETARGET}))
CFLAGS += -DACML -DACML4
endif

ifneq (,$(findstring acml5,${MAKETARGET}))
CFLAGS += -DACML -DACML5
endif

ifneq (,$(findstring resc,${MAKETARGET}))
CFLAGS += -DRESCALE
endif


# added    to enable KIM and NOLIMITS in the preprocessor
ifneq (,$(findstring kim,${MAKETARGET}))
  CFLAGS += -DKIM
  ifneq (,$(findstring nolimits,${MAKETARGET}))
    CFLAGS += -DNOLIMITS -DEAM
  else
    CFLAGS += -DAPOT -DPAIR
  endif
endif

# compability with KIM 
ifneq (,$(findstring kim,${MAKETARGET}))
  ifneq (,$(findstring dist,${MAKETARGET}))
    ERROR += "KIM potential does not support dist option!"
  endif
  ifneq (,$(findstring resc,${MAKETARGET}))
    ERROR += "KIM potential does not support resc option!"
  endif
endif 

ifneq (,$(findstring kim,${MAKETARGET}))
  ifneq (,$(findstring adp,${MAKETARGET}))
    ERROR += "Need not specify adp to use KIM potential!"
  else ifneq (,$(findstring coulomb,${MAKETARGET}))
    ERROR += "Need not to specify coulomb to use KIM potential!"
  else ifneq (,$(findstring dipole,${MAKETARGET}))
    ERROR += "Need not to specify dipole to use KIM potential!"
  else ifneq (,$(findstring meam,${MAKETARGET}))
    ERROR += "Need not to specify meam to use KIM potential!"
  else ifneq (,$(findstring stiweb,${MAKETARGET}))
    ERROR += "Need not to specify stiweb to use KIM potential!"
  else ifneq (,$(findstring tbeam,${MAKETARGET}))
    ERROR += "Need not to specify tbeam to use KIM potential!"
  else ifneq (,$(findstring tersoff,${MAKETARGET}))
    ERROR += "Need not to specify tersoff to use KIM potential!"
  else ifneq (,$(findstring tersoffmod,${MAKETARGET}))
    ERROR += "Need not to specify tersoffmod to use KIM potential!"
  endif
endif 
#added ends


# Substitute .o for .c to get the names of the object files
OBJECTS := $(subst .c,.o,${SOURCES})

###########################################################################
#
#	 Rules
#
###########################################################################

# all objects depend on headers
${OBJECTS}: ${HEADERS}

# How to compile *.c files
# special rules for force computation
powell_lsq.o: powell_lsq.c
	@echo " [CC] powell_lsq.c"
	@${CC} ${CFLAGS} ${CINCLUDE} -c $< || { \
		echo -e "The following command failed with the above error:\n"; \
		echo -e ${CC} ${CFLAGS} ${CINCLUDE} -c $<"\n"; \
		exit 1; \
		}

# special rules for function evaluation
utils.o: utils.c
	@echo " [CC] utils.c"
	@${CC} ${CFLAGS} ${CINCLUDE} -c $< || { \
		echo -e "The following command failed with the above error:\n"; \
		echo -e ${CC} ${CFLAGS} ${CINCLUDE} -c $<"\n"; \
		exit 1; \
		}

# generic compilation rule
%.o: %.c
ifeq (,${MAKETARGET})
	@echo -e "Usage:"
	@echo -e "  make potfit_[interaction]_[options]\n"
	@echo "For more details on compiling potfit please look at the Makefile"
	@exit
else
	@echo " [CC] ${@:.o=.c}"
	@${CC} ${CFLAGS} -c $< || { \
		echo -e "The following command failed with the above error:\n"; \
		echo -e ${CC} ${CFLAGS} -c $<"\n"; \
		exit 1; \
		}
endif

# How to link
${MAKETARGET}: ${OBJECTS}
	@echo " [LD] $@"
ifeq (,${BIN_DIR})
	@${CC} ${LFLAGS_${PARALLEL}} -o $@ ${OBJECTS} ${LIBS}
else
	@${CC} ${LFLAGS_${PARALLEL}} -o ${BIN_DIR}/$@ ${OBJECTS} ${LIBS}
endif
ifneq (,${STRIP})
  ifeq (,$(findstring prof,${MAKETARGET}))
    ifeq (,$(findstring debug,${MAKETARGET}))
      ifeq (,${BIN_DIR})
	@${STRIP} --strip-unneeded -R .comment $@
      else
	@${STRIP} --strip-unneeded -R .comment ${BIN_DIR}/$@
      endif
    endif
  endif
endif
	@echo -e "Building $@ was sucessfull."

# First recursion only set the MAKETARGET Variable
.DEFAULT:
ifneq (,${CC})
	@${MAKE} --no-print-directory MAKETARGET='$@' STAGE2
else
	@echo "There is no compiler defined for this option."
	@echo -e "Please adjust the Makefile.\n"
	@exit
endif


# Second recursion sets MAKETARGET variable and compiles
# An empty MAKETARGET variable would create an infinite recursion, so we check
STAGE2:
ifneq (,${ERROR})
	@echo -e "\nError: ${ERROR}\n"
else
ifneq (,${MAKETARGET})
ifneq (,${WARNING})
	@echo "${WARNING}"
endif
	@${MAKE} MAKETARGET='${MAKETARGET}' ${MAKETARGET}
else
	@echo 'No TARGET specified.'
endif
endif

###########################################################################
#
#	 Misc. TARGETs
#
###########################################################################

clean:
	@rm -f *.o *.u *~ \#* *.V *.T *.O *.il

distclean:
	@rm -f *.o *.u *~ \#* *.V *.T *.O *.il

help:
	@echo "Usage: make potfit[_<parallel>][_<option>[_<option>...]]"

