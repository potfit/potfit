############################################################################
#
# potfit -- The ITAP Force Matching Program
# 	Copyright 2002-2010
#
# 	Institute for Theoretical and Applied Physics,
# 	University of Stuttgart, D-70550 Stuttgart, Germany
# 	http://www.itap.physik.uni-stuttgart.de/
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
#
# Beware: This Makefile works only with GNU make (gmake)!
#
# Usage:  make <target>
#
# <target> has the form
#
#    potfit[_<parallel>][_<option>[_<option>...]]
#
# The parallelization method <parallel> can be:
#
#    mpi   compile for parallel execution, using MPI
#
###########################################################################
#
# Customizing this Makefile
#
# As potfit supports a large number of compile options, you will have to
# compile potfit freqently. Before doing so, however, you must check whether
# the settings in this Makefile fit your needs. You possibly have to
# customize these setttings. Before you can do that, we have to explain
# a bit how the compilation process works.
#
# The compilation process requires the SYSTEM variable in the Makefile to be
# set to any of the predefined values. It specifies what system you have, and
# what compiler you are using. The flags for the compiler and the linker
# are then selected as a function of this variable.
#
# Another important ingredient is the parallelization method, which is
# determined from the make target. The parallelization method is stored
# in the variable PARALLEL, which takes as value SERIAL or MPI.
#
# Depending on the value of ${SYSTEM}, a number of variables must be
# set, from which everything else is constructed.
#
# CC_SERIAL defines the compiler for serial compilation, CC_MPI the one
# to be used for parallelization
#
# BIN_DIR defines the directory where the potfit binary is put. Note that
# this directory must exist.
#
# MV defines the program used to move the potfit binary to ${BIN_DIR}.
# The default is mv, which is usually ok.
#
# The compilation options are stored in the variable CFLAGS.
# The initial value of CFLAGS is set to the variable FLAGS,
# which can be given on the command line.
#
# If the option debug was specified, ${DEBUG_FLAGS} is then appended
# to ${CFLAGS}, otherwise ${OPT_FLAGS}. If the option prof was specified
# (for profiling), ${PROF_FLAGS} is also appended to ${CFLAGS}. However,
# before appending ${OPT_FLAGS} or ${DEBUG_FLAGS} to ${CFLAGS}, some
# parallelization specific flags are appended to them:
#
#   OPT_FLAGS   += ${${PARALLEL}_FLAGS} ${OPT_${PARALLEL}_FLAGS}
#   DEBUG_FLAGS += ${${PARALLEL}_FLAGS} ${DEBUG_${PARALLEL}_FLAGS}
#
# If any of these variables is not defined, it is assumed to be empty.
# This setup should provide sufficient flexibility to set one's favorite
# flags, depending on parallelization, profiling, and optimization/debugging.
#
# Similarly, the link libraries are stored in the variable LIBS,
# to which ${${PARALLEL}_LIBS} and possibly ${PROF_LIBS} (for profiling)
# is appended.
#
# You may have to change the setting for an existing value of SYSTEM.
# or you have to add support for a new value of SYSTEM. The latter is
# best done by using the folloing template for SYSTEM=custom:
#
# ifeq (custom,${SYSTEM})
#   CC_SERIAL		= serial-compiler
#   CC_MPI		= MPI-compiler
#   OMPI_CC      	= compiler for mpicc
#   OMPI_CLINKER 	= linker for mpicc
#   OPT_FLAGS		+= generic flags for optimization
#   DEBUG_FLAGS		+= generic flags for debugging
#   PROF_FLAGS		+= flags for profiling
#   PROF_LIBS		+= libraries for profiling
#   LFLAGS_SERIAL 	+= flags for serial linking
#   LFLAGS_MPI 		+= flags for MPI linking
#   export        MPICH_CC MPICH_CLINKER
# endif
#
# Variables remaining empty need not be mentioned.

###########################################################################
#
#  Adjust these variables to your system
#
###########################################################################

# Currently the following systems are available:
# x86_64-icc  	64bit Intel Compiler
# x86_64-gcc    64bit GNU Compiler
# i586-icc 	32bit Intel Compiler
# i586-gcc  	32bit GNU Compiler
SYSTEM 		= x86_64-icc

# This is the directory where the potfit binary will be moved to
BIN_DIR 	= ${HOME}/bin

# Base directory of your installation of the MKL or ACML
MKLDIR          = /common/linux/paket/intel/compiler-11.0/cc/mkl
ACMLDIR  	= /common/linux/paket/acml4.4.0/ifort64
#ACMLDIR  	= /opt/acml4.4.0/gfortran64

###########################################################################
#
#  Defaults for some variables
#
###########################################################################

MV		= $(shell basename `which mv`)
STRIP 		= $(shell basename `which strip`)
LIBS		+= -lm
MPI_FLAGS	+= -DMPI
DEBUG_FLAGS	+= -DDEBUG
MKLPATH         = ${MKLDIR}/lib
RELEASE		= 0

###########################################################################
#
#  flags for 64bit
#
###########################################################################

ifeq (x86_64-icc,${SYSTEM})
# compiler
  CC_SERIAL     = icc
  CC_MPI        = mpicc

# general optimization flags
  OPT_FLAGS     += -fast -xHost -multiple-processes

  OMPI_CC      = icc
  OMPI_CLINKER = icc

# debug flags
  PROF_FLAGS    += -prof-gen
  PROF_LIBS 	+= -prof-gen
  DEBUG_FLAGS   += -g -Wall -wd981 -wd1572

# Intel Math Kernel Library
ifeq (,$(strip $(findstring acml,${MAKETARGET})))
  MKLPATH       = ${MKLDIR}/lib/em64t
  CINCLUDE 	+= -I${MKLDIR}/include
  LIBS 		+= -L${MKLPATH} ${MKLPATH}/libmkl_solver_lp64_sequential.a \
		   -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential \
		   -lmkl_core -Wl,--end-group -lpthread
endif

# AMD Core Math Library
ifneq (,$(strip $(findstring acml,${MAKETARGET})))
   ACMLPATH 	= ${ACMLDIR}/lib
   CINCLUDE     += -I$(ACMLDIR)/include
   LD_LIBRARY_PATH +=':$(ACMLPATH):'
   LIBS		:= -L${ACMLPATH} -lpthread -lacml -lacml_mv
endif

 export        OMPI_CC OMPI_CLINKER
endif

ifeq (x86_64-gcc,${SYSTEM})
# compiler
  CC_SERIAL     = gcc
  CC_MPI        = mpicc

# general optimization flags
  OPT_FLAGS     += -O3 -march=native -pipe -Wno-unused

  OMPI_CC      	= gcc
  OMPI_CLINKER 	= gcc

# debug flags
  PROF_FLAGS    += -pg
  DEBUG_FLAGS   += -g -Wall # -wd981 -wd1572

# Intel Math Kernel Library
ifeq (,$(strip $(findstring acml,${MAKETARGET})))
  MKLPATH       = ${MKLDIR}/lib/em64t/
  CINCLUDE      = -I${MKLDIR}/include
  LIBS 		+= -L${MKLPATH} ${MKLPATH}/libmkl_solver_lp64_sequential.a \
		   -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
		   -Wl,--end-group -lpthread -Wl,--as-needed
endif

# AMD Core Math Library
ifneq (,$(strip $(findstring acml,${MAKETARGET})))
  ACMLPATH      = ${ACMLDIR}/lib
  CINCLUDE     	+= -I$(ACMLDIR)/include
  LD_LIBRARY_PATH +=':$(ACMLPATH):'
  LIBS		+= -L${ACMLPATH} -lpthread -lacml -lacml_mv -Wl,--as-needed
endif

 export        OMPI_CC OMPI_CLINKER
endif


###########################################################################
#
#  flags for 32bit
#
###########################################################################

ifeq (i586-icc,${SYSTEM})
  CC_SERIAL	= icc
  CC_MPI	= mpicc
  OMPI_CC       = icc
  OMPI_CLINKER  = icc
  OPT_FLAGS	+= -fast -xHost -multiple-processes
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -prof-gen
  PROF_LIBS 	+= -prof-gen

# Intel Math Kernel Library
ifeq (,$(strip $(findstring acml,${MAKETARGET})))
  MKLPATH       = ${MKLDIR}/lib/32
  CINCLUDE      = -I${MKLDIR}/include
  LIBS 		+= -L${MKLPATH} ${MKLPATH}/libmkl_solver_sequential.a \
		   -Wl,--start-group -lmkl_intel -lmkl_sequential -lmkl_core \
		   -Wl,--end-group -lpthread
endif

# AMD Core Math Library
ifneq (,$(strip $(findstring acml,${MAKETARGET})))
  ACMLPATH      = ${ACMLDIR}/lib
  CINCLUDE     	+= -I$(ACMLDIR)/include
  LD_LIBRARY_PATH +=':$(ACMLPATH):'
  LIBS		+= -L${ACMLPATH} -lpthread -lacml -lacml_mv
endif

  export        OMPI_CC OMPI_CLINKER
endif

ifeq (i586-gcc,${SYSTEM})
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  OMPI_CC     	= gcc
  OMPI_CLINKER 	= gcc
  OPT_FLAGS	+= -O3 -march=native
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg

# Intel Math Kernel Library
ifeq (,$(strip $(findstring acml,${MAKETARGET})))
  MKLPATH       = ${MKLDIR}/lib/32
  CINCLUDE      = -I${MKLDIR}/include
  LIBS		+= -L${MKLPATH} ${MKLPATH}/libmkl_solver_sequential.a \
		   -Wl,--start-group -lmkl_intel -lmkl_sequential -lmkl_core \
		   -Wl,--end-group -lpthread -Wl,--as-needed
endif

# AMD Core Math Library
ifneq (,$(strip $(findstring acml,${MAKETARGET})))
  ACMLPATH      = ${ACMLDIR}/lib
  CINCLUDE     	+= -I$(ACMLDIR)/include
  LD_LIBRARY_PATH +=':$(ACMLPATH):'
  LIBS		+= -L${ACMLPATH} -lpthread -lacml -lacml_mv -Wl,--as-needed 
endif

  export        OMPI_CC OMPI_CLINKER
endif


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

POTFITHDR   	= bracket.h  potfit.h \
		  random-params.h  random.h  utils.h
POTFITSRC 	= bracket.c brent.c config.c linmin.c \
		  param.c potential.c potfit.c powell_lsq.c \
		  random.c simann.c splines.c utils.c

ifneq (,$(strip $(findstring pair,${MAKETARGET})))
POTFITSRC      += force_pair.c
endif

ifneq (,$(strip $(findstring eam,${MAKETARGET})))
POTFITSRC      += force_eam.c rescale.c
endif

ifneq (,$(strip $(findstring adp,${MAKETARGET})))
POTFITSRC      += force_adp.c rescale.c
endif

ifneq (,$(strip $(findstring apot,${MAKETARGET})))
POTFITSRC      += functions.c
ifneq (,$(strip $(findstring pair,${MAKETARGET})))
POTFITSRC      += chempot.c
endif
endif

ifneq (,$(strip $(findstring evo,${MAKETARGET})))
POTFITSRC      += diff_evo.c
endif

MPISRC          = mpi_utils.c

#########################################################
#
# potfit Configuration rules
#
#########################################################

HEADERS := ${POTFITHDR}

# serial or mpi
ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
SOURCES	:= ${POTFITSRC} ${MPISRC}
else
SOURCES	:= ${POTFITSRC}
endif

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
CFLAGS  += -DEAM
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

ifneq (,$(findstring 0,${INTERACTION}))
ERROR += No interaction model specified
endif

# EVO - for differential evolution
ifneq (,$(findstring evo,${MAKETARGET}))
CFLAGS += -DEVO
endif

# APOT - for analytic potentials
ifneq (,$(findstring apot,${MAKETARGET}))
CFLAGS += -DAPOT -DNORESCALE
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
CFLAGS += -DPARABEL
endif

ifneq (,$(findstring wzero,${MAKETARGET}))
CFLAGS += -DWZERO
endif

ifneq (,$(findstring dist,${MAKETARGET}))
ifeq (,$(findstring MPI,${PARALLEL}))
CFLAGS += -DPDIST
else
ERROR += "dist is not mpi parallelized -- "
endif
endif

ifneq (,$(findstring newscale,${MAKETARGET}))
ifeq (,$(findstring MPI,${PARALLEL}))
CFLAGS += -DNEWSCALE
else
ERROR += "newscale is not mpi parallelized -- "
endif
endif

ifneq (,$(findstring fweight,${MAKETARGET}))
CFLAGS += -DFWEIGHT
endif

ifneq (,$(findstring acml,${MAKETARGET}))
CFLAGS += -DACML
endif

ifneq (,$(findstring noresc,${MAKETARGET}))
CFLAGS += -DNORESCALE
endif

# Substitute .o for .c to get the names of the object files
OBJECTS := $(subst .c,.o,${SOURCES})

###########################################################################
#
# 	Check for bzr binary
#
###########################################################################

ifeq (Found,$(shell if `which bzr >& /dev/null`; then echo Found; fi))
	BAZAAR = 1
else
	BAZAAR = 0
endif

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
	${CC} ${CFLAGS} ${CINCLUDE} -c powell_lsq.c

# special rules for function evaluation
functions.o: functions.c
	${CC} ${CFLAGS} ${CINCLUDE} -c functions.c

# generic compilation rule
.c.o:
ifeq (,$@)
ifeq (,${MAKETARGET})
	@echo -e "Usage:"
	@echo -e "  make potfit_[interaction]_[options]\n"
	@echo "For more details on compiling potfit please look at the Makefile"
	@exit
endif
else
	${CC} ${CFLAGS} -c $<
endif

# How to link
${MAKETARGET}: ${OBJECTS}
	${CC} ${LIBS} ${LFLAGS_${PARALLEL}} -o $@ ${OBJECTS}
ifneq (,${STRIP})
ifeq (,$(findstring debug,${MAKETARGET}))
	${STRIP} --strip-unneeded -R .comment $@
endif
endif
	${MV} $@ ${BIN_DIR}; rm -f $@

# First recursion only set the MAKETARGET Variable
.DEFAULT:
ifneq (,${CC})
	${MAKE} MAKETARGET='$@' STAGE2
else
	@echo "There is no compiler defined for this option."
	@echo -e "Please adjust the Makefile.\n"
	@exit
endif

potfit:
	@echo -e "\nError:\tYou cannot compile potfit without any options."
	@echo -e "\tAt least an interaction is required.\n"

# Second recursion sets MAKETARGET variable and compiles
# An empty MAKETARGET variable would create an infinite recursion, so we check
STAGE2:
ifneq (,${ERROR})
	@echo -e "\nError: ${ERROR}\n"
else
ifneq (,${MAKETARGET})
	@echo "${WARNING}"
ifeq (0,${RELEASE})
ifeq (1,${BAZAAR})
	@echo -e "Writing bazaar data to version.h\n"
	@rm -f version.h
	@bzr version-info --custom \
	--template="#define VERSION_INFO \"potfit-{branch_nick} (r{revno})\"\n" > version.h
	@bzr version-info --custom \
	--template="#define VERSION_DATE \"{build_date}\"\n" >> version.h
else
	@echo -e "Writing fake bazaar data to version.h\n"
	@rm -f version.h
	@echo -e "#define VERSION_INFO \"potfit-`basename ${PWD}` (r ???)\"" > version.h
	@echo -e "#define VERSION_DATE \"`date +%Y-%m-%d\ %H:%M:%S\ %z`\"" >> version.h
endif
else
	@rm -f version.h
	@echo -e "#define VERSION_INFO \"potfit-setdatehere\"" > version.h
	@echo -e "#define VERSION_DATE \"`date +%Y-%m-%d\ %H:%M:%S\ %z`\"" >> version.h
endif
	${MAKE} MAKETARGET='${MAKETARGET}' ${MAKETARGET}
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
	rm -f *.o *.u *~ \#* *.V *.T *.O *.il

help:
	@echo "Usage: make potfit[_<parallel>][_<option>[_<option>...]]"

