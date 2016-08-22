###########################################################################
#
# potfit - open-source force-matching
#
# Copyright 2002-2016 - the potfit development team
#
# http://potfit.sourceforge.net/
#
###########################################################################
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
#
############################################################################

SHELL = /bin/bash

###########################################################################
#
#  Adjust these variables to your system
#
###########################################################################

# Currently the following systems are available:
# x86_64-icc    64bit Intel Compiler
# x86_64-gcc    64bit GNU Compiler
# x86_64-clang  64bit LLVM Compiler
# i686-icc      32bit Intel Compiler
# i686-gcc      32bit GNU Compiler
# i686-kim  	  32bit GNU Compiler

#SYSTEM = i686-gcc
#SYSTEM = x86_64-gcc
SYSTEM = i686-kim

# Base directory of your installation of the MKL or ACML

MKLDIR      = /opt/intel/composer_xe_2015.1.133/mkl
ACML4DIR    = /opt/acml4.4.0/gfortran64
ACML5DIR    = /opt/acml/gfortran64
LIBMDIR     = /opt/acml/libm

###########################################################################
#
#  Do not change anything below unless you know what you are doing
#
###########################################################################

############################################################################
#
#  flags for 64bit
#
############################################################################

ifeq (x86_64-icc,${SYSTEM})
# compiler
  CC_SERIAL     = icc
  CC_MPI        = mpicc
  OMPI_CC       = icc
  OMPI_CLINKER  = icc

# general optimization flags
  OPT_FLAGS     += -fast -xHost

# profiling and debug flags
  PROF_FLAGS    += --profile-functions
  PROF_LIBS     += --profile-functions
  DEBUG_FLAGS   += -g -Wall

# Intel Math Kernel Library
ifeq (,$(strip $(findstring acml,${MAKETARGET})))
  CINCLUDE 	+= -I${MKLDIR}/include
  LIBS 		+= -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential \
		   -lmkl_core -Wl,--end-group -lpthread
endif

# AMD Core Math Library
ifneq (,$(strip $(findstring acml4,${MAKETARGET})))
  CINCLUDE 	+= -I${ACML4DIR}/include
  LIBS		= -L${ACML4PATH} -lpthread -lacml -lacml_mv
endif
ifneq (,$(strip $(findstring acml5,${MAKETARGET})))
   LIBMPATH 	= ${LIBMDIR}/lib/dynamic
   CINCLUDE     += -I${ACML5DIR}/include -I${LIBMDIR}/include
   LIBS		+= -L${ACML5PATH} -L${LIBMPATH} -lpthread -lacml -lamdlibm
endif

 export        OMPI_CC OMPI_CLINKER
endif

ifeq (x86_64-gcc,${SYSTEM})
# compiler
  CC_SERIAL     = gcc
  CC_MPI        = mpicc
  OMPI_CC       = gcc
  OMPI_CLINKER  = gcc

# general optimization flags
  OPT_FLAGS     += -m64 -O3 -march=native -std=c99

# profiling and debug flags
  PROF_FLAGS    += -m64 -g3 -pg
  PROF_LIBS     += -m64 -g3 -pg
  DEBUG_FLAGS   += -m64 -g3 -Wall -Werror -pedantic -std=c99

# Intel Math Kernel Library
ifeq (,$(strip $(findstring acml,${MAKETARGET})))
  CINCLUDE      += -I${MKLDIR}/include
  LIBS 		+= -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
		   -Wl,--end-group -lpthread -Wl,--as-needed
endif

# AMD Core Math Library
ifneq (,$(strip $(findstring acml4,${MAKETARGET})))
  CINCLUDE     	+= -I${ACML4DIR}/include
  LIBS		+= -L${ACML4PATH} -lpthread -lacml -lacml_mv -Wl,--as-needed
endif
ifneq (,$(strip $(findstring acml5,${MAKETARGET})))
  LIBMPATH 	= ${LIBMDIR}/lib/dynamic
  CINCLUDE     	+= -I${ACML5DIR}/include -I${LIBMDIR}/include
  LIBS		+= -L${ACML5PATH} -L${LIBMPATH} -lpthread -lacml -lamdlibm -Wl,--as-needed
endif

 export        OMPI_CC OMPI_CLINKER
endif

ifeq (x86_64-clang,${SYSTEM})
# compiler
  CC_SERIAL     = clang
  CC_MPI        = mpicc
  OMPI_CC       = clang
  OMPI_CLINKER  = clang

# general optimization flags
  OPT_FLAGS     += -O3 -march=native -std=c99

# profiling and debug flags
  PROF_FLAGS    += -g
  PROF_LIBS     += -g
  DEBUG_FLAGS   += -g -Wall -Werror -pedantic -std=c99
  ASAN_FLAGS 	+= -g -fsanitize=address -fno-omit-frame-pointer
  ASAN_LFLAGS 	= -g -fsanitize=address

# Intel Math Kernel Library
ifeq (,$(strip $(findstring acml,${MAKETARGET})))
  CINCLUDE      += -I${MKLDIR}/include
  LIBS 		+= -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
		   -Wl,--end-group -lpthread -Wl,--as-needed
endif

# AMD Core Math Library
ifneq (,$(strip $(findstring acml4,${MAKETARGET})))
  CINCLUDE     	+= -I${ACML4DIR}/include
  LIBS		+= -L${ACML4PATH} -lpthread -lacml -lacml_mv -Wl,--as-needed
endif
ifneq (,$(strip $(findstring acml5,${MAKETARGET})))
  LIBMPATH 	= ${LIBMDIR}/lib/dynamic
  CINCLUDE     	+= -I${ACML5DIR}/include -I${LIBMDIR}/include
  LIBS		+= -L${ACML5PATH} -L${LIBMPATH} -lpthread -lacml -lamdlibm -Wl,--as-needed
endif

 export        OMPI_CC OMPI_CLINKER
endif

###########################################################################
#
#  flags for 32bit
#
###########################################################################

ifeq (i686-icc,${SYSTEM})
# compiler
  CC_SERIAL	= icc
  CC_MPI	= mpicc
  OMPI_CC       = icc
  OMPI_CLINKER  = icc

# general optimization flags
  OPT_FLAGS	+= -fast -xHost

# profiling and debug flags
  PROF_FLAGS	+= -prof-gen
  PROF_LIBS 	+= -prof-gen
  DEBUG_FLAGS	+= -g -Wall -wd981 -wd1572

# Intel Math Kernel Library
ifeq (,$(strip $(findstring acml,${MAKETARGET})))
  CINCLUDE      += -I${MKLDIR}/include
  LIBS 		+= -Wl,--start-group -lmkl_intel -lmkl_sequential -lmkl_core \
		   -Wl,--end-group -lpthread
endif

# AMD Core Math Library
ifneq (,$(strip $(findstring acml4,${MAKETARGET})))
  CINCLUDE     	+= -I$(ACML4DIR)/include
  LIBS		+= -L${ACML4PATH} -lpthread -lacml
endif
ifneq (,$(strip $(findstring acml5,${MAKETARGET})))
  LIBMPATH 	= ${LIBMDIR}/lib/dynamic
  CINCLUDE     	+= -I$(ACML5DIR)/include -I${LIBMDIR}/include
  LIBS		+= -L${ACML5PATH} -L${LIBMPATH} -lpthread -lacml
endif

  export        OMPI_CC OMPI_CLINKER
endif

ifeq (i686-gcc,${SYSTEM})
# compiler
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  OMPI_CC     	= gcc
  OMPI_CLINKER 	= gcc

# general optimization flags
  OPT_FLAGS	+= -m32 -O3 -march=native -Wno-unused
  LIBS 		+= -m32

# profiling and debug flags
  PROF_FLAGS	+= -m32 -g3 -pg
  PROF_LIBS	+= -m32 -g3 -pg
  DEBUG_FLAGS	+= -m32 -g3 -Wall

# Intel Math Kernel Library
ifeq (,$(strip $(findstring acml,${MAKETARGET})))
  CINCLUDE      += -I${MKLDIR}/include
  LIBS		+= -Wl,--start-group -lmkl_intel -lmkl_sequential -lmkl_core \
		   -Wl,--end-group -lpthread -Wl,--as-needed
endif

# AMD Core Math Library
ifneq (,$(strip $(findstring acml4,${MAKETARGET})))
  CINCLUDE     	+= -I$(ACML4DIR)/include
  LIBS		+= -L${ACML4PATH} -lpthread -lacml -Wl,--as-needed
endif
ifneq (,$(strip $(findstring acml5,${MAKETARGET})))
  LIBMPATH 	= ${LIBMDIR}/lib/dynamic
  CINCLUDE     	+= -I$(ACML5DIR)/include -I${LIBMDIR}/include
  LIBS		+= -L${ACML5PATH} -L${LIBMPATH} -lpthread -lacml -Wl,--as-needed
endif

  export        OMPI_CC OMPI_CLINKER
endif


################################################################################
# KIM
################################################################################
ifeq (i686-kim,${SYSTEM})
  ifeq ($(shell kim-api-build-config --version),)
	    $(error The utility 'kim-api-build-config' not found on path; Update PATH environment variable appropriately)
  endif

  CC_SERIAL      = $(shell kim-api-build-config --cc)
  CINCLUDE      += $(shell kim-api-build-config --includes)
  FLAGS         += $(shell kim-api-build-config --cflags) -std=c99
  LFLAGS_SERIAL += $(shell kim-api-build-config --ldflags)
  LIBS          += $(shell kim-api-build-config --ldlibs)

# Intel Math Kernel Library
ifeq (,$(strip $(findstring acml,${MAKETARGET})))
  CINCLUDE      += -I${MKLDIR}/include
  LIBS 		+= -L${MKLDIR}/lib/intel64 
  LIBS    += -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
					   -Wl,--end-group -lpthread -Wl,--as-needed
endif

# AMD Core Math Library
ifneq (,$(strip $(findstring acml4,${MAKETARGET})))
  CINCLUDE     	+= -I$(ACML4DIR)/include
  LIBS		+= -L${ACML4PATH} -lpthread -lacml -Wl,--as-needed
endif
ifneq (,$(strip $(findstring acml5,${MAKETARGET})))
  LIBMPATH 	= ${LIBMDIR}/lib/dynamic
  CINCLUDE     	+= -I$(ACML5DIR)/include -I${LIBMDIR}/include
  LIBS		+= -L${ACML5PATH} -L${LIBMPATH} -lpthread -lacml -Wl,--as-needed
endif

  export        OMPI_CC OMPI_CLINKER
endif


################################################################################
# LEVMAR 
################################################################################
ifneq (,$(strip $(findstring levmar,${MAKETARGET})))
  LIBS   += -L../lib/geodesicLMv1.1/geodesicLM -lgeodesiclm -lgfortran 
endif  #levmar

