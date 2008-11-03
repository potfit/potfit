############################################################################
#
# potfit -- The ITAP Force Matching Program
#
# Copyright 2002-2008 Institute for Theoretical and Applied Physics,
# University of Stuttgart, D-70550 Stuttgart
#
# $Revision: 1.41 $
# $Date: 2008/11/03 10:28:32 $
#
############################################################################
#
#     This file is part of potfit.
#
#     potfit is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
#
#     potfit is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with potfit; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin St, Fifth Floor,
#     Boston, MA 02110-1301, USA
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
# The parallelization method <parallel> can be one of:
#
#    mpi   compile for parallel execution, using MPI
#    omp   compile for parallel execution, using OpenMP
#    ompi  compile for parallel execution, using OpenMP and MPI
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
# The compilation process requires the environment variable IMDSYS to
# be set to a recognized value. It specifies what system you have, and
# what compiler you are using. The flags for the compiler and the linker
# are then selected as a function of this variable. It is also possible
# to pass the value of IMDSYS on the command line, e.g.:
#
#   make IMDSYS=P4-icc potfit_mpi_eam
#
# Another important ingredient is the parallelization method, which is
# determined from the make target. The parallelization method is stored
# in the variable PARALLEL, which takes as value one of SERIAL, MPI,
# OMP, OMPI, or PACX.
#
# Depending on the value of ${IMDSYS}, a number of variables must be
# set, from which everything else is constructed.
#
# CC_${PARALLEL} defines the compiler to be used for parallelization
# method ${PARALLEL}. If not defined, the parallelization method
# ${PARALLEL} is not available.
#
# BIN_DIR defines the directory where the potfit binary is put. Note that
# this directory must exist.
#
# MV defines the program used to move the potfit binary to ${BIN_DIR}.
# The default is mv, which is usually ok.
#
# The compilation options are stored in the variable CFLAGS.
# The initial value of CFLAGS is set to the variable FLAGS,
# which can be given on the command line as explained above for
# IMDSYS, although this is usually not necessary.
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
# You may have to change the setting for an existing value of IMDSYS.
# or you have to add support for a new value of IMDSYS. The latter is
# best done by using the folloing template for IMDSYS=sys-cc:
#
# ifeq (sys-cc,${IMDSYS})
#   CC_SERIAL		= serial-compiler
#   CC_OMP		= OpenMP-compiler
#   CC_MPI		= MPI-compiler
#   CC_OMPI		= OpenMP+MPI-compiler
#   BIN_DIR		= ${HOME}/bin
#   OPT_FLAGS		+= generic flags for optimization
#   OPT_MPI_FLAGS	+= MPI-specific flags for optimization
#                          similar variables for other parallelizations
#   MPI_FLAGS		+= MPI-specific flags
#                          similar variables for other parallelizations
#   DEBUG_FLAGS		+= generic flags for debugging
#   DEBUG_MPI_FLAGS	+= MPI-specific flags for debugging
#                          similar variables for other parallelizations
#   PROF_FLAGS		+= flags for profiling
#   LIBS		+= generically needed libraries
#   MPI_LIBS		+= MPI-specific libraries
#                          similar variables for other parallelizations
#   PROF_LIBS		+= libraries for profiling
# endif
#
# Variables remaining empty need not be mentioned.


###########################################################################
#
#  Defaults for some variables
#
###########################################################################

MV		= mv      # program to move imd binary to ${BIN_DIR}
LIBS		+= -lm
MPI_FLAGS	+= -DMPI
OMP_FLAGS	+= -DOMP
OMPI_FLAGS	+= -DMPI -DOMP
DEBUG_FLAGS	+= -DDEBUG
MKLDIR          =  /common/linux/paket/intel/mkl91/
MKLPATH         =  ${MKLDIR}/lib/
CINCLUDE        =  -I${MKLDIR}/include/

###########################################################################
#
#  flags for 64bit-Linux
#
###########################################################################

# AMD Opteron, gcc3
ifeq (x86_64-gcc3,${IMDSYS})
  CC_SERIAL     = gcc
  CC_MPI        = mpicc
  MPICH_CC      = gcc
#  MPICH_CLINKER = gcc
  BIN_DIR       = ${HOME}/bin/${HOSTTYPE}
#  FFTW_DIR     = /common/linux/paket/fftw-3.0.1
  OPT_FLAGS     += -O -march=opteron -Wno-unused
  DEBUG_FLAGS   += -g -O -Wall
  PROF_FLAGS    += -g3 -pg
  LFLAGS        +=  -static
#  ACMLPATH      = /common/linux/paket/acml3.5.0/gnu64
  MKLPATH       = ${MKLDIR}/lib/em64t/
#  CINCLUDE     += -I$(ACMLPATH)/include
#  CINCLUDE      = -I${MKLDIR}/include
#  LD_LIBRARY_PATH +=':$(ACMLPATH)/lib:'
  export        MPICH_CC # MPICH_CLINKER
  export        LD_LIBRARY_PATH
# acml
#  LIBS		:= $(ACMLPATH)/lib/libacml.a \
#		   -L${ACMLPATH}/lib -lpthread -lacml -lg2c
# intel mkl
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_em64t.a \
		   -L${MKLPATH} -lguide -lpthread
endif

# Intel EM64T "AMD inside", icc
ifeq (x86_64-icc,${IMDSYS})
  CC_SERIAL     = icc
  CC_MPI        = mpicc
  MPICH_CC      = icc
  MPICH_CLINKER = icc
  BIN_DIR       = ${HOME}/bin/${HOSTTYPE}
  OPT_FLAGS     += -O3 -ip # -fno-builtin # -axP # remove -axP for Opteron
  MPI_FLAGS     +=
  OMP_FLAGS     += -openmp
  OMPI_FLAGS    += -openmp
  DEBUG_FLAGS   += -g  -Wall
  PROF_FLAGS    += -prof_gen
  RCD_FLAGS     += # -DRCD -rcd
  MPI_LIBS      +=
  LFLAGS        += -i-static -openmp
# acml
#   ACMLPATH      = /common/linux/paket/acml3.5.0/gnu64
#   CINCLUDE     += -I$(ACMLPATH)/include
#   LD_LIBRARY_PATH +=':$(ACMLPATH)/lib:'
#   LIBS		:= $(ACMLPATH)/lib/libacml.a \
# 		   -L${ACMLPATH}/lib -lpthread -lacml -lg2c
# intel mkl
  MKLPATH       = ${MKLDIR}/lib/em64t/
  LIBS		+= ${MKLPATH}/libmkl_lapack.a  ${MKLPATH}/libmkl_em64t.a \
		   -L${MKLPATH} -lguide -lpthread
  export        MPICH_CC MPICH_CLINKER
endif

# Athlon MP or XP, icc
ifeq (AthlonMP-icc,${IMDSYS})
  CC_SERIAL	= icc
  CC_OMP	= icc
  CC_MPI	= mpicc
  CC_OMPI	= mpicc
  MPICH_CC      = icc
  MPICH_CLINKER = icc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  FFTW_DIR 	= /common/linux/paket/fftw-3.0.1
  OPT_FLAGS	+= -O -ip -tpp6 -axK #-static
  OMP_FLAGS	+= -openmp
  OMPI_FLAGS	+= -openmp
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -prof_gen
  RCD_FLAGS	+= -DRCD -rcd
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread

  export MPICH_CC MPICH_CLINKER
endif

# Athlon MP or XP, gcc3
ifeq (AthlonMP-gcc3,${IMDSYS})
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  MPICH_CC      = gcc
  MPICH_CLINKER = gcc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  FFTW_DIR 	= /common/linux/paket/fftw-3.0.1
  OPT_FLAGS	+= -O -march=athlon-mp # -static
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
  export        MPICH_CC #MPICH_CLINKER
endif

# Pentium 4 or Xeon, icc
ifeq (P4-icc,${IMDSYS})
  CC_SERIAL	= icc
  CC_OMP	= icc
  CC_MPI	= mpicc
  CC_OMPI	= mpicc
  MPICH_CC      = icc
  MPICH_CLINKER = icc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  FFTW_DIR 	= /common/linux/paket/fftw-3.0.1
  OPT_FLAGS	+= -O -ip -tpp7 # -static
  OMP_FLAGS	+= -openmp
  OMPI_FLAGS	+= -openmp
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -prof_gen
  RCD_FLAGS	+= -DRCD -rcd
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
  export        MPICH_CC MPICH_CLINKER
endif

# Pentium 4 or Xeon, gcc3
ifeq (P4-gcc3,${IMDSYS})
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  MPICH_CC      = gcc
  MPICH_CLINKER = gcc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  FFTW_DIR 	= /common/linux/paket/fftw-3.0.1
  OPT_FLAGS	+= -O -march=pentium4 # -static
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
  export        MPICH_CC # MPICH_CLINKER
endif

# Pentium III, icc
ifeq (P3-icc,${IMDSYS})
  CC_SERIAL	= icc
  CC_OMP	= icc
  CC_MPI	= mpicc
  CC_OMPI	= mpicc
  MPICH_CC      = icc
  MPICH_CLINKER = icc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  FFTW_DIR 	= /common/linux/paket/fftw-3.0.1
  OPT_FLAGS	+= -O -ip -tpp6 -axK # -static
  OMP_FLAGS	+= -openmp
  OMPI_FLAGS	+= -openmp
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -prof_gen
  RCD_FLAGS	+= -DRCD -rcd
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
  export        MPICH_CC MPICH_CLINKER
endif

# Pentium III, gcc3
ifeq (P3-gcc3,${IMDSYS})
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  MPICH_CC      = gcc
  MPICH_CLINKER = gcc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  FFTW_DIR 	= /common/linux/paket/fftw-3.0.1
  OPT_FLAGS	+= -O -march=pentium3 # -static
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
  export        MPICH_CC # MPICH_CLINKER
endif

# generic ia32 CPU, gcc 2.95 - slow!
ifeq (ia32-gcc2,${IMDSYS})
  CC_SERIAL	= gcc
  CC_MPI	= mpicc
  MPICH_CC      = gcc
  MPICH_CLINKER = gcc
  BIN_DIR	= ${HOME}/bin/${HOSTTYPE}
  FFTW_DIR 	= /common/linux/paket/fftw-3.0.1
  OPT_FLAGS	+= -O3
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -g3 -pg
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
  export        MPICH_CC MPICH_CLINKER
endif

###########################################################################
#
#  flags for IA64-Linux
#
###########################################################################

# Itanium 2, ecc
ifeq (ia64-ecc,${IMDSYS})
  CC_SERIAL	= ecc
  CC_OMP	= ecc
  CC_MPI	= mpicc
  CC_OMPI	= mpicc
  BIN_DIR	= ${HOME}/bin/ia64
  OPT_FLAGS	+= -O -ipo #-static
  OMP_FLAGS	+= -openmp
  OMPI_FLAGS	+= -openmp
  PROF_FLAGS	+= -prof_gen
  DEBUG_FLAGS	+= -g
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
endif

###########################################################################
#
#  flags for alpha
#
###########################################################################

# alpha ev56 or higher, cc
ifeq (alpha-cc,${IMDSYS})
  CC_SERIAL	= cc
  CC_OMP	= cc
  BIN_DIR	= ${HOME}/bin/alpha
  OPT_FLAGS	+= -DALPHA -O3 -float -fp_reorder -arch ev56 -tune host
  OMP_FLAGS	+= -mp
  OMPI_FLAGS	+= -mp
  PROF_FLAGS	+= -prof_gen
  DEBUG_FLAGS	+= -g3 -pg
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
endif

# alpha, gcc2 - slow!
ifeq (alpha-gcc2,${IMDSYS})
  CC_SERIAL	= gcc
  BIN_DIR	= ${HOME}/bin/$alpha
  OPT_FLAGS	+= -DALPHA -O3
  PROF_FLAGS	+= -g3 -pg
  DEBUG_FLAGS	+= -g
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
endif

###########################################################################
#
#  flags for IRIX
#
###########################################################################

# irix-mips3, cc
ifeq (irix-cc,${IMDSYS})
  CC_SERIAL	= cc
  CC_OMP	= cc
  BIN_DIR	= ${HOME}/bin/iris4d
  OPT_FLAGS	+= -Dsgi -O3 -n32 -mips3 -xansi -woff 1174
  OMP_FLAGS	+= -mp
  OMPI_FLAGS	+= -mp
  PROF_FLAGS	+= -g3
  DEBUG_FLAGS	+= -Dsgi -g  -n32 -mips3 -xansi -woff 1174
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
endif

# irix, gcc2 - slow!
ifeq (irix-gcc2,${IMDSYS})
  CC_SERIAL	= gcc
  BIN_DIR	= ${HOME}/bin/iris4d
  OPT_FLAGS	+= -Dsgi -O3
  PROF_FLAGS	+= -g3 -pg
  DEBUG_FLAGS	+= -Dsgi -g
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
endif

###########################################################################
#
#  flags for sparc
#
###########################################################################

# UltraSparc III, cc
ifeq (USparc3-cc,${IMDSYS})
  CC_SERIAL	= cc
  CC_OMP	= cc
  CC_MPI	= mpcc
  BIN_DIR	= ${HOME}/bin/sparc
  MPI_LIBS	+= -lmpi
  OPT_FLAGS	+= -fast -xtarget=ultra3 -xarch=v9b
  OMP_FLAGS	+= -xopenmp
  OMPI_FLAGS	+= -xopenmp
  PROF_FLAGS	+= -p
  DEBUG_FLAGS	+= -g -xO3
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
endif

# sparc, gcc2 - slow!
ifeq (sparc-gcc2,${IMDSYS})
  CC_SERIAL	= gcc
  BIN_DIR	= ${HOME}/bin/sparc
  OPT_FLAGS	+= -O3
  PROF_FLAGS	+= -g3 -pg
  DEBUG_FLAGS	+= -g
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
endif

###########################################################################
#
#  flags for T3E
#
###########################################################################

# Cray T3E, cc
ifeq (T3E-cc,${IMDSYS})
  CC_MPI	= cc
  CC_PACX	= cc
  BIN_DIR	= ${HOME}/bin/t3e
  MPI_LIBS	+= -lmpi
  PROF_LIBS	+= -lapp
  OPT_FLAGS	+= -Dt3e -O3 -htolerant,aggress,report=isf
  PROF_FLAGS	+= -Gf -happrentice
  DEBUG_FLAGS	+= -Dt3e -g
  PACX_DIR	= ${HOME}/WORK/PACX
  PACX_LIBS	+= -L ${PACX_DIR}/lib -lpacx -llzo -lmpi
  PACX_FLAGS	+= -I${PACX_DIR}/include
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
endif

###########################################################################
#
#  flags for Hitachi SR8000
#
###########################################################################

# Hitachi sr8k, cc
ifeq (sr8k-cc,${IMDSYS})
  CC_MPI		= mpicc
  CC_OMP		= cc
  CC_OMPI		= mpicc
  BIN_DIR		= ${HOME}/bin/SR8000
  OPT_FLAGS		+= -O4 +Op -msg=e
  OPT_MPI_FLAGS		+= -noparallel
  OPT_OMP_FLAGS		+= -omp -par
  OPT_OMPI_FLAGS	+= -omp -par
  DEBUG_FLAGS		+= -g
  DEBUG_OMP_FLAGS	+= -omp -par=1 -O2
  DEBUG_OMPI_FLAGS	+= -omp -par   -O2
  PROF_FLAGS		+= -Xfuncmonitor
  PROF_LIBS		+= -lpl
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
endif

# Hitachi sr8k, xcc cross compiler
ifeq (sr8k-xcc,${IMDSYS})
  CC_MPI		= xmpicc
  CC_OMP		= xcc
  CC_OMPI		= xmpicc
  BIN_DIR		= hwwfs1:sr8k/bin/SR8000
  MV			= scp   # we move binary to different machine
  OPT_FLAGS		+= -O4 +Op -msg=e
  OPT_MPI_FLAGS		+= -noparallel
  OPT_OMP_FLAGS		+= -omp -par
  OPT_OMPI_FLAGS	+= -omp -par
  DEBUG_FLAGS		+= -g
  DEBUG_OMP_FLAGS	+= -omp -par=1 -O2
  DEBUG_OMPI_FLAGS	+= -omp -par   -O2
  PROF_FLAGS		+= -Xfuncmonitor
  PROF_LIBS		+= -lpl
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
endif

###########################################################################
#
#  flags for IBM SP
#
###########################################################################

# Power4 regatta
ifeq (Power4-cc,${IMDSYS})
  CC_SERIAL	= xlc
  CC_MPI	= mpcc
  CC_OMP	= xlc_r
  CC_OMPI	= mpcc_r
  BIN_DIR	= ${HOME}/bin/powerpc
  OPT_FLAGS	+= -O4 -w
  OMP_FLAGS	+= -DUSE_WALLTIME -qsmp=omp
  OMPI_FLAGS	+= -qsmp=omp
  DEBUG_FLAGS	+= -g
  PROF_FLAGS	+= -p
  PROF_LIBS	+= -lpl
  LIBS		+= ${MKLPATH}/libmkl_lapack.a ${MKLPATH}/libmkl_ia32.a \
		   -L${MKLPATH} -lguide -lpthread
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
# OpenMP
ifneq (,$(strip $(findstring omp,${MAKETARGET})))
PARALLEL = OMP
endif
# MPI + OpenMP
ifneq (,$(strip $(findstring ompi,${MAKETARGET})))
PARALLEL = OMPI
endif
# PACX
ifneq (,$(strip $(findstring pacx,${MAKETARGET})))
PARALLEL = PACX
endif


###########################################################################
#
#  Compiler, flags, libraries
#
###########################################################################

# compiler; if empty, we issue an error later
CC = ${CC_${PARALLEL}}

# optimization flags
OPT_FLAGS   += ${${PARALLEL}_FLAGS} ${OPT_${PARALLEL}_FLAGS}
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

POTFITHDR   	= potfit.h powell_lsq.h utils.h
POTFITSRC 	= utils.c bracket_r.c powell_lsq.c brent_r.c \
		  linmin_r.c force.c \
		  config.c param.c potential.c potfit.c \
		  splines.c simann.c rescale.c functions.c smooth.c
MPISRC          = mpi_utils.c

#########################################################
#
# potfit Configuration rules
#
#########################################################

HEADERS := ${POTFITHDR}

# 3d, serial or mpi
ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
SOURCES	:= ${POTFITSRC} ${MPISRC}
else
SOURCES	:= ${POTFITSRC}
endif

###  INTERACTION  #######################################



########### HERE COMES POTFIT
# EAM2 or EAM  -  this is now the same
ifneq (,$(strip $(findstring eam,${MAKETARGET})))
CFLAGS  += -DEAM
endif

# APOT - for analytic potentials
ifneq (,$(findstring apot,${MAKETARGET}))
ifeq (,$(findstring eam,${MAKETARGET}))
CFLAGS += -DAPOT
else
ERROR += "apot and eam are not compatible yet -- "
endif
ifneq (,$(findstring mpi,${MAKETARGET}))
ERROR += "apot and mpi are not compatible yet -- "
endif
endif

# Stress
ifneq (,$(findstring stress,${MAKETARGET}))
CFLAGS += -DSTRESS
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
	${CC} ${CFLAGS} -c $<

# How to link
${MAKETARGET}: ${OBJECTS}
	${CC} ${LFLAGS} -o $@ ${OBJECTS} ${LIBS}
	${MV} $@ ${BIN_DIR}; rm -f $@

# First recursion only set the MAKETARGET Variable
.DEFAULT:
ifneq (,${CC})
	${MAKE} MAKETARGET='$@' STAGE2
else
ifneq (,${IMDSYS})
	@echo "IMDSYS variable ${IMDSYS} is not recognized"
else
	@echo "IMDSYS variable is not set"
endif
endif

# Second recursion sets MAKETARGET variable and compiles
# An empty MAKETARGET variable would create an infinite recursion, so we check
STAGE2:
ifneq (,${ERROR})
	@echo "${ERROR}"
else
ifneq (,${MAKETARGET})
	@echo "${WARNING}"
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
	@echo "Usage: gmake potfit[_<parallel>][_<option>[_<option>...]]"

