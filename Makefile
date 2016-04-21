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

include Makefile.inc

###########################################################################
#
#  Catch default targets
#
###########################################################################

potfit:
	@echo -e "\nError:\tYou cannot compile potfit without any options."
	@echo -e "\tAt least an interaction is required.\n"

STAGE2:
	@echo -e ""
	@echo -e ""

.DEFAULT:
ifneq (,${CC})
	@${MAKE} -C src --no-print-directory MAKETARGET='$@' STAGE2
else
	@echo "There is no compiler defined for this option."
	@echo -e "Please adjust the Makefile.\n"
	@exit
endif


###########################################################################
#
# Basic targets for help and cleanup
#
###########################################################################

clean:
	@make -C src --no-print-directory $@
	@rm -f *.o *.u *~ \#* *.V *.T *.O *.il

distclean:
	@make -C src --no-print-directory $@
	@rm -f *.o *.u *~ \#* *.V *.T *.O *.il
	@rm -f bin/*

help:
	@echo "Usage: make potfit[_<parallel>][_<option>[_<option>...]]"

