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

include Makefile.inc

###########################################################################
#
# catch default targets
#
###########################################################################

.DEFAULT:
ifneq (,${CC})
	@${MAKE} -C src --no-print-directory MAKETARGET=${@} $@
else
	@echo -e "\nError: There is no compiler defined for your SYSTEM variable (${SYSTEM})."
	@echo -e "Please adjust the SYSTEM variable in the Makefile.inc file.\n"
endif

potfit:
	@echo -e "\nError:\tYou cannot compile potfit without any options."
	@echo -e "\tAt least an interaction is required.\n"
	@${MAKE} --no-print-directory help

###########################################################################
#
# basic targets for help and cleanup
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
