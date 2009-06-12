#!/bin/sh
####################################################################
#
#   Copyright 2009 Daniel Schopf
#             Institute for Theoretical and Applied Physics
#             University of Stuttgart, D-70550 Stuttgart, Germany
#             http://www.itap.physik.uni-stuttgart.de/
#
####################################################################
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
#   along with potfit; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin St, Fifth Floor,
#   Boston, MA  02110-1301  USA
#
#/****************************************************************
#* $Revision: 1.2 $
#* $Date: 2009/06/12 09:10:40 $
#*****************************************************************/

if [ $# -eq 0 ]; then
	echo "Usage: get_config FILENAME"
	exit
elif [ $# -ge 2 ]; then
	echo "More then 1 file given on command line. Only $1 will be interpreted"
	echo
fi
n_conf=`grep -e \#E $1 | wc -l`
n_names=`grep -e generated $1 | wc -l`
if [ $n_conf -gt $n_names ]; then
grep -e generated -e \#E $1 | awk -v file=$1 'BEGIN{i=0}
{
	printf "%d ",i;
	if ($1=="#E") {
		i++;
		printf "no information found in %s\n",file;
	}
	else {
		printf "%s\n",$7;
		getline;
		i++;
	}
}'
else
grep -e generated $1 | awk 'BEGIN{i=0} {print i" "$7;i++; }'
fi
