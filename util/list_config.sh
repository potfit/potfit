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
#* $Revision: 1.5 $
#* $Date: 2010/01/11 09:04:21 $
#*****************************************************************/

if [ $# -eq 0 ]; then
	echo "Usage: list_config.sh FILENAME"
	exit
elif [ $# -ge 2 ]; then
	echo "More then 1 file given on command line. Only $1 will be interpreted"
	echo
fi
if [ ! -e "$1" ]; then
	echo "File $1 does not exists, aborting"
	echo
	exit 2;
fi
n_conf=`grep -e \#E $1 | wc -l`
n_names=`grep -e generated $1 | wc -l`
n_weights=`grep -e \#W $1 | wc -l`
if [ $n_conf -gt $n_names ]; then
grep -n -e \#N -e generated -e \#E $1 | awk -v file=$1 'BEGIN{i=0}
{
	gsub(":"," ");
	line=$1;
	n_atoms=$3;
	getline;
	gsub(":"," ");
	printf "%d ",i;
	if ($2=="#E") {
		i++;
		printf "no information found in %s\n",file;
		printf "\t with %d atoms, starting at line %d\n",n_atoms,line;
	}
	else {
		printf "generated from %s\n",$8;
		printf "\t with %d atoms, starting at line %d\n",n_atoms,line;
		getline;
		i++;
	}
}'
else
grep -n -e \#N -e generated $1 | awk 'BEGIN{i=0} 
{
	gsub(":"," ");
	line=$1;
	n_atoms=$3;
	getline;
	printf i" generated from "$7"\n\t with "n_atoms" atoms, starting at line "line"\n";i++;
}'
fi
