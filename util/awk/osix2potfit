#!/usr/bin/awk -f
#################################################################
#
# osix2potfit: convert osix potential format to potfit format 0
#
#################################################################
#
#   Copyright 2009 Daniel Schopf
#             Institute for Theoretical and Applied Physics
#             University of Stuttgart, D-70550 Stuttgart, Germany
#             http://www.itap.physik.uni-stuttgart.de/
#
#################################################################
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
####################################################################

{
	species = $1
	pairs = $2
	if ( (species*(species+1)/2) != pairs ) {
		print "Number of pairs is wrong, aborting";
		exit 2;
	}
	getline;
	for (i=0;i<species;i++) {
		element[i]=$4;
		getline;
	}
	# get comment line
	name[0]=$4;
	name[1]=$5;
	name[2]=$6;
	name[3]=$7;
	name[4]=$8;
	name[5]=$9;
	getline;
	for (i=0;i<pairs;i++) {
		pot_i[i]=$1;
		pot_j[i]=$2;
		par_0[i]=$3;
		par_1[i]=$4;
		par_2[i]=$5;
		par_3[i]=$6;
		par_4[i]=$7;
		par_5[i]=$8;
		getline;
	}
}

END {

	printf "#F 0 %d\n",pairs;
	printf "#C";
	for (i=0;i<species;i++)
		printf " %s",element[i];
	printf "\n##";
	for (i=0;i<pairs;i++)
		printf " %s-%s",element[pot_i[i]-1],element[pot_j[i]-1];
	printf "\n#I";
	for (i=0;i<pairs;i++)
		printf " 0"
	printf "\n#E\n\n";

	for (i=0;i<pairs;i++) {
		printf "type eopp\n";
		printf "cutoff 7\n";
		printf "%s %f %d %d\n",name[0],par_0[i],1,10000;
		printf "%s %f %d %d\n",name[1],par_1[i],1,20;
		printf "%s %f %d %d\n",name[2],par_2[i],-100,100;
		printf "%s %f %d %d\n",name[3],par_3[i],2,20;
		printf "%s %f %d %d\n",name[4],par_4[i],0,6;
		printf "%s %f %d %f\n",name[5],par_5[i],0,6.3;
		if (i!=(pairs-1)) 
			printf "\n";
	}
}
