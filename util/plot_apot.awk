#!/usr/bin/awk -f
#####################################################################
#
# plot_pot.awk: create gnuplot readable potential from
#		from analytic potential file format
#
####################################################################
#
#   Copyright 2008 Daniel Schopf
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
####################################################################
# $Revision: 1.3 $
# $Date: 2008/11/26 14:18:56 $
####################################################################
#
# Usage: plot_pot.awk potfit_apot.pot
#
# The resulting potential is written to standard output.
#
# ATTENTION plot_pot.awk is PRE-ALPHA. No valdiation whatsoever!!!
#
####################################################################

BEGIN {
	count=0;
	maxdist=0;
	ORS="";
}

{
	for (a=0;a<(ARGC-1);a++) {
		while (substr($0,2,1)!="F") getline;
		if ($2 != 0) {
			print "Error - wrong potential format of " ARGV[ARGIND]  ;
			exit 2;
		}
		total_pots=$3;
		if (int(total_pots)!=total_pots) {
			print "ERROR - incorrect parameter file " ARGV[ARGIND]  ;
			exit 2;
		}
		for (i=count;i<(count+total_pots);i++){
			while (substr($0,1,4)!="type") getline;
			pot_name[i] = $2;
			if ((x=index(pot_name[i],"_"))>0)
				pot_name[i] = substr(pot_name[i],1,x-1);
			if (pot_name[i]=="eopp")
				n_param[i]=6;
			getline;
			if ($2>maxdist)
				maxdist=$2;
			getline;
			if (substr($0,1,1)=="#") getline;
			for (l=1;l<=n_param[i];l++) {
				params[i "," l ] = $2;
				getline
			}
		}
		count = count + total_pots;
	}
	nextfile;
}

END {
	print "reset;" > "plot";
	print "set grid;" > "plot";
	print "set arrow 1 from " maxdist ",.2 to " maxdist ",-.2 nohead size 2,15,10 lw 2;" > "plot";
	print "set label \"cutoff\" at " maxdist*0.95 ",0.23;" > "plot";
	print "pl [0.1:" maxdist*1.1 "][-0.5:1.5]" > "plot";
	for (i=0;i<count;i++) {
		if (pot_name[i] == "eopp") {
			printf "%f/x**%f+%f/x**%f*cos(%f*x+%f) w l",params[i","1],params[i","2],params[i","3],params[i","4],params[i","5],params[i","6] > "plot";
		}
		if (i!=(count-1))
			print "," > "plot";
	}
	print ";" > "plot";
	system("gnuplot -persist plot");
}
