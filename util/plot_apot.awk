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
# $Revision: 1.1 $
# $Date: 2008/11/20 10:06:51 $
####################################################################
#
# Usage: plot_pot.awk potfit_apot.pot
#
# The resulting potential is written to standard output.
#
# ATTENTION plot_pot.awk is PRE-ALPHA. No valdiation whatsoever!!!
#
####################################################################

{
  while (substr($0,2,1)!="F") getline;
  if ($2 != 0) {
    print "Error - wrong potential format of " ARGV[ARGIND]  ;
    exit 2;
  }
  total_pots=$3;
  if (int(total_pots)!=total_pots) { 
    printf "ERROR - incorrect parameter file " ARGV[ARGIND]  ;
    exit 2;
  }
  while (substr($0,1,1)!="#") getline; 
  maxdist = 0;
  for (i=1;i<=total_pots;i++){
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
	getline;
	 for (j=1;j<=n_param[i];j++) {
		params[i "," j ] = $2;
		getline	
		 }
 }
 nextfile;
}

END {
	print "set grid;";
	print "pl [0:" maxdist*1.1 "][:3]";
	for (i=1;i<=total_pots;i++) {
		if (pot_name[i] == "eopp") {
		print params[i","1] "/x**" params[i","2] "+" params[i","3] "/x**" params[i","4] "*cos(" params[i","5] "*x+" params[i","6] ") w l"
		}
		if (i!=total_pots)
			print ",";
	}
}
  
