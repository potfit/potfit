#!/usr/bin/awk -f
#####################################################################
#
# combine_eam.awk: combine pair, transfer, embedding to EAM potential.
#
####################################################################
# 
#   Copyright 2008 Peter Brommer
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
# $Date: 2008/10/09 18:03:01 $
####################################################################
#
# Usage: combine_eam.awk potential1.pt potential2.pt ....
#
# Combines all potentials into one single potential file written to 
#   standard output. Can be used to combine three potentials to form
#   an EAM potential, but will combine any files containing only 
#   a single potential.
#
# ATTENTION: combine_eam.awk is PRE-ALPHA! No validation whatsoever!!!
#
####################################################################

BEGIN { gstring="#G" }
$1=="#F" { #beginning of new file
  fn=ARGIND
  if ($2 != 3 || $3 != 1 ) {
    print "Error - wrong potential format" ;
    exit 2;
  }
  while (substr($0,1,1)=="#") { 
    getline; 
    if ($1=="#G") gstring = gstring " " $2;
  }
  range[fn]=$0;
  steps[fn]=$3;
  getline; getline;
  for (i=0;i<=steps[fn];i++) {
    table[i,fn]=$0;
    getline;
  }
}

END {
  print "#F 3 " fn;
  print gstring;
  print "## EAM potential collated with combine_eam.awk at " strftime();
  printf("## from files");
  for (i=1;i<=fn;i++) printf(" %s",ARGV[i]);
  print "";
  print "#E";
  for (i=1;i<=fn;i++) print range[i];
  for (i=1;i<=fn;i++) {
    print "";
    for (j=0;j<=steps[i];j++)
      print table[j,i];
  }
}
  
