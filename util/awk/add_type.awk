#!/usr/bin/awk -f
#################################################################
#
# add_type.awk: add one more atom type to EAM potential.
#
#################################################################
# 
#   Copyright 2008 Peter Brommer
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
#################################################################
#
# Usage: add_type.awk orig_eam.pt add_eam.pt cross_pair1.pt ...
#
# add_type.awk needs the following arguments
#   orig_eam.pt    - original EAM potential for n atom types.
#   add_eam.pt     - EAM potential for additional atom type (3 columns).
#   cross_pair*.pt - old-new pair interaction potential, each 1 column.
#                    (1 file for each original atom type, n in total).
#
# The resulting potential is written to standard output.
#
# ATTENTION add_type.awk is PRE-ALPHA. No valdiation whatsoever!!!
# 
#################################################################


ARGIND==1 { # Original potential 
  while (substr($0,2,1)!="F") getline;
  if ($2 != 3 ) {
    print "Error - wrong potential format of " ARGV[ARGIND]  ;
    exit 2;
  }
  ncols_old=$3;
  ntypes=(-5+sqrt(25+8*ncols_old))/2;
  if (int(ntypes)!=ntypes) { 
    printf "ERROR - incorrect parameter file " ARGV[ARGIND]  ;
    exit 2;
  }
  while (substr($0,1,2)!="#E") { 
    getline; 
    if ($1=="#I") {
      for (i=1;i<=ncols_old;i++) {
	istring[i] = $(i+1);
      }
    }
    if ($1=="#G") {
      for (i=1;i<=ncols_old;i++) {
	gstring[i] = $(i+1);
      }
    }
  }
  for (i=1;i<=ncols_old;i++) {
    getline;
    range[i]=$0;
    steps[i]=$3;
  }
  for (j=1;j<=ncols_old;j++) {
    getline; 			# blank line between tables
    for (i=0;i<=steps[j];i++) {
      getline;
      table[i,j]=$0;
    }
  }  
  nextfile;
}

ARGIND==2 { # EAM potential of additional atom 
  while (substr($0,2,1)!="F") getline;
  if ($2 != 3 || $3 != 3 ) {
    print "Error - wrong potential format of " ARGV[ARGIND]  ;
    exit 2;
  }
  startcol=ncols_old+1;
  while (substr($0,1,2)!="#E") { 
    getline; 
    if ($1=="#G") {
      for (i=0;i<3;i++) {
	gstring[i+startcol] = $(i+2);
      }
    }
    if ($1=="#I") {
      for (i=0;i<3;i++) {
	istring[i+startcol] = $(i+2);
      }
    }
  }
  for (i=startcol;i<startcol+3;i++) {
    getline;
    range[i]=$0;
    steps[i]=$3;
  }
  for (j=startcol;j<startcol+3;j++) {
    getline;
    for (i=0;i<=steps[j];i++) {
      getline;
      table[i,j]=$0;
    }
  }
  startcol += 3;
  nextfile;
}  
ARGIND>2 { # new mixed pair potentials
  while (substr($0,2,1)!="F") getline;
  if ($2 != 3 || $3 != 1 ) {
    print "Error - wrong potential format of " ARGV[ARGIND]  ;
    exit 2;
  }  
  while (substr($0,1,2)!="#E") {
    getline;
    if ($1=="#G") {
      gstring[startcol] = $2;
    }
    if ($1=="#I") {
      istring[startcol] = $2;
    }
  }
  getline;
  range[startcol]=$0;
  steps[startcol]=$3;
  getline;
  for (i=0;i<=steps[startcol];i++) {
    getline;
    table[i,startcol]=$0;
  }
  startcol += 1;
  nextfile;
}

END {
# New order
  ntypes_new = ntypes+1;
  ncols = (ntypes_new * (ntypes_new + 5) )/ 2;
  if (ncols != startcol-1 ) {
    print "Error - not enough/too many columns read";
    exit;
  }
  # Pair potentials
  new_cols=1;
  old_cols=1;
  for (i=1;i<=ntypes;i++) {
    for (j=i; j<=ntypes; j++) {
      nc[new_cols++]=old_cols++;
    }
    # insert new column here
    nc[new_cols++]=ncols_old+3+i;
  }
  nc[new_cols++]=ncols_old+1;
  # Transfer fns
  for (i=1;i<=ntypes; i++) {
    nc[new_cols++]=old_cols++;
  }
  nc[new_cols++]=ncols_old+2;
  # Embedding fns
  for (i=1;i<=ntypes; i++) {
    nc[new_cols++]=old_cols++;
  }
  nc[new_cols]=ncols_old+3;
  print "#F 3 " new_cols;
  printf ("#G");
  for (i=1;i<=new_cols;i++) {
    printf(" %i",gstring[nc[i]]);
  }
  print "";
  printf ("#I");
  for (i=1;i<=new_cols;i++) {
    if (nc[i] in istring) {
      printf(" %i",istring[nc[i]]);
    } else {
      printf(" 0");
    }
  }
  print "";

  print "## EAM potential " ARGV[1] " with additional type ";
  print "## from " ARGV[2] ", cross-pair potentials used: ";
  printf("##");
  for (i=1;i<=ntypes;i++) printf(" %s",ARGV[i+2]);
  print "";
  print "#E";
  for (i=1;i<=new_cols;i++) print range[nc[i]];
  for (i=1;i<=new_cols;i++) {
    print "";
    for (j=0;j<=steps[nc[i]];j++)
      print table[j,nc[i]];
  }
}
  
