#!/usr/bin/awk -f
####################################################################
#
# genpot.awk - generates an interpolated representation of an analytic 
#              potential with an arbitrary number of parameters.
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
# $Date: 2008/10/09 10:08:36 $
####################################################################
#
# Usage: genpot.awk 
#    takes no parameters, however the parameters of the potential 
#    may be piped into the script (each parameter on one line)
#    mindist, maxdist, nsamp, type, space separated list of params.
#
# ATTENTION: genpot.awk is pre-alpha. There is no verification of
#    ANY parameter.
#
###################################################################
#
# Implemented potential types:
#
#    1) f(r) = a1 * r^(-12) 
#       repulsive part of Lennard-Jones potential.
#
BEGIN {
  # define functions
  # Repulsive part of Lennard-Jones
  fn[1] ="f(r)=a[1]*r^-12";

  OFMT="%.12g";
  par1=3;
  print "Enter minimal distance" > "/dev/stderr";
}
NR==1 {
  mindist=$1;
  print "Enter maximal distance" > "/dev/stderr";
}
NR==2 {
  maxdist=$1;
  print "Enter number of sampling points"  > "/dev/stderr";
}
NR==3 {
  nsamp=$1;
  i = 1;
  while ( i in fn ) {
    print i ") " fn[i] "\n" > "/dev/stderr";
    i++;
  }
  print "Which potential type do you want?" > "/dev/stderr";
}
NR==4 {
  type=$1;
  print "Enter space separated list of parameters" > "/dev/stderr";
}
NR>4 {
  numpar=split($0,a);
  exit;
}
END {
#header
  print "#F 3 1";
  print "#G 2";
  print "## Interpolated potential from  " fn[type]  ;
  print "## generated with genpot.awk at " strftime();
  print "## Parameters:";
  for (i=1;i<=numpar;i++) {
    print "## a[" i "]= " a[i];
  }
  print "#E" ;

# Header line
  printf("%.10f %.10f %i\n\n", mindist,maxdist,nsamp);
  step=(maxdist-mindist)/(nsamp-1.0); 

# Gradient
  r=mindist;
  if (type==1) {
    grad=1.0/(r*r*r);
    grad *= grad;
    grad *= grad;
    grad *= -12.0*a[1]/r;
  }
  printf("%.10e %.10e\n", grad,0.0000000);
     
# Table
  for (i=0;i<nsamp-1;i++) {
    r=mindist+i*step; 
    if (type==1) {
      val = 1.0/(r*r*r); 
      val *= val; 
      val *= val; 
      val *= a[1];
    }
    printf("%.10e\n", val);
  }
  print "0.000000000\n";
}
	
