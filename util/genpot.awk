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
# $Revision: 1.2 $
# $Date: 2008/10/09 18:03:01 $
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
#    0) f(r) = 0
#       constant potential = 0
#    1) f(r) = a1 * r^(-12) 
#       repulsive part of Lennard-Jones potential.
#    2) f(r) = a1 * exp(-a2 * r)
#       exponential decay
#    3) f(r) = a1 * (r - a[2])^2 + a[3]
#       parabola with vertex (a[2],a[3])
#
BEGIN {
  # define functions
  # Repulsive part of Lennard-Jones
  fn[0] ="f(r)=0";
  cmt[0]="Constant potential, zero everywhere."
  fn[1] ="f(r)=a[1]*r^-12";
  cmt[1]="Repulsive part of Lennard-Jones potential.";
  fn[2] ="f(r)=a[1]*exp(-a[2]*r)";
  cmt[2]="Exponential decay";
  fn[3] ="f(r)=a[1]*(r-a[2])^2 + a[3]";
  cmt[3]="Parabola with vertex (a[2],a[3])";
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
  i = 0;
  while ( i in fn ) {
    print i ") " fn[i] "\n" cmt[i] > "/dev/stderr";
    i++;
  }
  print "Which potential type do you want?" > "/dev/stderr";
}
NR==4 {
  type=$1;
  gstring=(type==3)?"3":"2";
  print "Enter space separated list of parameters" > "/dev/stderr";
}
NR>4 {
  numpar=split($0,a);
  exit;
}
END {
#header
  print "#F 3 1";
  print "#G " gstring;
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
  grad = 0;
  if (type==1) {
    grad=1.0/(r*r*r);
    grad *= grad;
    grad *= grad;
    grad *= -12.0*a[1]/r;
  } else if (type==2) {
    grad = -a[1]*a[2]*exp(-a[2]*r);
  } else if (type==3) {
    grad = 2.*a[1]*(r-a[2]);
  }
  if (type==3) {
# Not a pair, but an embedding potential!
    printf("%.10e %.10e\n", grad,2.*a[1]*(maxdist-a[2]));
  } else {
    printf("%.10e %.10e\n", grad,0.0000000);
  }
# Table
  val=0;
  for (i=0;i<nsamp-1;i++) {
    r=mindist+i*step; 
    if (type==1) {
      val = 1.0/(r*r*r); 
      val *= val; 
      val *= val; 
      val *= a[1];
    } else if (type==2) {
      val = a[1]*exp(-a[2]*r);
    } else if (type==3) {
      val = a[1]*(r-a[2])*(r-a[2])+a[3];
    }
    printf("%.10e\n", val);
  }
  if (type==3) {
# not a pair, but an embedding potential
    printf("%.10e\n", a[1]*(maxdist-a[2])*(maxdist-a[2])+a[3]);
  } else {
    print "0.000000000\n";
  }
}
	
