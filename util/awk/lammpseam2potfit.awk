#!/usr/bin/awk -f
#################################################################
#
# lammpseam2potfit.awk:  convert lammps eam/alloy potential to potfit format 3
#
#################################################################
#
#   Copyright 2011 Peter Brommer
#             Departement de physique et RQMP
#             Universite de Montreal, Montreal, QC, Canada and
#             ITAP, Universitaet Stuttgart, Stuttgart, Germany
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
# Usage: lammpseam2potfit.awk <potential>
#
#        <potential> is a lammps eam/alloy potential as described
#        in http://lammps.sandia.gov/doc/pair_eam.html.
#        Command writes an potfit tabulated EAM potential (type 3)
#        to standard output. Atom types are sorted alphabetically.   
# 
#        ATTENTION: USE AT YOUR OWN RISK. VERIFY POTENTIALS!
#
#################################################################


# First 3 lines: Comments
# Line 4..5: Metadata
NR<=3 { print "##",$0 } 
NR==4 { 
    ntypes=$1
    ntypepairs=(ntypes*(ntypes+1))/2
    for (i=1;i<=ntypes;i++) {
	types[i]=$(i+1)
	order[$(i+1)] = i
    }
    asorti(order,sorttypes)
    # now types contains list of Element names as in LAMMPS file
    # order assigns Element names the original Lammps index
    # sorttypes is the list of Element names as in potfit

}
NR==5 { 
    nrho=$1
    drho=$2
    nr=$3
    dr=$4
    cutoff=$5
    #HEADER
    print "#F 3 " (ntypes*(ntypes+5))/2
    print "#T EAM"
    printf "#C"
    for (i=1;i<=ntypes;i++) printf " %s",sorttypes[i]
    print "\n#E"
    #Heads
    # Pair potentials
    for (i=1;i<=ntypepairs;i++) {
	printf "%17.10f %17.10f %10d\n",dr,(nr-1.)*dr,nr-1
    }
    # Transfer functions:
    for (i=1;i<=ntypes;i++) printf "%17.10f %17.10f %10d\n",0.0,(nr-1.)*dr,nr
    # Embedding function
    for (i=1;i<=ntypes;i++) printf "%17.10f %17.10f %10d\n",0.0,(nrho-1.)*drho,nrho
    print ""
# Now let's read the data
    for (i=1;i<=ntypes;i++) {
	getline # header line of Element, ignore
	nval=0
	while (nval < nrho) {
	    getline # Embed
	    for (k=1;k<=NF;k++) embed[i,nval+k]=$k
	    nval+=NF
	}
	nval=0
	while (nval < nr ) {
	    getline # Transfer
	    for (k=1;k<=NF;k++) transf[i,nval+k]=$k
	    nval+=NF
	}
    }
    for (i=1;i<=ntypes;i++) {
	for (j=1;j<=i;j++) {
	    nval=0
	    while (nval < nr) {
		getline # Pair
		for (k=1;k<=NF;k++) {
		    # ignore automatic zero
		    r=(nval+k-1)*dr
		    if (nval+k>1) pair[j,i,nval+k-1]=$k/r
		}
		nval+=NF
	    }
	}
    }
    # and write out:
    for (i=1;i<=ntypes;i++) {
	for (j=i;j<=ntypes;j++) {
	    for (k=1;k<=nr-1;k++) {
		m=order[sorttypes[i]]
		n=order[sorttypes[j]]
		if (n<m) { t=m;m=n;n=t }
		printf "%17.10e\n", pair[m,n,k]
	    }
	    print ""
	}
    }
    for (i=1;i<=ntypes;i++) {
	for (k=1;k<=nr;k++) 
	    printf "%17.10e\n", transf[order[sorttypes[i]],k]
	print ""
    }
    for (i=1;i<=ntypes;i++) {
	for (k=1;k<=nr;k++) 
	    printf "%17.10e\n", embed[order[sorttypes[i]],k]
	print ""
    }
}

