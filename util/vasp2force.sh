#!/bin/sh
####################################################################
#
#   Copyright 2003-2009 Peter Brommer, Daniel Schopf
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
#* $Revision: 1.17 $
#* $Date: 2009/05/12 08:31:03 $
#*****************************************************************/

wdir=`pwd`
list_types="0";
new_list="0";
f_arg="0";
s_arg="0";
poscar=`[ -f POSCAR ]`;
e_file="../single_atom_energies";
saeng="0 0 0";
mycat="cat";
declare -a type_list;

while getopts 'e:lfr:s:?h' OPTION
do
    case $OPTION in
	e) e_file="$OPTARG";
            ;;
	l) list_types="1";
      	    ;;
	f) f_arg="1";
	    ;;
	s) s_arg="$OPTARG";
	    ;;
	r) new_list="$OPTARG";
            ;;
	?) printf "\nUsage: %s: [-e file] [-f] [-l] [-r list] [-s list] <OUTCAR files>\n" $(basename $0) >&2
            printf "\n <OUTCAR files> is an optional list of files, if not given" >&2
	    printf "\n all files starting with OUTCAR will be scanned" >&2
	    printf "\n (it is possible to read gzipped files ending with .gz)\n" >&2
            printf "\n -e <file>\t\tspecify file to read single atom energies from\n" >&2
	    printf "\t\t\tif not found, \"0\" will be used for every atom type\n" >&2
	    printf " -f\t\t\tonly use the final configuration from OUTCAR\n" >&2
	    printf " -l\t\t\tlist all chemical species found in OUTCAR and exit\n" >&2
	    printf " -r <list>\t\trewrite atom types according to <list>\n" >&2
	    printf "\t\t\t(same order as listed with -l; -r 1,0 will swap 0 and 1)\n" >&2
	    printf " -s <list>\t\tcomma separated list of configurations to use\n" >&2
	    exit 2
	    ;;
    esac
done

shift $(($OPTIND-1))
outcars="$*";

if [ "$outcars" == "" ]; then
    echo "No files specified on the command line. Searching for OUTCAR* ..." >&2;
    for i in `find . -maxdepth 1 -name OUTCAR\* -print`; do
	outcars="${outcars} `basename $i`";
    done
    if [ "$outcars" != "" ]; then
	echo "Found the following files: $outcars" >&2;
    else
	echo "Could not find any OUTCAR files in this directory, aborting."
	break;
    fi
fi

for file in $outcars; do
    if [ ! -f $file ]; then
	continue;
    fi
    mycat="cat";
    zipped=`echo $file | awk 'sub(/\.gz$/,""){print $0}'`
    if [ "$zipped" != "" ]; then
	mycat="zcat";
    fi
    count=`$mycat $file | grep -c "TOTAL-FORCE"`;
    pr_conf="0";
    if [ "$f_arg" == "1" ]; then
	pr_conf="${pr_conf},$count";
    fi
    if [ "$s_arg" != "0" ]; then
	pr_conf="${pr_conf},$s_arg";
    fi
    if [ "$list_types" == "1" ]; then
	types=`$mycat $file | grep VRHFIN | wc -l`;
	if [ $types != "0" ]; then
	    name=(`$mycat $file | grep VRHFIN | awk '{ sub("=",""); sub(":",""); print $2; }'`);
	    echo "Found $types atom types in $file:";
	    for (( i=0; $i<$types; i++ )); do
		echo ${name[$i]} "= "$i;
	    done
	else
	    types=`$mycat $file | grep TITEL | wc -l`;
	    if [ $types != "0" ]; then
		name=(`$mycat $file | grep TITEL | awk '{print $4; }'`);
		echo "Found $types atom types in $file:";
		for (( i=0; $i<$types; i++ )); do
		    echo ${name[$i]} "= "$i;
		done
	    else
		types=`$mycat $file | grep POTCAR | wc -l`;
		name=(`$mycat $file | grep POTCAR | awk '{print $3; }'`);
		echo "Found $(($types/2)) atom types in $file:";
		for (( i=0; $i<$(($types/2)); i++ )); do
		    echo ${name[$i]} "= "$i;
		done
	    fi
	fi
    else
	echo "There are $count configurations in $file" >&2
	if [ -f $e_file ]; then
	    saeng=`cat $e_file`;
	else
	    printf "Warning: $e_file could not be found - using \"0 0 0\"\n" >&2;
	fi
	
	if [ "X$pr_conf" == "X0" ]; then
	    pr_conf=1;
	    for (( i=2; $i<=$count; i++ )); do
		pr_conf="${pr_conf},$i";
	    done
	fi
	
	$mycat $file | awk -v pr_conf="${pr_conf}" -v wdir="${wdir}" -v poscar="${poscar}" -v new_list="$new_list" -v saeng="$saeng" '  BEGIN {
    OFMT="%11.7g"
#Select confs to print
    count=0;
    split(pr_conf,pr_arr,",");
    for (i in pr_arr) pr_flag[pr_arr[i]]++;
#pr_flag now is set for the configurations to be printed.
#    getline saeng < "../single_atom_energies";
    if (poscar) {
    getline < "POSCAR"; getline scale < "POSCAR";
    getline boxx < "POSCAR"; getline boxy < "POSCAR"; getline boxz < "POSCAR";
    getline < "POSCAR"; ntypes = split($0,a);
    }
    single_energy=0.;
    split(saeng,sae);
    if (new_list!="0") {
	    split(new_list,number,",");
    } else {
	number[1]=0;number[2]=1;number[3]=2;
    }
  };
  /ions per type/ {
    ntypes = split($0,a) - 4;
    for (i=1;i<=ntypes;i++) a[i]=a[i+4];
    for (i=1; i<=ntypes; i++) single_energy += a[i]*sae[i];
    for (i=2; i<=ntypes; i++) a[i]+=a[i-1];
    j=1; for (i=1; i<=a[ntypes]; i++) { if (i>a[j]) j++; b[i]=j-1; }
    }
# Marker for energy is two spaces between energy and without.
# Correct energy is energy(sigma->0)
  /energy  without/ {
    energy=($7-single_energy)/a[ntypes];
    delete stress;
  }
# Find box vectors
  /VOLUME and BASIS/ {
    for (i=1;i<=4;i++) getline;
    getline boxx; getline boxy; getline boxz;
    gsub("-"," -",boxx);
    gsub("-"," -",boxy);
    gsub("-"," -",boxz);
    split(boxx,boxx_v);
    split(boxy,boxy_v);
    split(boxz,boxz_v);
    scale=1.0;
}
  ($2=="kB") {
     for (i=1;i<=6;i++) stress[i]=$(i+2)/1602.;
}
  ($2=="TOTAL-FORCE") {
     count++;
     if (count in pr_flag) {
       print "#N",a[ntypes],1; #flag indicates whether to use forces or not
       print "## force file generated from directory " wdir;
       printf "#X %13.8f %13.8f %13.8f\n",boxx_v[1]*scale,\
                   boxx_v[2]*scale,boxx_v[3]*scale;
       printf "#Y %13.8f %13.8f %13.8f\n",boxy_v[1]*scale,\
                   boxy_v[2]*scale,boxy_v[3]*scale;
       printf "#Z %13.8f %13.8f %13.8f\n",boxz_v[1]*scale,\
                   boxz_v[2]*scale,boxz_v[3]*scale;
       printf("#E %.10f\n",energy) ;
       if ( 1 in stress )
         print "#S",stress[1],stress[2],stress[3],stress[4],stress[5],stress[6];
       print "#F";
     }
     getline; getline;
     for (i=1; i<=a[ntypes]; i++) {
       if (count in pr_flag)
	   printf("%d %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g\n",
	       number[b[i]+1],$1,$2,$3,$4,$5,$6);
     getline; }
  };'
    fi
done
