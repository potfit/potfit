#!/usr/bin/awk -f
#################################################################
#
# compnodes:    calculates compositions of configurations
# 		and proposes some composition nodes
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
#################################################################
#
# Usage: compnodes <option> config.file
#
# <option> can be one of the following:
# 	-l 	lists all configurations with different atom types
# 	-la 	lists all configurations
# 	-ls 	lists all configurations with different atom types sorted
# 	-le 	lists all effective configurations
#
#################################################################

function reduce_nodes(n,nodes,weights,new_nodes)
{
    for (i=0;i<n;i++) {
	nodes[i]=(nodes[i]*weights[i]+nodes[i+1]*weights[i+1])/(weights[i]+weights[i+1]);
	weights[i]+=weights[i+1];
    }
}

function sort(to_sort,to_change,n,nr)
{
    false=0
    true=1
    readyflag=false;
    
    j=0;
    while (readyflag==false)
    {
	readyflag=true;
	for (i=1;i<n;i++) {
	    if ((to_sort[i-1]/to_change[i-1])>(to_sort[i]/to_change[i]))
	    {
		temp=to_sort[i-1];
		temp2=to_change[i-1];
		to_sort[i-1]=to_sort[i];
		to_change[i-1]=to_change[i];
		to_sort[i]=temp;
		to_change[i]=temp2;
		readyflag=false;
	    }
	    
	}
    }
}


function list_conf(ntot,conc,n,conf) {
    if (listall)
	printf "Listing all "n" configurations with 2 atom types:\n";
    if (listsorted)
	printf "Listing "n" sorted configurations with 2 atom types:\n";
    if (list)
	printf "Listing "n" configurations with 2 atom types:\n";
    if (listeff)
	printf "Listing "n" effective configurations with 2 atom types:\n";
    printf "conf\tc[0]\tc[1]\tn[0]\tn[1]\n";
    for (i=0;i<n;i++)
	printf "%d\t%.3f\t%.3f\t%d\t%d\n",conf[i]+1,conc[conf[i]]/ntot[conf[i]],1-conc[conf[i]]/ntot[conf[i]],conc[conf[i]],ntot[conf[i]]-conc[conf[i]];
}

BEGIN {
    n=0; 			# number of configuration
    max=0; 			# ntypes
    if (ARGC==1) {
	print "Configuration file missing." > "/dev/stderr";
	abort=1;
	exit 1;
    }
    
    if (ARGV[1]=="-l") {
	list=1;
	delete ARGV[1];
    } else if (ARGV[1]=="-la") {
	listall=1;
	delete ARGV[1];
    } else if (ARGV[1]=="-ls") {
	listsorted=1;
	delete ARGV[1];
    } else if (ARGV[1]=="-le") {
	listeff=1;
	delete ARGV[1];
    } else if (ARGV[1]=="-n") {
	only_n=1;
	print_n=ARGV[2];
	delete ARGV[1];
	delete ARGV[2];
    } else list=0;
}

{
    if (substr($1,1,1)=="#") {
	nconf[n,"tot"]=$2;
	while (substr($1,2,1)!="F") getline;
	getline;
	for (i=0;i<(nconf[n,"tot"]-1);i++) {
	    if ($1>max) max=$1;
	    nconf[n,$1]++;
	    getline;
	}
	if ($1>max) max=$1;
	nconf[n,$1]++;
	n++;
    } else {
	nconf[n,"tot"]=$1;
	getline; getline; getline; getline; getline; getline;
	for (i=0;i<(nconf[n,"tot"]-1);i++) {
	    if ($1>max) max=$1;
	    nconf[n,$1]++;
	    getline;
	}
	if ($1>max) max=$1;
	nconf[n,$1]++;
	n++;
    }
    
}

END {
    if (abort)
	exit 1;
    if (max>=2) {
	print "More than 2 different atom types detected.";
	print "Showing only the configurations, compnodes disabled.";
	listall=1;
	list=0;
	listsorted=0;
	listeff=0;
    }
    
    if (!only_n)
	print "Found "n" configurations with "max+1" different atom types.";
    if (listall) {
	print "Listing all "n" configurations:";
	printf "conf";
	for (i=0;i<=max;i++)
	    printf "\tc["i"]";
	printf "\t";
	for (i=0;i<=max;i++)
	    printf "\tn["i"]";
    }
    for (i=0;i<n;i++) {
	if (listall)
	    printf "\n"i+1;
	for (j=0;j<=max;j++){
	    c[i,j]=nconf[i,j]/nconf[i,"tot"];	
	    if (listall)
		printf "\t%.3f",c[i,j];
	}
	printf "\t";
	for (j=0;j<=max;j++){
	    if (listall)
		printf "\t%d",nconf[i,j];
	}
    }
    if (listall)
	printf "\n";
    
    avg=0;
    count=0;
    zeroconf=0;
    j=0;
    for (i=0;i<n;i++) {
	if (c[i,0]!=0 && c[i,0]!=1) {
	    avg=avg+c[i,0]*nconf[i,"tot"];
	    count=count+nconf[i,"tot"];
	    nonzeroconf[j]=i;
	    j++;
	} else 
	    zeroconf++;
    }
    
    len=n-zeroconf;
    
    for (i=0;i<=len;i++) {
	concentration[i]=nconf[nonzeroconf[i],0];
	config[i]=nconf[nonzeroconf[i],"tot"];
    }
    
    if (list || listsorted)
	for (i=0;i<=len;i++)
	    nonzeroconf[i]=i;
    if (list)
	list_conf(config,concentration,len,nonzeroconf);
    
    sort(concentration,config,len,nonzeroconf);
    
    if (listsorted)
	list_conf(config,concentration,len,nonzeroconf);
    
# filter out multiple configurations with same concentration
    k=0;
    if (len>1) {
	for (i=0;i<(len-1);i++) {
	    j=1;
	    done=false
	    while (done==false && i+j < (len-1)) {
		done=true;
		if ((concentration[i]/config[i])==(concentration[i+j]/config[i+j])) {
		    concentration[i]+=concentration[i+j];
		    config[i]+=config[i+j];
		    done=false;
		    j++;
		} 
	    }
	    effconfig[k]=i;
	    k++;
	    i=i+j-1;
	}
	
	
	if ((concentration[len-2]/config[len-2])!=(concentration[len-1]/config[len-1])) {
	    effconfig[k]=len-1;
	    k++
	}
    } else {
	effconfig[0]=0;
	k=1;
    }
    
    if (listeff)
	list_conf(config,concentration,k,effconfig);
    
    
    for (i=0;i<k;i++) {
	concentration[i]=concentration[effconfig[i]]/config[effconfig[i]];
	config[i]=config[effconfig[i]];
    }
    
    if (!only_n && max<2) {
	print "Found "k" effective configurations.";
	if (k>5) {
	    print "All composition nodes:";
	    for (i=0;i<k;i++)
		printf "%.3f ",concentration[i];
	    printf "\n";
	}
	
	print "Proposing the following composition nodes:";
	if (k<6) {
	    printf k" node(s): ";
	    for (i=0;i<k;i++)
		printf "%.3f ",concentration[i];
	    printf "\n";
	}
    }
    if (only_n && max<2) {
	if (print_n>k) {
	    print "There are not enough configurations for "print_n" composition nodes" > "/dev/stderr";
	    exit 1;
	}
	if (print_n==k) {
	    for (i=0;i<k;i++)
		printf "%.3f ",concentration[i];
	    printf "\n";
	}
    }
    if (max<2)
    for (i=(k-1);i>=0;i--) {
	reduce_nodes(i,concentration,config);	
	show = false;
	if (!only_n){
	    prefix=i" node(s): ";
	    show=true;
	}
	if (only_n)
	    if (i==print_n)
		show=true;
	if (i==5 && show) {
	    printf "%s%.3f %.3f %.3f %.3f %.3f\n",prefix,concentration[0],concentration[1],concentration[2],concentration[3],concentration[4]; 
	}
	if (i==4 && show) {
	    printf "%s%.3f %.3f %.3f %.3f\n",prefix,concentration[0],concentration[1],concentration[2],concentration[3]; 
	}
	if (i==3 && show) {
	    printf "%s%.3f %.3f %.3f\n",prefix,concentration[0],concentration[1],concentration[2];
	}
	if (i==2 && show) {
	    printf "%s%.3f %.3f\n",prefix,concentration[0],concentration[1];
	}
	if (i==1 && show) {
	    printf "%s%.3f\n",prefix,concentration[0];
	}
    }
}
