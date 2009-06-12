#!/bin/bash
if [ $# -eq 0 ]; then
	echo "Usage: get_config FILENAME"
	exit
elif [ $# -ge 2 ]; then
	echo "More then 1 file given on command line. Only $1 will be interpreted"
	echo
fi
n_conf=`grep -e \#E $1 | wc -l`
n_names=`grep -e generated $1 | wc -l`
if [ $n_conf -gt $n_names ]; then
grep -e generated -e \#E $1 | awk -v file=$1 'BEGIN{i=0}
{
	printf "%d ",i;
	if ($1=="#E") {
		i++;
		printf "no information found in %s\n",file;
	}
	else {
		printf "%s\n",$7;
		getline;
		i++;
	}
}'
else
grep -e generated $1 | awk 'BEGIN{i=0} {print i" "$7;i++; }'
fi
