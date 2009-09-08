#!/bin/bash
rm -rf formulas *.png;
latex2html -no_reuse formulas.tex;
declare -a names;
names=(`grep ^\%\% formulas.tex | awk '{print $2}'`);
j=0;
for i in "${names[@]}"
	do
		cp formulas/img$((2*j+1)).png ${names[$j]}.png;
		cp formulas/img$((2*j+2)).png ${names[$j]}-order.png;
		j=$((j+1));
	done
rm -rf formulas
