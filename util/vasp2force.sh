#!/bin/sh
#/****************************************************************
#* $Revision: 1.2 $
#* $Date: 2004/03/23 10:07:00 $
#*****************************************************************/

[ -f ../single_atom_energies ] || { echo file ../single_atom_energies not found; exit;}
cat OUTCAR | awk '
  BEGIN { getline saeng < "../single_atom_energies"; 
    getline < "POSCAR"; getline scale < "POSCAR";
    getline boxx < "POSCAR"; getline boxy < "POSCAR"; getline boxz < "POSCAR";
    getline < "POSCAR"; ntypes = split($0,a);single_energy=0.;
    split(saeng,sae);
    #sae[1]=-0.000219; sae[2]=-0.993872; sae[3]=-0.855835;
    for (i=1; i<=ntypes; i++) single_energy += a[i]*sae[i];
     for (i=2; i<=ntypes; i++) a[i]=a[i-1]+$i;
    j=1; for (i=1; i<=a[ntypes]; i++) { if (i>a[j]) j++; b[i]=j-1; }
    split(boxx,boxx_v);
    split(boxy,boxy_v);
    split(boxz,boxz_v);
    print a[ntypes]; 
    print boxx_v[1]*scale " " boxx_v[2]*scale " " boxx_v[3]*scale; 
    print boxy_v[1]*scale " " boxy_v[2]*scale " " boxy_v[3]*scale; 
    print boxz_v[1]*scale " " boxz_v[2]*scale " " boxz_v[3]*scale; 
  }
  ($2=="ENERGIE") {getline; getline; getline; getline;
     print ($4-single_energy)/a[ntypes] }
  ($2=="kB") {print $3/1602., $4/1602., $5/1602., $6/1602., $7/1602., $8/1602.}
  ($2=="TOTAL-FORCE") { getline; getline;
     for (i=1; i<=a[ntypes]; i++) { print b[i],$1,$2,$3,$4,$5,$6; getline; } 
  }'
