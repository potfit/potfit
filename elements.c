/****************************************************************
 *
 * elements.c: data and routines for periodic table
 *
 ****************************************************************
 *
 * Copyright 2012
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.itap.physik.uni-stuttgart.de/
 *
 ****************************************************************
 *
 *   This file is part of potfit.
 *
 *   potfit is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   potfit is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

#include "potfit.h"

#include "elements.h"

#define N_ELEMENTS 110

element_t *element_table;

void init_elements()
{
  element_table = (element_t *) malloc((N_ELEMENTS + 1) * sizeof(element_t));

/* 1 - Hydrogen */
  strcpy(element_table[1].name, "Hydrogen");
  strcpy(element_table[1].short_name, "H");
  element_table[1].mass = 1.0079;
/* 2 - Helium */
  strcpy(element_table[2].name, "Helium");
  strcpy(element_table[2].short_name, "He");
  element_table[2].mass = 4.0026;
/* 3 - Lithium */
  strcpy(element_table[3].name, "Lithium");
  strcpy(element_table[3].short_name, "Li");
  element_table[3].mass = 6.941000;
/* 4 - Beryllium */
  strcpy(element_table[4].name, "Beryllium");
  strcpy(element_table[4].short_name, "Be");
  element_table[4].mass = 9.012180;
/* 5 - Boron */
  strcpy(element_table[5].name, "Boron");
  strcpy(element_table[5].short_name, "B");
  element_table[5].mass = 10.810000;
/* 6 - Carbon */
  strcpy(element_table[6].name, "Carbon");
  strcpy(element_table[6].short_name, "C");
  element_table[6].mass = 12.011000;
/* 7 - Nitrogen */
  strcpy(element_table[7].name, "Nitrogen");
  strcpy(element_table[7].short_name, "N");
  element_table[7].mass = 14.006740;
/* 8 - Oxygen */
  strcpy(element_table[8].name, "Oxygen");
  strcpy(element_table[8].short_name, "O");
  element_table[8].mass = 15.999400;
/* 9 - Fluorine */
  strcpy(element_table[9].name, "Fluorine");
  strcpy(element_table[9].short_name, "F");
  element_table[9].mass = 18.998403;
/* 10 - Neon */
  strcpy(element_table[10].name, "Neon");
  strcpy(element_table[10].short_name, "Ne");
  element_table[10].mass = 20.179000;
/* 11 - Sodium */
  strcpy(element_table[11].name, "Sodium");
  strcpy(element_table[11].short_name, "Na");
  element_table[11].mass = 22.989770;
/* 12 - Magnesium */
  strcpy(element_table[12].name, "Magnesium");
  strcpy(element_table[12].short_name, "Mg");
  element_table[12].mass = 24.305000;
/* 13 - Aluminum */
  strcpy(element_table[13].name, "Aluminum");
  strcpy(element_table[13].short_name, "Al");
  element_table[13].mass = 26.981540;
/* 14 - Silicon */
  strcpy(element_table[14].name, "Silicon");
  strcpy(element_table[14].short_name, "Si");
  element_table[14].mass = 28.086000;
/* 15 - Phosphorus */
  strcpy(element_table[15].name, "Phosphorus");
  strcpy(element_table[15].short_name, "P");
  element_table[15].mass = 30.973760;
/* 16 - Sulfur */
  strcpy(element_table[16].name, "Sulfur");
  strcpy(element_table[16].short_name, "S");
  element_table[16].mass = 32.060000;
/* 17 - Chlorine */
  strcpy(element_table[17].name, "Chlorine");
  strcpy(element_table[17].short_name, "Cl");
  element_table[17].mass = 35.453000;
/* 18 - Argon */
  strcpy(element_table[18].name, "Argon");
  strcpy(element_table[18].short_name, "Ar");
  element_table[18].mass = 39.948000;
/* 19 - Potassium */
  strcpy(element_table[19].name, "Potassium");
  strcpy(element_table[19].short_name, "K");
  element_table[19].mass = 39.098000;
/* 20 - Calcium */
  strcpy(element_table[20].name, "Calcium");
  strcpy(element_table[20].short_name, "Ca");
  element_table[20].mass = 40.080000;
/* 21 - Scandium */
  strcpy(element_table[21].name, "Scandium");
  strcpy(element_table[21].short_name, "Sc");
  element_table[21].mass = 44.955900;
/* 22 - Titanium */
  strcpy(element_table[22].name, "Titanium");
  strcpy(element_table[22].short_name, "Ti");
  element_table[22].mass = 47.900000;
/* 23 - Vanadium */
  strcpy(element_table[23].name, "Vanadium");
  strcpy(element_table[23].short_name, "V");
  element_table[23].mass = 50.941400;
/* 24 - Cromium */
  strcpy(element_table[24].name, "Cromium");
  strcpy(element_table[24].short_name, "Cr");
  element_table[24].mass = 51.996000;
/* 25 - Manganese */
  strcpy(element_table[25].name, "Manganese");
  strcpy(element_table[25].short_name, "Mn");
  element_table[25].mass = 54.938000;
/* 26 - Iron */
  strcpy(element_table[26].name, "Iron");
  strcpy(element_table[26].short_name, "Fe");
  element_table[26].mass = 55.847000;
/* 27 - Cobalt */
  strcpy(element_table[27].name, "Cobalt");
  strcpy(element_table[27].short_name, "Co");
  element_table[27].mass = 58.933200;
/* 28 - Nickel */
  strcpy(element_table[28].name, "Nickel");
  strcpy(element_table[28].short_name, "Ni");
  element_table[28].mass = 58.700000;
/* 29 - Copper */
  strcpy(element_table[29].name, "Copper");
  strcpy(element_table[29].short_name, "Cu");
  element_table[29].mass = 63.546000;
/* 30 - Zinc */
  strcpy(element_table[30].name, "Zinc");
  strcpy(element_table[30].short_name, "Zn");
  element_table[30].mass = 65.380000;
/* 31 - Gallium */
  strcpy(element_table[31].name, "Gallium");
  strcpy(element_table[31].short_name, "Ga");
  element_table[31].mass = 69.720000;
/* 32 - Germanium */
  strcpy(element_table[32].name, "Germanium");
  strcpy(element_table[32].short_name, "Ge");
  element_table[32].mass = 72.590000;
/* 33 - Arsenic */
  strcpy(element_table[33].name, "Arsenic");
  strcpy(element_table[33].short_name, "As");
  element_table[33].mass = 74.921600;
/* 34 - Selenium */
  strcpy(element_table[34].name, "Selenium");
  strcpy(element_table[34].short_name, "Se");
  element_table[34].mass = 78.960000;
/* 35 - Bromine */
  strcpy(element_table[35].name, "Bromine");
  strcpy(element_table[35].short_name, "Br");
  element_table[35].mass = 79.904000;
/* 36 - Krypton */
  strcpy(element_table[36].name, "Krypton");
  strcpy(element_table[36].short_name, "Kr");
  element_table[36].mass = 83.800000;
/* 37 - Rubidium */
  strcpy(element_table[37].name, "Rubidium");
  strcpy(element_table[37].short_name, "Rb");
  element_table[37].mass = 85.467800;
/* 38 - Strontium */
  strcpy(element_table[38].name, "Strontium");
  strcpy(element_table[38].short_name, "Sr");
  element_table[38].mass = 87.620000;
/* 39 - Yttrium */
  strcpy(element_table[39].name, "Yttrium");
  strcpy(element_table[39].short_name, "Y");
  element_table[39].mass = 88.905900;
/* 40 - Zirconium */
  strcpy(element_table[40].name, "Zirconium");
  strcpy(element_table[40].short_name, "Zr");
  element_table[40].mass = 91.220000;
/* 41 - Niobium */
  strcpy(element_table[41].name, "Niobium");
  strcpy(element_table[41].short_name, "Nb");
  element_table[41].mass = 92.906400;
/* 42 - Molybdenum */
  strcpy(element_table[42].name, "Molybdenum");
  strcpy(element_table[42].short_name, "Mo");
  element_table[42].mass = 95.940000;
/* 43 - Technetium */
  strcpy(element_table[43].name, "Technetium");
  strcpy(element_table[43].short_name, "Tc");
  element_table[43].mass = 97.000000;
/* 44 - Ruthenium */
  strcpy(element_table[44].name, "Ruthenium");
  strcpy(element_table[44].short_name, "Ru");
  element_table[44].mass = 101.070000;
/* 45 - Rhodium */
  strcpy(element_table[45].name, "Rhodium");
  strcpy(element_table[45].short_name, "Rh");
  element_table[45].mass = 102.905500;
/* 46 - Palladium */
  strcpy(element_table[46].name, "Palladium");
  strcpy(element_table[46].short_name, "Pd");
  element_table[46].mass = 106.400000;
/* 47 - Silver */
  strcpy(element_table[47].name, "Silver");
  strcpy(element_table[47].short_name, "Ag");
  element_table[47].mass = 107.868000;
/* 48 - Cadmium */
  strcpy(element_table[48].name, "Cadmium");
  strcpy(element_table[48].short_name, "Cd");
  element_table[48].mass = 112.400000;
/* 49 - Indium */
  strcpy(element_table[49].name, "Indium");
  strcpy(element_table[49].short_name, "In");
  element_table[49].mass = 114.820000;
/* 50 - Tin */
  strcpy(element_table[50].name, "Tin");
  strcpy(element_table[50].short_name, "Sn");
  element_table[50].mass = 118.690000;
/* 51 - Antimony */
  strcpy(element_table[51].name, "Antimony");
  strcpy(element_table[51].short_name, "Sb");
  element_table[51].mass = 121.750000;
/* 52 - Tellurium */
  strcpy(element_table[52].name, "Tellurium");
  strcpy(element_table[52].short_name, "Te");
  element_table[52].mass = 127.600000;
/* 53 - Iodine */
  strcpy(element_table[53].name, "Iodine");
  strcpy(element_table[53].short_name, "I");
  element_table[53].mass = 126.904500;
/* 54 - Xenon */
  strcpy(element_table[54].name, "Xenon");
  strcpy(element_table[54].short_name, "Xe");
  element_table[54].mass = 131.300000;
/* 55 - Cesium */
  strcpy(element_table[55].name, "Cesium");
  strcpy(element_table[55].short_name, "Cs");
  element_table[55].mass = 132.905400;
/* 56 - Barium */
  strcpy(element_table[56].name, "Barium");
  strcpy(element_table[56].short_name, "Ba");
  element_table[56].mass = 137.340000;
/* 57 - Lanthanum */
  strcpy(element_table[57].name, "Lanthanum");
  strcpy(element_table[57].short_name, "La");
  element_table[57].mass = 138.905500;
/* 58 - Cerium */
  strcpy(element_table[58].name, "Cerium");
  strcpy(element_table[58].short_name, "Ce");
  element_table[58].mass = 140.120000;
/* 59 - Praseodynium */
  strcpy(element_table[59].name, "Praseodynium");
  strcpy(element_table[59].short_name, "Pr");
  element_table[59].mass = 140.907700;
/* 60 - Neodymium */
  strcpy(element_table[60].name, "Neodymium");
  strcpy(element_table[60].short_name, "Nd");
  element_table[60].mass = 144.240000;
/* 61 - Promethium */
  strcpy(element_table[61].name, "Promethium");
  strcpy(element_table[61].short_name, "Pm");
  element_table[61].mass = 145.000000;
/* 62 - Samarium */
  strcpy(element_table[62].name, "Samarium");
  strcpy(element_table[62].short_name, "Sm");
  element_table[62].mass = 150.400000;
/* 63 - Europium */
  strcpy(element_table[63].name, "Europium");
  strcpy(element_table[63].short_name, "Eu");
  element_table[63].mass = 151.960000;
/* 64 - Gadolinium */
  strcpy(element_table[64].name, "Gadolinium");
  strcpy(element_table[64].short_name, "Gd");
  element_table[64].mass = 157.250000;
/* 65 - Terbium */
  strcpy(element_table[65].name, "Terbium");
  strcpy(element_table[65].short_name, "Tb");
  element_table[65].mass = 158.925400;
/* 66 - Dysprosium */
  strcpy(element_table[66].name, "Dysprosium");
  strcpy(element_table[66].short_name, "Dy");
  element_table[66].mass = 162.500000;
/* 67 - Holmium */
  strcpy(element_table[67].name, "Holmium");
  strcpy(element_table[67].short_name, "Ho");
  element_table[67].mass = 164.930400;
/* 68 - Erbium */
  strcpy(element_table[68].name, "Erbium");
  strcpy(element_table[68].short_name, "Er");
  element_table[68].mass = 167.260000;
/* 69 - Thulium */
  strcpy(element_table[69].name, "Thulium");
  strcpy(element_table[69].short_name, "Tm");
  element_table[69].mass = 168.934200;
/* 70 - Ytterbium */
  strcpy(element_table[70].name, "Ytterbium");
  strcpy(element_table[70].short_name, "Yb");
  element_table[70].mass = 173.040000;
/* 71 - Lutetium */
  strcpy(element_table[71].name, "Lutetium");
  strcpy(element_table[71].short_name, "Lu");
  element_table[71].mass = 174.970000;
/* 72 - Hafnium */
  strcpy(element_table[72].name, "Hafnium");
  strcpy(element_table[72].short_name, "Hf");
  element_table[72].mass = 178.490000;
/* 73 - Tantalum */
  strcpy(element_table[73].name, "Tantalum");
  strcpy(element_table[73].short_name, "Ta");
  element_table[73].mass = 180.947900;
/* 74 - Tungsten */
  strcpy(element_table[74].name, "Tungsten");
  strcpy(element_table[74].short_name, "W");
  element_table[74].mass = 183.500000;
/* 75 - Rhenium */
  strcpy(element_table[75].name, "Rhenium");
  strcpy(element_table[75].short_name, "Re");
  element_table[75].mass = 186.207000;
/* 76 - Osmium */
  strcpy(element_table[76].name, "Osmium");
  strcpy(element_table[76].short_name, "Os");
  element_table[76].mass = 190.200000;
/* 77 - Iridium */
  strcpy(element_table[77].name, "Iridium");
  strcpy(element_table[77].short_name, "Ir");
  element_table[77].mass = 192.220000;
/* 78 - Platinum */
  strcpy(element_table[78].name, "Platinum");
  strcpy(element_table[78].short_name, "Pt");
  element_table[78].mass = 195.090000;
/* 79 - Gold */
  strcpy(element_table[79].name, "Gold");
  strcpy(element_table[79].short_name, "Au");
  element_table[79].mass = 196.966500;
/* 80 - Mercury */
  strcpy(element_table[80].name, "Mercury");
  strcpy(element_table[80].short_name, "Hg");
  element_table[80].mass = 200.590000;
/* 81 - Thallium */
  strcpy(element_table[81].name, "Thallium");
  strcpy(element_table[81].short_name, "Tl");
  element_table[81].mass = 204.370000;
/* 82 - Lead */
  strcpy(element_table[82].name, "Lead");
  strcpy(element_table[82].short_name, "Pb");
  element_table[82].mass = 207.200000;
/* 83 - Bismuth */
  strcpy(element_table[83].name, "Bismuth");
  strcpy(element_table[83].short_name, "Bi");
  element_table[83].mass = 208.980400;
/* 84 - Polonium */
  strcpy(element_table[84].name, "Polonium");
  strcpy(element_table[84].short_name, "Po");
  element_table[84].mass = 209.000000;
/* 85 - Astatine */
  strcpy(element_table[85].name, "Astatine");
  strcpy(element_table[85].short_name, "At");
  element_table[85].mass = 210.000000;
/* 86 - Radon */
  strcpy(element_table[86].name, "Radon");
  strcpy(element_table[86].short_name, "Rn");
  element_table[86].mass = 222.000000;
/* 87 - Francium */
  strcpy(element_table[87].name, "Francium");
  strcpy(element_table[87].short_name, "Fr");
  element_table[87].mass = 223.000000;
/* 88 - Radium */
  strcpy(element_table[88].name, "Radium");
  strcpy(element_table[88].short_name, "Ra");
  element_table[88].mass = 226.025400;
/* 89 - Actinium */
  strcpy(element_table[89].name, "Actinium");
  strcpy(element_table[89].short_name, "Ac");
  element_table[89].mass = 227.000000;
/* 90 - Thorium */
  strcpy(element_table[90].name, "Thorium");
  strcpy(element_table[90].short_name, "Th");
  element_table[90].mass = 232.038100;
/* 91 - Protactinium */
  strcpy(element_table[91].name, "Protactinium");
  strcpy(element_table[91].short_name, "Pa");
  element_table[91].mass = 231.035900;
/* 92 - Uranium */
  strcpy(element_table[92].name, "Uranium");
  strcpy(element_table[92].short_name, "U");
  element_table[92].mass = 238.029000;
/* 93 - Neptunium */
  strcpy(element_table[93].name, "Neptunium");
  strcpy(element_table[93].short_name, "Np");
  element_table[93].mass = 237.048200;
/* 94 - Plutonium */
  strcpy(element_table[94].name, "Plutonium");
  strcpy(element_table[94].short_name, "Pu");
  element_table[94].mass = 244.000000;
/* 95 - Americium */
  strcpy(element_table[95].name, "Americium");
  strcpy(element_table[95].short_name, "Am");
  element_table[95].mass = 243.000000;
/* 96 - Curium */
  strcpy(element_table[96].name, "Curium");
  strcpy(element_table[96].short_name, "Cm");
  element_table[96].mass = 247.000000;
/* 97 - Berkelium */
  strcpy(element_table[97].name, "Berkelium");
  strcpy(element_table[97].short_name, "Bk");
  element_table[97].mass = 247.000000;
/* 98 - Californium */
  strcpy(element_table[98].name, "Californium");
  strcpy(element_table[98].short_name, "Cf");
  element_table[98].mass = 251.000000;
/* 99 - Einsteinium */
  strcpy(element_table[99].name, "Einsteinium");
  strcpy(element_table[99].short_name, "Es");
  element_table[99].mass = 252.000000;
/* 100 - Fermium */
  strcpy(element_table[100].name, "Fermium");
  strcpy(element_table[100].short_name, "Fm");
  element_table[100].mass = 257.000000;
/* 101 - Mendelevium */
  strcpy(element_table[101].name, "Mendelevium");
  strcpy(element_table[101].short_name, "Md");
  element_table[101].mass = 258.000000;
/* 102 - Nobelium */
  strcpy(element_table[102].name, "Nobelium");
  strcpy(element_table[102].short_name, "No");
  element_table[102].mass = 259.000000;
/* 103 - Lawrencium */
  strcpy(element_table[103].name, "Lawrencium");
  strcpy(element_table[103].short_name, "Lr");
  element_table[103].mass = 262.000000;
/* 104 - Rutherfordium */
  strcpy(element_table[104].name, "Rutherfordium");
  strcpy(element_table[104].short_name, "Rf");
  element_table[104].mass = 261.000000;
/* 105 - Dubnium */
  strcpy(element_table[105].name, "Dubnium");
  strcpy(element_table[105].short_name, "Db");
  element_table[105].mass = 262.000000;
/* 106 - Seaborgium */
  strcpy(element_table[106].name, "Seaborgium");
  strcpy(element_table[106].short_name, "XX");
  element_table[106].mass = 263.000000;
/* 107 - Bohrium */
  strcpy(element_table[107].name, "Bohrium");
  strcpy(element_table[107].short_name, "Bf");
  element_table[107].mass = 262.000000;
/* 108 - Hassium */
  strcpy(element_table[108].name, "Hassium");
  strcpy(element_table[108].short_name, "Hs");
  element_table[108].mass = 265.000000;
/* 109 - Meitnerium */
  strcpy(element_table[109].name, "Meitnerium");
  strcpy(element_table[109].short_name, "Mt");
  element_table[109].mass = 265.000000;
/* 110 - Mix */
  strcpy(element_table[110].name, "Mix");
  strcpy(element_table[110].short_name, "X");
  element_table[110].mass = 1.0;
}

real ele_mass_from_number(int num)
{
  if (num > N_ELEMENTS || num < 1) {
    return 0.0;
  } else {
    return element_table[num].mass;
  }

}

real ele_mass_from_name(char *name)
{
  int   i = 0;
  if (strlen(name) < 3) {
    for (i = 1; i < N_ELEMENTS; i++)
      if (strcmp(name, element_table[i].short_name) == 0)
	return element_table[i].mass;
  } else {
    for (i = 1; i < N_ELEMENTS; i++)
      if (strcmp(name, element_table[i].name) == 0)
	return element_table[i].mass;
  }
  return 0;
}

int ele_number_from_name(char *name)
{
  int   i = 0;
  if (strlen(name) < 3) {
    for (i = 1; i < N_ELEMENTS; i++)
      if (strcmp(name, element_table[i].short_name) == 0)
	return i;
  } else {
    for (i = 1; i < N_ELEMENTS; i++)
      if (strcmp(name, element_table[i].name) == 0)
	return i;
  }
  return 0;
}
