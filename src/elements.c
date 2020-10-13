/****************************************************************
 *
 * elements.c: data and routines for periodic table
 *
 ****************************************************************
 *
 * Copyright 2002-2018 - the potfit development team
 *
 * https://www.potfit.net/
 *
 ****************************************************************
 *
 * This file is part of potfit.
 *
 * potfit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * potfit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

#include "potfit.h"

#include "elements.h"

const struct element_table {
  const char* name;
  const char* short_name;
  double mass;
} elements[] = {
    {"Hydrogen", "H", 1.0079},          {"Helium", "He", 4.0026},
    {"Lithium", "Li", 6.941000},        {"Beryllium", "Be", 9.012180},
    {"Boron", "B", 10.810000},          {"Carbon", "C", 12.011000},
    {"Nitrogen", "N", 14.006740},       {"Oxygen", "O", 15.999400},
    {"Fluorine", "F", 18.998403},       {"Neon", "Ne", 20.179000},
    {"Sodium", "Na", 22.989770},        {"Magnesium", "Mg", 24.305000},
    {"Aluminum", "Al", 26.981540},      {"Silicon", "Si", 28.086000},
    {"Phosphorus", "P", 30.973760},     {"Sulfur", "S", 32.060000},
    {"Chlorine", "Cl", 35.453000},      {"Argon", "Ar", 39.948000},
    {"Potassium", "K", 39.098000},      {"Calcium", "Ca", 40.080000},
    {"Scandium", "Sc", 44.955900},      {"Titanium", "Ti", 47.900000},
    {"Vanadium", "V", 50.941400},       {"Cromium", "Cr", 51.996000},
    {"Manganese", "Mn", 54.938000},     {"Iron", "Fe", 55.847000},
    {"Cobalt", "Co", 58.933200},        {"Nickel", "Ni", 58.700000},
    {"Copper", "Cu", 63.546000},        {"Zinc", "Zn", 65.380000},
    {"Gallium", "Ga", 69.720000},       {"Germanium", "Ge", 72.630000},
    {"Arsenic", "As", 74.921600},       {"Selenium", "Se", 78.960000},
    {"Bromine", "Br", 79.904000},       {"Krypton", "Kr", 83.800000},
    {"Rubidium", "Rb", 85.467800},      {"Strontium", "Sr", 87.620000},
    {"Yttrium", "Y", 88.905900},        {"Zirconium", "Zr", 91.220000},
    {"Niobium", "Nb", 92.906400},       {"Molybdenum", "Mo", 95.940000},
    {"Technetium", "Tc", 97.000000},    {"Ruthenium", "Ru", 101.070000},
    {"Rhodium", "Rh", 102.905500},      {"Palladium", "Pd", 106.400000},
    {"Silver", "Ag", 107.868000},       {"Cadmium", "Cd", 112.400000},
    {"Indium", "In", 114.820000},       {"Tin", "Sn", 118.690000},
    {"Antimony", "Sb", 121.750000},     {"Tellurium", "Te", 127.600000},
    {"Iodine", "I", 126.904500},        {"Xenon", "Xe", 131.300000},
    {"Cesium", "Cs", 132.905400},       {"Barium", "Ba", 137.340000},
    {"Lanthanum", "La", 138.905500},    {"Cerium", "Ce", 140.120000},
    {"Praseodynium", "Pr", 140.907700}, {"Neodymium", "Nd", 144.240000},
    {"Promethium", "Pm", 145.000000},   {"Samarium", "Sm", 150.400000},
    {"Europium", "Eu", 151.960000},     {"Gadolinium", "Gd", 157.250000},
    {"Terbium", "Tb", 158.925400},      {"Dysprosium", "Dy", 162.500000},
    {"Holmium", "Ho", 164.930400},      {"Erbium", "Er", 167.260000},
    {"Thulium", "Tm", 168.934200},      {"Ytterbium", "Yb", 173.040000},
    {"Lutetium", "Lu", 174.970000},     {"Hafnium", "Hf", 178.490000},
    {"Tantalum", "Ta", 180.947900},     {"Tungsten", "W", 183.500000},
    {"Rhenium", "Re", 186.207000},      {"Osmium", "Os", 190.200000},
    {"Iridium", "Ir", 192.220000},      {"Platinum", "Pt", 195.090000},
    {"Gold", "Au", 196.966500},         {"Mercury", "Hg", 200.590000},
    {"Thallium", "Tl", 204.370000},     {"Lead", "Pb", 207.200000},
    {"Bismuth", "Bi", 208.980400},      {"Polonium", "Po", 209.000000},
    {"Astatine", "At", 210.000000},     {"Radon", "Rn", 222.000000},
    {"Francium", "Fr", 223.000000},     {"Radium", "Ra", 226.025400},
    {"Actinium", "Ac", 227.000000},     {"Thorium", "Th", 232.038100},
    {"Protactinium", "Pa", 231.035900}, {"Uranium", "U", 238.029000},
    {"Neptunium", "Np", 237.048200},    {"Plutonium", "Pu", 244.000000},
    {"Americium", "Am", 243.000000},    {"Curium", "Cm", 247.000000},
    {"Berkelium", "Bk", 247.000000},    {"Californium", "Cf", 251.000000},
    {"Einsteinium", "Es", 252.000000},  {"Fermium", "Fm", 257.000000},
    {"Mendelevium", "Md", 258.000000},  {"Nobelium", "No", 259.000000},
    {"Lawrencium", "Lr", 262.000000},   {"Rutherfordium", "Rf", 261.000000},
    {"Dubnium", "Db", 262.000000},      {"Seaborgium", "XX", 263.000000},
    {"Bohrium", "Bf", 262.000000},      {"Hassium", "Hs", 269.000000},
    {"Meitnerium", "Mt", 278.000000},   {"Darmstadtium", "Ds", 281.000000},
};

/****************************************************************
  ele_mass_from_number
****************************************************************/

double ele_mass_from_number(int num)
{
  const size_t size = sizeof(elements) / sizeof(elements[0]);
  if (num > size || num < 1)
    return 0.0;
  else
    return elements[num - 1].mass;
}

/****************************************************************
  ele_mass_from_name
****************************************************************/

double ele_mass_from_name(const char* name)
{
  const size_t size = sizeof(elements) / sizeof(elements[0]);
  if (strlen(name) < 3) {
    for (size_t i = 0; i < size; ++i)
      if (strcmp(name, elements[i].short_name) == 0)
        return elements[i].mass;
  } else {
    for (size_t i = 0; i < size; ++i)
      if (strcmp(name, elements[i].name) == 0)
        return elements[i].mass;
  }

  return 0.0;
}

/****************************************************************
  ele_number_from_name
****************************************************************/

int ele_number_from_name(const char* name)
{
  const size_t size = sizeof(elements) / sizeof(elements[0]);
  if (strlen(name) < 3) {
    for (size_t i = 0; i < size; ++i)
      if (strcmp(name, elements[i].short_name) == 0)
        return i + 1;
  } else {
    for (size_t i = 0; i < size; ++i)
      if (strcmp(name, elements[i].name) == 0)
        return i + 1;
  }

  return 0;
}
