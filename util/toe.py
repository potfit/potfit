#!/usr/bin/env python
################################################################
#
# toe.py:
#   table of elements, containing names, symbols and masses
#
################################################################
#
#   Copyright 2013
#             Institute for Theoretical and Applied Physics
#             University of Stuttgart, D-70550 Stuttgart, Germany
#             https://www.potfit.net/
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


class Element:
	def __init__(self, number, name, symbol, weight):
		self.number = number
		self.name = name
		self.symbol = symbol
		self.weight = weight

el_table = [
    Element(1, 'Hydrogen', 'H', 1.0079),
    Element(2, 'Helium', 'He', 4.0026),
    Element(3, 'Lithium', 'Li', 6.941000),
    Element(4, 'Beryllium', 'Be', 9.012180),
    Element(5, 'Boron', 'B', 10.810000),
    Element(6, 'Carbon', 'C', 12.011000),
    Element(7, 'Nitrogen', 'N', 14.006740),
    Element(8, 'Oxygen', 'O', 15.999400),
    Element(9, 'Fluorine', 'F', 18.998403),
    Element(10, 'Neon', 'Ne', 20.179000),
    Element(11, 'Sodium', 'Na', 22.989770),
    Element(12, 'Magnesium', 'Mg', 24.305000),
    Element(13, 'Aluminum', 'Al', 26.981540),
    Element(14, 'Silicon', 'Si', 28.086000),
    Element(15, 'Phosphorus', 'P', 30.973760),
    Element(16, 'Sulfur', 'S', 32.060000),
    Element(17, 'Chlorine', 'Cl', 35.453000),
    Element(18, 'Argon', 'Ar', 39.948000),
    Element(19, 'Potassium', 'K', 39.098000),
    Element(20, 'Calcium', 'Ca', 40.080000),
    Element(21, 'Scandium', 'Sc', 44.955900),
    Element(22, 'Titanium', 'Ti', 47.900000),
    Element(23, 'Vanadium', 'V', 50.941400),
    Element(24, 'Cromium', 'Cr', 51.996000),
    Element(25, 'Manganese', 'Mn', 54.938000),
    Element(26, 'Iron', 'Fe', 55.847000),
    Element(27, 'Cobalt', 'Co', 58.933200),
    Element(28, 'Nickel', 'Ni', 58.700000),
    Element(29, 'Copper', 'Cu', 63.546000),
    Element(30, 'Zinc', 'Zn', 65.380000),
    Element(31, 'Gallium', 'Ga', 69.720000),
    Element(32, 'Germanium', 'Ge', 72.630000),
    Element(33, 'Arsenic', 'As', 74.921600),
    Element(34, 'Selenium', 'Se', 78.960000),
    Element(35, 'Bromine', 'Br', 79.904000),
    Element(36, 'Krypton', 'Kr', 83.800000),
    Element(37, 'Rubidium', 'Rb', 85.467800),
    Element(38, 'Strontium', 'Sr', 87.620000),
    Element(39, 'Yttrium', 'Y', 88.905900),
    Element(40, 'Zirconium', 'Zr', 91.220000),
    Element(41, 'Niobium', 'Nb', 92.906400),
    Element(42, 'Molybdenum', 'Mo', 95.940000),
    Element(43, 'Technetium', 'Tc', 97.000000),
    Element(44, 'Ruthenium', 'Ru', 101.070000),
    Element(45, 'Rhodium', 'Rh', 102.905500),
    Element(46, 'Palladium', 'Pd', 106.400000),
    Element(47, 'Silver', 'Ag', 107.868000),
    Element(48, 'Cadmium', 'Cd', 112.400000),
    Element(49, 'Indium', 'In', 114.820000),
    Element(50, 'Tin', 'Sn', 118.690000),
    Element(51, 'Antimony', 'Sb', 121.750000),
    Element(52, 'Tellurium', 'Te', 127.600000),
    Element(53, 'Iodine', 'I', 126.904500),
    Element(54, 'Xenon', 'Xe', 131.300000),
    Element(55, 'Cesium', 'Cs', 132.905400),
    Element(56, 'Barium', 'Ba', 137.330000),
    Element(57, 'Lanthanum', 'La', 138.905500),
    Element(58, 'Cerium', 'Ce', 140.120000),
    Element(59, 'Praseodynium', 'Pr', 140.907700),
    Element(60, 'Neodymium', 'Nd', 144.240000),
    Element(61, 'Promethium', 'Pm', 145.000000),
    Element(62, 'Samarium', 'Sm', 150.400000),
    Element(63, 'Europium', 'Eu', 151.960000),
    Element(64, 'Gadolinium', 'Gd', 157.250000),
    Element(65, 'Terbium', 'Tb', 158.925400),
    Element(66, 'Dysprosium', 'Dy', 162.500000),
    Element(67, 'Holmium', 'Ho', 164.930400),
    Element(68, 'Erbium', 'Er', 167.260000),
    Element(69, 'Thulium', 'Tm', 168.934200),
    Element(70, 'Ytterbium', 'Yb', 173.040000),
    Element(71, 'Lutetium', 'Lu', 174.970000),
    Element(72, 'Hafnium', 'Hf', 178.490000),
    Element(73, 'Tantalum', 'Ta', 180.947900),
    Element(74, 'Tungsten', 'W', 183.500000),
    Element(75, 'Rhenium', 'Re', 186.207000),
    Element(76, 'Osmium', 'Os', 190.200000),
    Element(77, 'Iridium', 'Ir', 192.220000),
    Element(78, 'Platinum', 'Pt', 195.090000),
    Element(79, 'Gold', 'Au', 196.966500),
    Element(80, 'Mercury', 'Hg', 200.590000),
    Element(81, 'Thallium', 'Tl', 204.370000),
    Element(82, 'Lead', 'Pb', 207.200000),
    Element(83, 'Bismuth', 'Bi', 208.980400),
    Element(84, 'Polonium', 'Po', 209.000000),
    Element(85, 'Astatine', 'At', 210.000000),
    Element(86, 'Radon', 'Rn', 222.000000),
    Element(87, 'Francium', 'Fr', 223.000000),
    Element(88, 'Radium', 'Ra', 226.025400),
    Element(89, 'Actinium', 'Ac', 227.000000),
    Element(90, 'Thorium', 'Th', 232.038100),
    Element(91, 'Protactinium', 'Pa', 231.035900),
    Element(92, 'Uranium', 'U', 238.029000),
    Element(93, 'Neptunium', 'Np', 237.048200),
    Element(94, 'Plutonium', 'Pu', 244.000000),
    Element(95, 'Americium', 'Am', 243.000000),
    Element(96, 'Curium', 'Cm', 247.000000),
    Element(97, 'Berkelium', 'Bk', 247.000000),
    Element(98, 'Californium', 'Cf', 251.000000),
    Element(99, 'Einsteinium', 'Es', 252.000000),
    Element(100, 'Fermium', 'Fm', 257.000000),
    Element(101, 'Mendelevium', 'Md', 258.000000),
    Element(102, 'Nobelium', 'No', 259.000000),
    Element(103, 'Lawrencium', 'Lr', 262.000000),
    Element(104, 'Rutherfordium', 'Rf', 261.000000),
    Element(105, 'Dubnium', 'Db', 262.000000),
    Element(106, 'Seaborgium', 'XX', 263.000000),
    Element(107, 'Bohrium', 'Bf', 262.000000),
    Element(108, 'Hassium', 'Hs', 265.000000),
    Element(109, 'Meitnerium', 'Mt', 265.000000),
    Element(110, 'Mix', 'X', 1.0)
        ]

def getMassByNumber(number):
    if number > 110:
        return 0.0
    else:
        return el_table[number-1].weight

def getMassByName(name):
    if len(name)<3:
        for i in el_table:
            if name.strip() == i.symbol:
                return i.weight
    else:
        for i in el_table:
            if name.strip() == i.name:
                return i.weight
    return 0.0

def getNumberByMass(mass):
    for i in el_table:
        if abs(mass-i.weight)<=1e-2:
            return i.number
    return 0

def getNumberByName(name):
    if len(name)<3:
        for i in el_table:
            if name.strip() == i.symbol:
                return i.number
    else:
        for i in el_table:
            if name.strip() == i.name:
                return i.number
    return 0

def getNameByMass(mass):
    for i in el_table:
        if abs(mass-i.weight)<=1e-2:
            return i.name
    return ''

def getSymbolByMass(mass):
    for i in el_table:
        if abs(mass-i.weight)<=1e-2:
            return i.symbol
    return ''

def getNameByNumber(number):
    if number > 110:
        return ''
    else:
        return el_table[number-1].name

def getSymbolByNumber(number):
    if number > 110:
        return ''
    else:
        return el_table[number-1].symbol

if __name__ == "__main__":
    print "\nPlease do not call this script directly."
