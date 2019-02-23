################################################################
#
# kim_compare_lammps
#
################################################################
#
# Copyright 2018 the potfit development team
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the “Software”), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall
# be included in all copies or substantial portions of the
# Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY
# KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE
# AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# https://www.potfit.net/
#
#################################################################

import random

from itertools import product

LATTICE_TYPES = ['Primitive', 'BCC', 'FCC']

def copy_atom(a, box, disp):
  res = []
  for i in range(len(a)):
    res.append(round(a[i] + disp[0] * box[0][i] + disp[1] * box[1][i] + disp[2] * box[2][i], 5))
  return res

class atom_config(object):
  def __init__(self, model):
    self.atoms = []
    self.atom_types = []
    self.box = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    self.model = model
    self.num_atom_types = model['NUM_SPECIES']
    self.scale = [1, 1, 1]
    self.__generate(LATTICE_TYPES[random.randint(0, len(LATTICE_TYPES) - 1)])

  def num_atoms(self):
    return len(self.atoms)

  def __generate(self, lattice_type):
    if lattice_type == 'Primitive':
      base_atoms = [[0, 0, 0]]
    elif lattice_type == 'BCC':
      base_atoms = [[0, 0, 0], [0.5, 0.5, 0.5]]
    elif lattice_type == 'FCC':
      base_atoms = [[0, 0, 0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]
    elif lattice_type == 'Diamond':
      base_atoms = [[0, 0, 0], [0.5, 0.5, 0.5]]
    return self.__generate_Primitive(base_atoms)

  def __generate_Primitive(self, base_atoms):
    self.scale = [random.randint(3, 10) for x in range(3)]
    lattice_constant = round(random.uniform(1.5, self.model['INFLUENCE_DISTANCE']), 5)
    for i in range(len(self.box)):
      self.box[i] = [x * lattice_constant for x in self.box[i]]
    for i in range(len(base_atoms)):
      base_atoms[i] = [x * lattice_constant for x in base_atoms[i]]
    for disp in product(range(self.scale[0]), range(self.scale[1]), range(self.scale[2])):
      for atom in base_atoms:
        self.atoms.append(copy_atom(atom, self.box, list(disp)))
        self.atom_types.append(random.randint(1, self.num_atom_types))

if __name__ == '__main__':
  print('Please do not run this script directly, use kim_compare_lammps.py instead!')
  sys.exit(-1)
