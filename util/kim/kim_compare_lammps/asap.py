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

import os
import sys

from subprocess import run

from numpy import *
from asap3 import Atoms, OpenKIMcalculator
from ase import Atom

class asap_run(object):
  def __init__(self, model, config, directory):
    self.config = config
    self.directory = directory
    self.model = model
    self.energy = None
    self.forces = []
    self.__write_config()
    self.__write_input()

  def __write_config(self):
    self.atoms = []
    filename = os.path.join(self.directory, 'config')
    with open(filename, 'w') as f:
      f.write('[\n')
      for i, atom in enumerate(self.config.atoms):
        type_str = self.model['SPECIES'] if self.config.num_atom_types == 1 else self.model['SPECIES'][self.config.atom_types[i]]
        f.write('  Atom(\'{}\', ({}, {}, {})),\n'.format(type_str, atom[0], atom[1], atom[2]))
        self.atoms.append(Atom(type_str, (atom[0], atom[1], atom[2])))
      f.write(']')

  def __write_input(self):
    filename = self.directory / 'input.py'
    with open(filename, 'w') as f:
      f.write('from numpy import *\n')
      f.write('from asap3 import Atoms, OpenKIMcalculator\n')
      f.write('from ase import Atom\n')
      f.write('with open(\'config\', \'r\') as f:\n')
      f.write('  atom_list = eval(f.read())\n')

      self.cell_vectors = []
      f.write('cell_vectors = [\n')
      for i in range(3):
        self.cell_vectors.append([x * self.config.scale[i] for x in self.config.box[i]])
        f.write('  [{:.08f}, {:.08f}, {:.08f}],\n'.format(*self.cell_vectors[-1]))
      f.write(']\n\n')

      f.write('atoms = Atoms(symbols=atom_list, cell=cell_vectors, pbc=True)\n')
      f.write('atoms.set_calculator(OpenKIMcalculator(\'{}\'))\n'.format(self.model['NAME']))
      f.write('asap_epot = atoms.get_potential_energy() / len(atoms)\n')
      f.write('asap_forces = atoms.get_forces()\n')

  def run(self):
    atoms = Atoms(symbols=self.atoms, cell=self.cell_vectors, pbc=True)
    atoms.set_calculator(OpenKIMcalculator(self.model['NAME']))
    self.energy = atoms.get_potential_energy()
    self.forces = atoms.get_forces()
    return self.energy, self.forces

  def cleanup(self):
    pass

if __name__ == '__main__':
  print('Please do not run this script directly, use kim_compare_lammps.py instead!')
  sys.exit(-1)
