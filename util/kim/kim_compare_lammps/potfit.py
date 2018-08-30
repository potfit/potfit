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

import logging
import math
import os
import random

from subprocess import run

logger = logging.getLogger('kim_compare_lammps')

class potfit_run(object):
  def __init__(self, binary, model, config, directory):
    self.binary = binary
    self.config = config
    self.directory = directory
    self.model = model
    self.energy = None
    self.forces = []
    self.__write_parameter_file()
    self.__write_config_file()
    self.__write_potential_file()

  def __write_parameter_file(self):
    filename = os.path.join(self.directory, 'param')
    with open(filename, 'w') as f:
      f.write('ntypes {}\n'.format(self.config.num_atom_types))
      f.write('config config\n')
      f.write('startpot pot\n')
      f.write('endpot pot.end\n')
      f.write('tempfile pot.temp\n')
      f.write('output_prefix out\n')
      f.write('eng_weight 1\n')
      f.write('kim_model_name {}\n'.format(self.model['NAME']))
      f.write('kim_model_params use_default\n')

  def __write_config_file(self):
    filename = os.path.join(self.directory, 'config')
    with open(filename, 'w') as f:
      f.write('#N {} 1\n'.format(self.config.num_atoms()))
      f.write('#C {}\n'.format(' '.join(self.model['SPECIES']) if self.config.num_atom_types > 1 else self.model['SPECIES']))
      f.write('#X {:.08f} {:.08f} {:.08f}\n'.format(*[x * self.config.scale[0] for x in self.config.box[0]]))
      f.write('#Y {:.08f} {:.08f} {:.08f}\n'.format(*[x * self.config.scale[1] for x in self.config.box[1]]))
      f.write('#Z {:.08f} {:.08f} {:.08f}\n'.format(*[x * self.config.scale[2] for x in self.config.box[2]]))
      f.write('#E 0.0\n')
      f.write('#F\n')
      for i in range(len(self.config.atoms)):
        f.write('{:5}\t{:2.8f}\t{:2.8f}\t{:2.8f}'.format(self.config.atom_types[i] - 1, self.config.atoms[i][0], self.config.atoms[i][1], self.config.atoms[i][2]))
        f.write('\t0.0\t0.0\t0.0\n')

  def __write_potential_file(self):
    filename = os.path.join(self.directory, 'pot')
    with open(filename, 'w') as f:
      f.write('#F 5 1\n')
      f.write('#C {}\n'.format(' '.join(self.model['SPECIES']) if self.config.num_atom_types > 1 else self.model['SPECIES']))
      f.write('#E\n\n')

  def run(self):
    res = run([self.binary, 'param'], capture_output=True, cwd=self.directory)
    if res.returncode:
      logger.error('Error running potfit: {}'.format(res.returncode))
      print(self.directory)
      print(res.args)
      print(res.stdout.decode())
      print(res.stderr.decode())
      raise Exception('Error running potfit')
    filename = os.path.join(self.directory, 'out.energy')
    with open(filename, 'r') as f:
      for line in f:
        if line[0] == '#':
          continue
        items = line.split()
        if items[0] == '0':
          self.energy = float(items[3])
          break
    filename = os.path.join(self.directory, 'out.force')
    with open(filename, 'r') as f:
      for line in f:
        if line[0] == '#':
          continue
        items = line.split()
        atom = items[1].split(':')
        atom_idx = int(atom[0])
        if atom_idx >= len(self.forces):
          self.forces.append([999, 999, 999])
        if atom[1] == 'x':
          atom_coord = 0
        elif atom[1] == 'y':
          atom_coord = 1
        else:
          atom_coord = 2
        self.forces[atom_idx][atom_coord] = float(items[4])
    return self.energy, self.forces

  def cleanup(self):
    pass

if __name__ == '__main__':
  print('Please do not run this script directly, use kim_compare_potfit.py instead!')
  sys.exit(-1)
