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

from atom_config import atom_config

logger = logging.getLogger('kim_compare_lammps')

class compare(object):
  def __init__(self, energy_l, energy_p, forces_l, forces_p):
    self.energy_l = energy_l
    self.energy_p = energy_p
    self.forces_l = forces_l
    self.forces_p = forces_p

  def run(self):
    logger.info('Comparing energy and {} forces ...'.format(len(self.forces_l)))
    if self._cmp(self.energy_l, self.energy_p):
      logger.error('Energy mismatch: LAMMPS = {} potfit = {}'.format(self.energy_l, self.energy_p))
      raise Exception('Energy mismatch!')
    if len(self.forces_l) != len(self.forces_p):
      logger.error('Number of forces does not match: LAMMPS {} vs potfit {}'.format(len(self.forces_l), len(self.forces_p)))
      raise Exception('Force count mismatch!')
    for i in range(len(self.forces_l)):
      for j in range(3):
        if self._cmp(self.forces_l[i][j], self.forces_p[i][j]):
          logger.error('Force mismatch for atom {} index {}: {} <-> {}'.format(i, j, self.forces_l[i][j], self.forces_p[i][j]))
          raise Exception('Force mismatch!')
    return False

  def _cmp(self, a, b):
    if math.fabs(a) < 1e-10 and math.fabs(b) < 1e-10:
        return False
    return not math.isclose(a, b, rel_tol=0.001)

if __name__ == '__main__':
  print('Please do not run this script directly, use kim_compare_lammps.py instead!')
  sys.exit(-1)
