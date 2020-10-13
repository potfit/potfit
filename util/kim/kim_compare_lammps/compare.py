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
  def __init__(self, energies, forces):
    self.energies = energies
    self.forces = forces

  def run(self):
    logger.info('Comparing energy and {} forces ...'.format(len(self.forces[0])))
    if self._cmp(self.energies):
      if len(self.energies) == 3:
        logger.error('Energy mismatch: ASAP = {} LAMMPS = {} potfit = {}'.format(self.energies[0], self.energies[1], self.energies[2]))
      else:
        logger.error('Energy mismatch: LAMMPS = {} potfit = {}'.format(self.self.energies[0], self.self.energies[1]))
      raise Exception('Energy mismatch!')
    if len(self.forces[0]) != len(self.forces[1]) or (len(self.forces) == 3 and (len(self.forces[1]) != len(self.forces[2]))):
      logger.error('Number of forces does not match: {}'.format([len(x) for x in self.forces]))
      raise Exception('Force count mismatch!')
    for i in range(len(self.forces[0])):
      for j in range(3):
        if self._cmp([x[i][j] for x in self.forces]):
          logger.error('Force mismatch for atom {} index {}: {}'.format(i, j, [x[i][j] for x in self.forces]))
          raise Exception('Force mismatch!')
    return False

  def _cmp(self, e):
    if len(e) == 3:
      return self._cmp3(e)
    if math.fabs(e[0]) < 1e-10 and math.fabs(e[1]) < 1e-10:
        return False
    if math.fabs(e[0]) < 1e-7 or math.fabs(e[1]) < 1e-7:
        return not math.isclose(e[0], e[1], rel_tol=0.01)
    return not math.isclose(e[0], e[1], rel_tol=0.001)

  def _cmp3(self, e):
    if math.fabs(e[0]) < 1e-10 and math.fabs(e[1]) < 1e-10 and math.fabs(e[2]) < 1e-10:
        return False
    if math.fabs(e[0]) < 1e-7 or math.fabs(e[1]) < 1e-7 or math.fabs(e[2]) < 1e-7:
        return not math.isclose(e[0], e[1], rel_tol=0.01) or not math.isclose(e[1], e[2], rel_tol=0.01)
    return not math.isclose(e[0], e[1], rel_tol=0.001) or not math.isclose(e[1], e[2], rel_tol=0.001)

if __name__ == '__main__':
  print('Please do not run this script directly, use kim_compare_lammps.py instead!')
  sys.exit(-1)
