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
import sys
import tempfile

from asap import asap_run
from atom_config import atom_config
from compare import compare
from lammps import lammps_run
from pathlib import Path
from potfit import potfit_run

logger = logging.getLogger('kim_compare_lammps')

class testrun(object):
  def __init__(self, potfit, lammps, have_asap, model, run_index, basedir):
    self.potfit = potfit
    self.model = model
    self.config = atom_config(model)
    self.basedir = basedir / 'run_{:03d}'.format(run_index)
    Path.mkdir(self.basedir, exist_ok=True)
    self.have_asap = have_asap
    if self.have_asap:
      self.asapdir = self.basedir / 'asap'
      Path.mkdir(self.asapdir, exist_ok=True)
      logger.info('Generating asap input data ...')
      self.asap = asap_run(model, self.config, self.asapdir)
    self.lammpsdir = self.basedir / 'lammps'
    Path.mkdir(self.lammpsdir, exist_ok=True)
    logger.info('Generating LAMMPS input data ...')
    self.lammps = lammps_run(lammps, model, self.config, self.lammpsdir)
    self.potfitdir = self.basedir / 'potfit'
    Path.mkdir(self.potfitdir, exist_ok=True)
    logger.info('Generating potfit input data ...')
    self.potfit = potfit_run(potfit, model, self.config, self.potfitdir)

  def run(self):
    res = True

    energies = []
    forces = []

    if self.have_asap:
      try:
        logger.info('Running ASAP calculation ...')
        energy_a, forces_a = self.asap.run()
        energies.append(energy_a)
        forces.append(forces_a)
      except:
        pass

    try:
      logger.info('Running LAMMPS calculation ...')
      energy_l, forces_l = self.lammps.run()
      energies.append(energy_l)
      forces.append(forces_l)
    except:
      pass

    try:
      logger.info('Running potfit calculation ...')
      energy_p, forces_p = self.potfit.run()
      energies.append(energy_p)
      forces.append(forces_p)
    except:
      pass

    try:
      compare(energies, forces).run()
    except Exception as e:
      logger.error(e)
      res = False
    finally:
      if self.have_asap:
        self.asap.cleanup()
      self.lammps.cleanup()
      self.potfit.cleanup()

    return res

if __name__ == '__main__':
  print('Please do not run this script directly, use kim_compare_lammps.py instead!')
  sys.exit(-1)
