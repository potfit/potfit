#!/usr/bin/env python
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

import argparse
import importlib
import logging
import os
import random
import shutil
import sys
import testrun
import time

from pathlib import Path
from subprocess import run

DISABLED_MODELS = [
  'LennardJones612_UniversalShifted__MO_959249795837_003'
]

VERSION = '0.1'

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger('kim_compare_lammps')

def print_banner(line):
  buffer = '#'
  linewidth = 70
  logger.info(buffer * linewidth)
  padding = int((linewidth - 2 * len(buffer) - len(line)) / 2)
  skip = 1 if (len(line) % 2 == 1) else 0
  logger.info(buffer + ' ' * padding + line + ' ' * (padding + skip) + buffer)
  logger.info(buffer * linewidth + '\n')

def check_for_tools(args):
  if not os.path.isfile(args.potfit):
    raise Exception('Provided potfit binary does not exist or is not a file.')
  if not os.access(args.potfit, os.X_OK):
    raise Exception('Provided potfit binary is not executable.')
  if not args.lammps:
    lammps = shutil.which('lmp')
    if not lammps:
      raise Exception('No LAMMPS binary found in PATH, please provide one using the --lammps command line argument')
  else:
    if not os.path.isfile(args.lammps):
      raise Exception('Provided LAMMPS binary does not exist or is not a file.')
    if not os.access(args.lammps, os.X_OK):
      raise Exception('Provided LAMMPS binary is not executable.')
    lammps = args.lammps
  props = shutil.which('kim_read_model_props')
  if not props:
    props = os.path.abspath(args.props)
    if not os.path.isfile(props):
      raise Exception('Provided kim_read_model_props binary does not exist or is not a file.')
    if not os.access(props, os.X_OK):
      raise Exception('Provided kim_read_model_props binary is not executable.')
  have_asap = (importlib.util.find_spec("numpy") and importlib.util.find_spec("asap3.Internal.OpenKIMcalculator") and importlib.util.find_spec("ase"))
  return os.path.abspath(args.potfit), lammps, have_asap, props

def decode_properties(model_props, output = None):
  if not output:
    output = dict()
  for prop in model_props:
    item = prop.split('=', maxsplit=1)
    token = item[1].split('#')
    if len(token) == 1:
      token = token[0]
      if token[0] == '"' and token[-1] == '"':
        token = token.strip('"')
      else:
        token = float(token) if token.find('.') != -1 else int(token)
      output[item[0]] = token
    else:
      output[item[0]] = []
      for sub in token:
        if sub[0] == '"' and sub[-1] == '"':
          sub = sub.strip('"')
        else:
          sub = float(sub) if sub.find('.') != -1 else int(sub)
        output[item[0]].append(sub)
  return output

def get_model_props(model_name, props_binary):
  res = run([props_binary, model_name], capture_output=True)
  if res.returncode != 0:
    return None
  return [model_name, decode_properties(res.stdout.decode().strip().split('\n'))]

def get_model_list(props):
  pkg_config_cmd = os.environ.get('PKG_CONFIG', None) or 'pkg-config'
  libexecdir = run([pkg_config_cmd, 'libkim-api', '--variable=libexecdir'], capture_output=True).stdout.decode().strip()
  models = run([os.path.join(libexecdir, 'kim-api', 'kim-api-collections-info'), 'models'], capture_output=True).stdout.decode().strip().split('\n')
  ret = []
  for line in models:
    p = get_model_props(line.split()[1], props)
    if p and p[1]['PARAMS'] > 0 and not p[1]['NAME'] in DISABLED_MODELS:
      ret.append(p)
  return ret

def main(args):
  print_banner('kim_compare_lammps.py version {}'.format(VERSION))
  potfit, lammps, have_asap, props = check_for_tools(args)
  logger.info('Using potfit binary {}'.format(potfit))
  logger.info('Using LAMMPS binary {}'.format(lammps))
  if have_asap:
    logger.info('Found usable asap installation with OpenKIM support')
  logger.info('Using props binary {}'.format(props))
  models = get_model_list(props)
  logger.info('Performing a total of {} runs using {} KIM models\n'.format(args.runs, len(models)))
  if args.seed:
    random.seed(args.seed)
  basedir = Path('/tmp/kim_compare_lammps_' + time.strftime("%Y%d%m_%H%M%S", time.localtime()))
  Path.mkdir(basedir, exist_ok=False)
  for i in range(args.runs):
    print_banner('Run {}/{}'.format(i + 1, args.runs))
    #logger.info(models)
    model = models[random.randrange(0, len(models))]
    logger.info('Randomly selected model is "{}"'.format(model[0]))
    if testrun.testrun(potfit, lammps, have_asap, model[1], i + 1, basedir).run() == False:
      logger.error('FAIL')
      sys.exit(-1)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Continuously run potfit and LAMMPS using KIM models and compare results.')
  parser.add_argument('--lammps', type=str, help='LAMMPS binary to use')
  parser.add_argument('--potfit', type=str, help='potfit binary to use', required=True)
  parser.add_argument('--props', type=str, help='kim_read_model_props binary', default='../kim_read_model_props/kim_read_model_props')
  parser.add_argument('-r', '--runs', type=int, help='Number of comparison runs (default=10)', default=10)
  parser.add_argument('-s', '--seed', type=int, help='RNG seed')
  parser.add_argument('-l', '--loglevel', type=int, help='Logging level')
  #try:
  logger.setLevel(logging.DEBUG)
  main(parser.parse_args())
  #except Exception as e:
    #print('Error: {}'.format(e))
    #sys.exit(-1)
