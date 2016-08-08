#!/usr/bin/env python3

import argparse
import sys

from itertools import chain, combinations
from subprocess import call

def all_subsets(ss):
  return chain(*map(lambda x: combinations(ss, x), range(0, len(ss)+1)))

parser = argparse.ArgumentParser(description='Compile potfit binary with all possible combinations of options.')
parser.add_argument('--start', metavar='N', default=0, type=int, required=False, help='integer starting position')
parser.add_argument('--acml', action='store_true')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--nproc', default=1, type=int, required=False, help='number of processer cores for parallel compilation')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--apot', action='store_true')
group.add_argument('--tab', action='store_true')
group.add_argument('--kim', action='store_true')
args = parser.parse_args()

tab_interactions = [ 'pair', 'eam', 'tbeam', 'adp', 'meam' ]
apot_interactions = [ 'pair', 'eam', 'tbeam', 'adp', 'meam', 'coulomb', 'dipole', 'stiweb', 'tersoff', 'tersoffmod', 'coulomb_eam', 'dipole_eam' ]
kim_interactions = [ 'kim' ]

tab_options = [ 'evo', 'mpi', 'stress', 'nopunish', 'dist', 'fweight', 'noresc', 'contrib' ]
apot_options = [ 'evo', 'mpi', 'stress', 'nopunish', 'fweight', 'contrib' ]
kim_options = [ 'evo', 'stress', 'fweight', 'contrib' ]

disable = [ [ 'dist', 'mpi' ] ]

counter = 0

executable = "potfit"

if args.acml:
    executable += "_acml5"
if args.debug:
    executable += "_debug"

interactions = []
array = []
model = ''

if args.tab:
    interactions = tab_interactions
    array = tab_options
    model = '_tab_'
elif args.apot:
    interactions = apot_interactions
    array = apot_options
    model = '_apot_'
else:
    interactions = kim_interactions
    array = kim_options
    model = '_'

for interaction in interactions:
    for subset in all_subsets(array):
        enable = 1
        for items in disable:
            nmax = (len(items))
            for item in items:
                if item in subset:
                    nmax-=1
            if (nmax==0):
                enable = 0
        if enable:
            if counter > args.start:
                target = executable+model+interaction+'_'+'_'.join(subset)+'_TEST'
                call (['echo', '-e', 'Compiling '+target+'\n' , 'Counter '+str(counter)+'\n' ])
                call (['make', 'clean'])
                retcode = call([ 'make', '-j'+str(args.nproc), target ])
                if retcode:
                    sys.exit(1)
            counter += 1
