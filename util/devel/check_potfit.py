#!/usr/bin/env python3

import argparse
import sys

from itertools import chain, combinations
from subprocess import call

tab_interactions = [ 'pair', 'eam', 'tbeam', 'adp', 'meam' ]
apot_interactions = [ 'pair', 'eam', 'tbeam', 'adp', 'meam', 'coulomb', 'dipole', 'stiweb', 'tersoff', 'tersoffmod', 'coulomb_eam', 'dipole_eam' ]
kim_interactions = [ 'kim' ]

tab_options = [ 'evo', 'mpi', 'stress', 'nopunish', 'dist', 'fweight', 'noresc', 'contrib' ]
apot_options = [ 'evo', 'mpi', 'stress', 'nopunish', 'fweight', 'contrib' ]
kim_options = [ 'evo', 'stress', 'fweight', 'contrib' ]

disable = [ \
    [ 'dist', 'mpi' ] \
]

def parse_arguments():
    parser = argparse.ArgumentParser(description='Compile potfit binary with all possible combinations of options.')
    parser.add_argument('--start', metavar='N', default=1, type=int, required=False, help='integer starting position')
    parser.add_argument('--acml', action='store_true')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--nproc', default=1, type=int, required=False, help='number of processer cores for parallel compilation')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--apot', action='store_true')
    group.add_argument('--tab', action='store_true')
    group.add_argument('--kim', action='store_true')
    return parser.parse_args()

def get_working_arrays(args):
    s = 'potfit'
    if args.tab:
        i = tab_interactions
        o = tab_options
        s += '_tab'
    elif args.apot:
        i = apot_interactions
        o = apot_options
        s += '_apot'
    else:
        i = kim_interactions
        o = kim_options
        s += ''

    if args.acml:
        s += '_acml5'

    if args.debug:
        s += '_debug'

    return [i, o, s]

def all_subsets(subset):
    return chain(*map(lambda x: combinations(subset, x), range(0, len(subset)+1)))

def get_number_of_targets(i, o):
    count = sum(1 for x in all_subsets(o))
    return count * len(i)

def options_are_supported(subset, disable_array):
    for disabled_set in disable_array:
        disabled_item_count = len(disabled_set)
        for disabled_item in disabled_set:
            if disabled_item in subset:
                disabled_item_count -= 1
        if disabled_item_count == 0:
            return False
    return True

def build_targets(args):
    i, o, p = get_working_arrays(args)
    num_targets = get_number_of_targets(i, o)
    print("Building {} possible targets for {}.".format(num_targets, p))

    count = 0
    for interaction in i:
        for options in all_subsets(o):
            count += 1
            target = p + '_' + interaction
            if len(options):
                target += '_' + '_'.join(options)
            if options_are_supported(options, disable):
                print("\nBuilding target {} ({}/{})".format(target, count, num_targets))
                call (['make', 'clean'])
                retcode = call([ 'make', '-j'+str(args.nproc), target ])
                if retcode:
                    sys.exit(1)
            else:
                print("\nTarget {} is disabled ({}/{})".format(target, count, num_targets))
    return

def main():
    args = parse_arguments()
    return build_targets(args)

if __name__ == "__main__":
    main()
