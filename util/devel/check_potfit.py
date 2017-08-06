#!/usr/bin/env python3

import argparse
import copy
import sys

from itertools import chain, combinations
from subprocess import call

# arrays for supported interactions
tab_interactions = ['pair', 'eam', 'tbeam', 'adp', 'meam']
apot_interactions = ['pair', 'pairang', 'eam', 'tbeam', 'adp', 'meam', 'coulomb', 'dipole',
                     'stiweb', 'tersoff', 'tersoffmod', 'coulomb_eam', 'dipole_eam']
                    # pairang_coulomb NYI
                    #'stiweb', 'tersoff', 'tersoffmod', 'pairang_coulomb', 'coulomb_eam', 'dipole_eam']
kim_interactions = ['dummy']

# arrays for supported options
tab_options = ['evo', 'mpi', 'stress', 'nopunish',
               'dist', 'fweight', 'noresc', 'contrib']
apot_options = ['evo', 'mpi', 'stress', 'nopunish', 'fweight', 'contrib']
kim_options = ['evo', 'stress', 'fweight', 'contrib']

# array of array of unsupported combinations of options
disable = [
    ['dist', 'mpi']
]

# array of options for only certain interactions
special_tab_options = []

special_apot_options = [
    ['coulomb', ['dsf']],
    ['coulomb_eam', ['dsf']],
]

special_kim_options = []


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Compile potfit binary with all possible combinations of options.')

    subs = parser.add_subparsers(
        title="Potential type", description="Choose one of the following types")

    apot = subs.add_parser('apot', help='analytic potentials')
    apot_list = copy.copy(apot_interactions)
    apot_list.sort()
    apot_list.insert(0, 'all')
    apot.add_argument('apot', choices=apot_list)

    tab = subs.add_parser('tab', help='tabulated potentials')
    tab_list = copy.copy(tab_interactions)
    tab_list.sort()
    tab_list.insert(0, 'all')
    tab.add_argument('tab', choices=tab_list)

    kim = subs.add_parser('kim', help='kim interactions')
    kim.add_argument('kim', choices=['all'])

    parser.add_argument('--start', metavar='N', default=1,
                        type=int, required=False, help='start at combination N')
    parser.add_argument('--count', metavar='N', default=-1,
                        type=int, required=False, help='only compile N versions')
    parser.add_argument('--acml', action='store_true')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--nproc', default=1, type=int, required=False,
                        help='number of processer cores for parallel compilation')
    return parser.parse_args()


def get_working_arrays(args):
    string = 'potfit'

    if hasattr(args, 'tab'):
        if (args.tab == 'all'):
            interactions = tab_interactions
        else:
            interactions = [args.tab]
        options = tab_options
        string += '_tab'
        special_options = special_tab_options
    elif hasattr(args, 'apot'):
        if (args.apot == 'all'):
            interactions = apot_interactions
        else:
            interactions = [args.apot]
        options = apot_options
        string += '_apot'
        special_options = special_apot_options
    elif hasattr(args, 'kim'):
        # kim only supports 'all'
        interactions = kim_interactions
        options = kim_options
        string += '_kim'
        special_options = special_kim_options
    else:
        raise RuntimeWarning('No interaction type defined!')

    if args.acml:
        string += '_acml5'

    if args.debug:
        string += '_debug'

    return [interactions, options, string, special_options]


def all_subsets(subset):
    return chain(*map(lambda x: combinations(subset, x), range(0, len(subset) + 1)))


def get_number_of_targets(interactions, options, special_options):
    count = 0
    for interaction in interactions:
        option_list = options[:]
        for special_option in special_options:
            if special_option[0] == interaction:
                option_list += special_option[1]
                break
        count += sum(1 for x in all_subsets(option_list))
    return count


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
    interactions, options, potfit_string, special_options = get_working_arrays(
        args)
    num_targets = get_number_of_targets(interactions, options, special_options)
    print("Building {} possible targets for {}.".format(
        num_targets, potfit_string))

    count = 0
    for interaction in interactions:
        option_list = list(options[:])
        for special_option in special_options:
            if special_option[0] == interaction:
                option_list += special_option[1]
                break
        for build_string in all_subsets(option_list):
            count += 1
            if count < args.start:
                continue
            if args.count != -1 and (count - args.start) >= args.count:
                print("\nMaximum build count ({}) reached.".format(args.count))
                return
            target = potfit_string + '_' + interaction
            if len(build_string):
                target += '_' + '_'.join(build_string)
            if options_are_supported(build_string, disable):
                print("\nBuilding target {} ({}/{})".format(target, count, num_targets))
                call(['make', 'clean'])
                retcode = call(['make', '-j' + str(args.nproc), target])
                if retcode:
                    sys.exit(1)
            else:
                print(
                    "\nTarget {} is disabled ({}/{})".format(target, count, num_targets))
    return


def main():
    args = parse_arguments()
    return build_targets(args)


if __name__ == "__main__":
    main()
