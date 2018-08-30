#!/usr/bin/env python3

import argparse
import copy
import sys

from itertools import chain, combinations
from subprocess import call

# arrays for supported interactions
tab_interactions = ['pair', 'eam', 'tbeam', 'adp', 'meam']
apot_interactions = ['pair', 'ang', 'eam', 'tbeam', 'adp', 'meam', 'coulomb', 'dipole',
                     'stiweb', 'tersoff', 'tersoffmod', 'ang_elstat', 'eam_coulomb', 'eam_dipole']
kim_interactions = ['kim']

# arrays for supported options
tab_options = ['evo', 'mpi', 'stress', 'nopunish',
               'dist', 'fweight', 'resc', 'contrib']
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
    ['eam_coulomb', ['dsf']],
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
    parser.add_argument('--compiler', choices=['clang', 'gcc', 'icc'], default='clang')
    return parser.parse_args()


def get_working_arrays(args):
    cmd = []
    target_name = 'potfit'

    if hasattr(args, 'tab'):
        if (args.tab == 'all'):
            interactions = tab_interactions
        else:
            interactions = [args.tab]
        options = tab_options
        cmd.extend(['-m', 'tab'])
        target_name += '_tab'
        special_options = special_tab_options
    elif hasattr(args, 'apot'):
        if (args.apot == 'all'):
            interactions = apot_interactions
        else:
            interactions = [args.apot]
        options = apot_options
        cmd.extend(['-m', 'apot'])
        target_name += '_apot'
        special_options = special_apot_options
    elif hasattr(args, 'kim'):
        # kim only supports 'all'
        interactions = kim_interactions
        options = kim_options
        target_name += '_apot'
        special_options = special_kim_options
    else:
        raise RuntimeWarning('No interaction type defined!')

    if args.acml:
        cmd.extend(['-l', 'acml'])
        target_name += '_acml'

    if args.debug:
        cmd.append('--debug')
        target_name += '_debug'

    return [interactions, options, cmd, target_name, special_options]


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
    try:
        interactions, options, cmd, target, special_options = get_working_arrays(args)
    except RuntimeWarning as e :
        print("Error: {}".format(e))
        sys.exit(-1)

    num_targets = get_number_of_targets(interactions, options, special_options)
    print("Building {} possible targets for {}.".format(num_targets, target))

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
            cmds = cmd + ['-i', interaction]
            target_str = target + '_' + interaction
            if len(build_string):
                cmds.extend(['--enable-{}'.format(x) for x in build_string])
                target_str += '_' + '_'.join(build_string)
            if options_are_supported(build_string, disable):
                print("\nBuilding target {} ({}/{})".format(target_str, count, num_targets))
                call(['./waf', 'clean'])
                retcode = call(['./waf', 'configure', *cmds, '--check-c-compiler={}'.format(args.compiler)])
                if retcode:
                    print("Error when calling ./waf configure ({})".format(retcode))
                    sys.exit(1)
                retcode = call(['./waf'])
                if retcode:
                    sys.exit(1)
            else:
                print(
                    "\nTarget {} is disabled ({}/{})".format(target_str, count, num_targets))
    return


def main():
    build_targets(parse_arguments())


if __name__ == "__main__":
    main()
