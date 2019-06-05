#!/usr/bin/env python3

import argparse
import copy
import os
import sys

from itertools import chain, combinations
from subprocess import run

# arrays for supported interactions

TAB_INTERACTIONS = [
    'pair',
    'eam',
    'tbeam',
    'adp',
    'meam'
]

APOT_INTERACTIONS = [
    'pair',
    'ang',
    'eam',
    'tbeam',
    'adp',
    'meam',
    'coulomb',
    'dipole',
    'stiweb',
    'tersoff',
    'tersoffmod',
    'ang_elstat',
    'eam_coulomb',
    'eam_dipole'
]

# arrays for supported options

TAB_OPTIONS = [
    'evo',
    'mpi',
    'stress',
    'nopunish',
    'bindist',
    'fweight',
    'resc',
    'contrib'
]

APOT_OPTIONS = [
    'evo',
    'mpi',
    'stress',
    'nopunish',
    'fweight',
    'contrib'
]

KIM_OPTIONS = [
    'evo',
    'stress',
    'fweight',
    'contrib'
]

# array of array of unsupported combinations of options

DISABLED_OPTIONS = [
    ['bindist', 'mpi']
]

# array of options for only certain interactions

SPECIAL_TAB_OPTIONS = []

SPECIAL_APOT_OPTIONS = [
    ['coulomb', ['dsf']],
    ['eam_coulomb', ['dsf']],
]

SPECIAL_KIM_OPTIONS = []


def all_subsets(subset):
    return chain(*map(lambda x: combinations(subset, x), range(0, len(subset) + 1)))


def options_are_supported(subset, disable_array):
    for disabled_set in disable_array:
        disabled_item_count = len(disabled_set)
        for disabled_item in disabled_set:
            if disabled_item in subset:
                disabled_item_count -= 1
        if disabled_item_count == 0:
            return False
    return True


class check_potfit:
    def __init__(self, args):
        self.cmd = []
        self.compiler = args.compiler
        self.max_count = args.count
        self.enable_debug = args.debug
        self.options = []
        self.special_options = []
        self.startpos = args.start

        if hasattr(args, 'enable_apot'):
            self.model = 'apot'
            self.interactions = args.interaction
            self.target_name = 'analytic potentials'
        elif hasattr(args, 'enable_kim'):
            self.model = 'kim'
            self.interactions = None
            self.target_name = 'KIM models'
        elif hasattr(args, 'enable_tab'):
            self.model = 'tab'
            self.interactions = args.interaction
            self.target_name = 'tabulated potentials'

        if args.debug:
            self.cmd.append('--debug')

    def run(self):
        try:
            if self.model == 'kim':
                self._run_kim()
            else:
                self._run()
        except RuntimeWarning as e:
            print('Error: {}'.format(e))
            sys.exit(-1)

    def _run(self):
        self._prepare_arrays()

        num_targets = 0
        for i in self.interactions:
            option_list = self.options[:]
            for s in self.special_options:
                if s[0] == i:
                    option_list.extend(s[1])
                    break
            num_targets += sum(1 for x in all_subsets(option_list))

        print('Building {} possible targets for {}\n'.format(num_targets, self.target_name))

        count = 0
        for i in self.interactions:
            option_list = list(self.options[:])
            for s in self.special_options:
                if s[0] == i:
                    option_list.extend(s[1])
                    break
            for build_string in all_subsets(option_list):
                count += 1
                if count < self.startpos:
                    continue
                if self.max_count != -1 and (count - self.startpos) >= self.max_count:
                    print('\nMaximum build count ({}) reached.'.format(self.max_count))
                    return
                cmds = self.cmd + ['-i', i]
                target_str = 'potfit_' + self.model + '_' + i
                if len(build_string):
                    cmds.extend(['--enable-{}'.format(x) for x in build_string])
                    target_str += '_' + '_'.join(build_string)

                if not options_are_supported(build_string, DISABLED_OPTIONS):
                    print('Target {} is disabled ({}/{})'.format(target_str, count, num_targets))
                    continue

                print('\nBuilding target {} ({}/{})'.format(target_str, count, num_targets))

                run(['./waf', 'clean'])
                my_env = os.environ.copy()
                run_cmd = ['./waf', 'configure', *cmds]
                if 'mpi' in build_string:
                    my_env['OMPI_CC'] = self.compiler
                else:
                    run_cmd.append('--check-c-compiler={}'.format(self.compiler))
                res = run(run_cmd, env=my_env)
                if res.returncode:
                    print('Error when calling ./waf configure ({})'.format(retcode))
                    sys.exit(1)
                res = run(['./waf'])
                if res.returncode:
                    sys.exit(1)

    def _prepare_arrays(self):
        if self.model == 'apot':
            if (self.interactions == 'all'):
                self.interactions = APOT_INTERACTIONS
            else:
                self.interactions = [self.interactions]
            self.options = APOT_OPTIONS
            self.cmd.extend(['-m', 'apot'])
            self.special_options = SPECIAL_APOT_OPTIONS
        elif self.model == 'tab':
            if (self.interactions == 'all'):
                self.interactions = TAB_INTERACTIONS
            else:
                self.interactions = [self.interactions]
            self.options = TAB_OPTIONS
            self.cmd.extend(['-m', 'tab'])
            self.special_options = SPECIAL_TAB_OPTIONS
        else:
            raise RuntimeWarning('No interaction type defined!')

    def _run_kim(self):
        num_targets = 0
        option_list = KIM_OPTIONS
        for s in SPECIAL_KIM_OPTIONS:
            if s[0] == i:
                option_list.extend(s[1])
                break
        num_targets += sum(1 for x in all_subsets(option_list))

        print('Building {} possible targets for {}\n'.format(num_targets, self.target_name))

        count = 0
        option_list = list(KIM_OPTIONS)
        target_str = 'potfit_' + self.model
        for s in SPECIAL_KIM_OPTIONS:
            if s[0] == i:
                option_list.extend(s[1])
                break
        for build_string in all_subsets(option_list):
            count += 1
            if count < self.startpos:
                continue
            if self.max_count != -1 and (count - self.startpos) >= self.max_count:
                print('Maximum build count ({}) reached.'.format(self.max_count))
                return
            cmds = self.cmd + ['-m', 'kim']
            if len(build_string):
                cmds.extend(['--enable-{}'.format(x) for x in build_string])
                target_str += '_' + '_'.join(build_string)

            if not options_are_supported(build_string, DISABLED_OPTIONS):
                print('Target {} is disabled ({}/{})'.format(target_str, count, num_targets))
                continue

            print('Building target {} ({}/{})\n'.format(target_str, count, num_targets))

            run(['./waf', 'clean'])
            my_env = os.environ.copy()
            run_cmd = ['./waf', 'configure', *cmds]
            if 'mpi' in build_string:
                my_env['OMPI_CC'] = self.compiler
            else:
                run_cmd.append('--check-c-compiler={}'.format(self.compiler))
            res = run(run_cmd, env=my_env)
            if res.returncode:
                print('Error when calling ./waf configure ({})'.format(retcode))
                sys.exit(1)
            res = run(['./waf'])
            if res.returncode:
                sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Compile potfit binary with all possible combinations of options.')

    subs = parser.add_subparsers(
        title='Potential type', description='Choose one of the following types')

    apot = subs.add_parser('apot', help='analytic potentials')
    apot_list = copy.copy(APOT_INTERACTIONS)
    apot_list.sort()
    apot_list.insert(0, 'all')
    apot.add_argument('--enable_apot', action='store_true', default=True, help=argparse.SUPPRESS)
    apot.add_argument('interaction', choices=apot_list)

    tab = subs.add_parser('tab', help='tabulated potentials')
    tab_list = copy.copy(TAB_INTERACTIONS)
    tab_list.sort()
    tab_list.insert(0, 'all')
    tab.add_argument('--enable_tab', action='store_true', default=True, help=argparse.SUPPRESS)
    tab.add_argument('interaction', choices=tab_list)

    kim = subs.add_parser('kim', help='kim interactions')
    kim.add_argument('--enable_kim', action='store_true', default=True, help=argparse.SUPPRESS)

    parser.add_argument('--start', metavar='N', default=1,
                        type=int, required=False, help='start at combination N')
    parser.add_argument('--count', metavar='N', default=-1,
                        type=int, required=False, help='only compile N versions')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--compiler', choices=['clang', 'gcc', 'icc'], default='clang')
    return parser.parse_args()


def main():
    check_potfit(parse_arguments()).run()


if __name__ == '__main__':
    main()
