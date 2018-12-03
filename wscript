#!/usr/bin/env python
# encoding: utf-8

from os import environ, path
from os import mkdir as _mkdir
from sys import platform as _platform
from shutil import copy as _copy
from waflib.Configure import conf
from waflib import Logs
from waflib.Tools.compiler_c import c_compiler

APPNAME = 'potfit'
VERSION = '0.8'
top = '.'
out = 'build'

# Add all options to this list
# [ name, description, defines ]
# Special requirements for options can be handled in
# the _check_enable_options function below

OPTIONS = [
    ['bindist', 'Write a binned radial distribution file', ['BINDIST']],
    ['contrib', 'Enable support for box of contributing particles', ['CONTRIB']],
    ['dsf', 'R|Use damped shifted force approach \n\t(coulomb-based interactions only)', ['DSF']],
    ['evo', 'Use evolutionary algorithm instead of simulated annealing', ['EVO']],
    ['fweight', 'Use modified weights for the forces', ['FWEIGHT']],
    ['mpi', 'Enable MPI parallelization', ['MPI']],
    ['nopunish', 'Disable punishments', ['NOPUNISH']],
    ['resc', 'Enable rescaling (use with care!)', ['RESCALE']],
    ['stress', 'Include stress in fitting process', ['STRESS']],
    ['uq', 'Generate potential ensemble for uncertainty quantification', ['UQ']]
]

# Add all potential models to this list
# [ name, description, supported models, defines, files to compile ]

POTENTIALS = [
    ['pair', 'pair potentials', ['apot', 'tab'], ['PAIR'], ['force_pair.c']],
    ['adp', 'angular dependent potentials', ['apot', 'tab'], ['ADP'], ['force_adp.c']],
    ['ang', 'angular pair potentials', ['apot', 'tab'], ['ANG'], ['force_ang.c']],
    ['ang_elstat', 'angular pair potentials with eletrostatics', [
        'apot', 'tab'], ['ANG', 'COULOMB'], ['force_ang_elstat.c']],
    ['coulomb', 'coulomb interactions', ['apot'], ['COULOMB'], ['force_elstat.c']],
    ['dipole', 'dipole interactions', ['apot'], ['COULOMB', 'DIPOLE'], ['force_elstat.c']],
    ['eam', 'embedded atom method', ['apot', 'tab'], ['EAM'], ['force_eam.c']],
    ['eam_coulomb', 'embedded atom method with coulomb interactions',
        ['apot', 'tab'], ['EAM', 'COULOMB'], ['force_eam_elstat.c']],
    ['eam_dipole', 'embedded atom method with dipole interactions', [
        'apot', 'tab'], ['EAM', 'COULOMB', 'DIPOLE'], ['force_eam_elstat.c']],
    ['kim', 'use OpenKIM framework for force calculation', [], ['KIM'], ['force_kim.c']],
    ['meam', 'modified embedded atom method', ['apot', 'tab'], ['MEAM'], ['force_meam.c']],
    ['stiweb', 'Stillinger-Weber potentials', ['apot'], ['STIWEB'], ['force_stiweb.c']],
    ['tbeam', 'two-band embedded atom method', ['apot', 'tab'], ['EAM', 'TBEAM'], ['force_eam.c']],
    ['tersoff', 'Tersoff potentials', ['apot'], ['TERSOFF'], ['force_tersoff.c']],
    ['tersoffmod', 'modified Tersoff potentials', ['apot'], ['TERSOFF', 'TERSOFFMOD'], ['force_tersoff.c']],
]

# Math libray settings

MKLDIR = '/opt/intel/mkl'

supported_math_libs = [
    ['mkl', 'Intel Math Kernel Library'],
]

if _platform == 'darwin':  # Mac OS
    supported_math_libs.append(['acc', 'Apple Accelerate Framework'])


def options(opt):
    """
    Function which creates the option parser / --help page

    Simple options and new interactions don't need to be added manually.
    Just use the tables above.
    """
    opt.load('compiler_c')
    opts = opt.add_argument_group('potfit general options',
                                  'Please check the explanations on the potfit homepage for more details.')
    # --enable-XXX options are generated automatically from wscript_potfit.py
    for option in OptList():
        opts.add_argument('--enable-{}'.format(option.name), action='store_true',
                          default=False, help=option.description)
    # interactions
    pots = opt.add_argument_group(
        'potfit potential options', 'available interactions in alphabetical order are:\n{}'.format(pl.get_pot_desc()))
    pots.add_argument('-i', '--interaction', action='store', type=str, choices=pl.get_pot_list(),
                      help='one of the interactions listed above', metavar='INTERACTION')
    pots.add_argument('-m', '--model', action='store', type=str, choices=[
                      'apot', 'tab'], help='R|analytic or tabulated potentials\nonly required if the interaction can support both')
    # math libraries
    libs = opt.add_argument_group('potfit math library options', 'available math libraries are:\n{}'.format(
        '\n'.join(sorted(['\t{:<16}{}'.format(x, v) for x, v in supported_math_libs]))))
    libs.add_argument('-l', '--math-lib', action='store', type=str, default='mkl', choices=[name for name, _ in supported_math_libs],
                      help='Select math library to use (default: %(default)s)', metavar='MATHLIB')
    libs.add_argument('--math-lib-base-dir', action='store', type=str,
                      help='Base directory of selected math library')
    # debug options
    debug = opt.add_argument_group('potfit debugging options')
    debug.add_argument('--debug', action='store_true',
                       default=False, help='Build binary with debug information')
    debug.add_argument('--asan', action='store_true', default=False,
                       help='Build binary with address sanitizer support')
    debug.add_argument('--profile', action='store_true',
                       default=False, help='Build binary with profiling support')


def configure(cnf):
    _check_interaction_options(cnf)
    _check_enable_options(cnf)
    _check_compiler_options(cnf)
    _check_math_lib_options(cnf)

    # pass some settings on to the build stage via cnf.env
    cnf.env.model = cnf.options.model
    cnf.env.interaction = cnf.options.interaction
    cnf.env.force_files = pl.get_pot(cnf.options.interaction).files
    cnf.env.append_value('DEFINES_POTFIT', pl.get_pot(cnf.options.interaction).defines)

    # some interactions need more files based on some conditions
    if cnf.options.interaction == 'pair' and cnf.options.model == 'apot':
        cnf.env.force_files.extend(['chempot.c'])

    cnf.env.target_name = _generate_target_name(cnf)

    print('\npotfit has been configured with the following options:')
    print('\n'.join(['{:20} = {}'.format(x, v) for x, v in [['potential model', cnf.options.model], [
          'interaction', cnf.options.interaction], ['math library', cnf.options.math_lib]]]))
    opts = [x[7:] for x in dir(cnf.options) if x.startswith('enable_') and getattr(cnf.options, x)]
    if len(opts):
        for i in opts:
            for j in OPTIONS:
                if i == j[0]:
                    cnf.env.append_value('DEFINES_POTFIT', j[2])
        print('{:20} = {}'.format('options', ', '.join(sorted(opts))))
    print("\nNow run './waf' to start building potfit\n")
    #cnf.define('POTFIT_VERSION', VERSION)
    #cnf.write_config_header('header.h')


def build(bld):
    bld.add_post_fun(_post)
    bld.recurse('src')


def _post(bld):
    """
    Post-build function to copy the binary to the bin/ directory
    """
    if not path.exists('bin'):
        try:
            _mkdir('bin')
        except:
            Logs.error('Could not create bin/ directory')

    try:
        _copy('build/src/' + bld.env.target_name, 'bin/')
        Logs.warn('\n---> Successfully moved {} to bin/ folder <---\n'.format(bld.env.target_name))
    except:
        Logs.error('Could not move potfit binary into bin/ folder')

@conf
def _check_interaction_options(cnf):
    if cnf.options.interaction is None:
        cnf.fatal('No interaction specified, please provide -i/--interaction')

    # common checks for interaction compatibility - don't change
    pot = pl.get_pot(cnf.options.interaction)
    if cnf.options.interaction == 'kim':
        if cnf.options.model != None:
            cnf.fatal('OpenKIM does not support setting a model!')
        cnf.options.model = 'kim'
    elif cnf.options.model is None:
        # if there is only one supported model we don't need the model explicitly
        if pot is not None and len(pot.supp_models) == 1:
            cnf.options.model = pot.supp_models[0]
        else:
            cnf.fatal('No potential model specified for interaction "{}", please provide -m/--model'.format(cnf.options.interaction))

    if not pot.supports_model(cnf.options.model):
        cnf.fatal('Interaction {} does not support model {}!'.format(pot.name, cnf.options.model))


@conf
def _check_enable_options(cnf):
    # array for storing additional source files
    cnf.env.option_files = []

    # check for incompatible options
    if cnf.options.enable_mpi and cnf.options.interaction == 'kim':
        cnf.fatal('KIM does currently not support MPI parallelization')
    if cnf.options.enable_mpi and cnf.options.enable_bindist:
        cnf.fatal('bindist option is not supported for MPI-enabled builds')
    if cnf.options.enable_dsf and cnf.options.interaction not in ['ang_elstat', 'coulomb', 'eam_coulomb']:
        cnf.fatal('DSF can only be used with COULOMB-based interactions and not with {}'.format(cnf.options.interaction))

    # rescale is only allowed with tab and needs additional source files
    if cnf.options.enable_resc:
        if cnf.options.model == 'apot':
            cnf.fatal('Analytic potentials are incompatible with the rescale option!')
        else:
            if cnf.options.interaction == 'meam':
                cnf.env.option_files.append('rescale_meam.c')
            else:
                cnf.env.option_files.append('rescale.c')

    # binned radial distribution is only allowed for tabulated potentials
    if cnf.options.enable_bindist and cnf.options.model == 'apot':
        cnf.fatal('Binned radial distribution files are not available for analytic potentials.')

    # array for storing additional source files
    # use '-...' to remove files added by default in src/wscript
    cnf.env.optimization_files = []

    # differential evolution has its own source file and deactivates simann.c
    if cnf.options.enable_evo:
        cnf.env.optimization_files.extend(['diff_evo.c', '-simann.c'])

    if cnf.options.enable_uq:
        cnf.env.option_files.extend(['uq.c'])


@conf
def _check_compiler_options(cnf):
    # try to detect a suitable compiler
    if cnf.options.interaction == 'kim':
        cnf.check_cfg(package='libkim-api-v2', args=['--cflags', '--libs'], uselib_store='KIM')
    else:
        c_compiler[_platform] = ['icc', 'clang', 'gcc']
    cnf.load('compiler_c')

    if cnf.options.asan:
        if cnf.env.CC_NAME == 'icc':
            cnf.fatal('Intel compiler does not support sanitizer features')

    # check MPI compiler
    if cnf.options.enable_mpi:
        _check_mpi_compiler(cnf)

    # potfit compiler flags
    if cnf.env.CC_NAME == 'icc':
        if cnf.options.debug:
            cnf.env.append_value('CFLAGS_POTFIT', ['-g', '-Wall', '-std=c99'])
        else:
            cnf.env.append_value('CFLAGS_POTFIT', ['-fast', '-xHost', '-std=c99', '-fPIC'])
            cnf.env.append_value('LINKFLAGS_POTFIT', ['-pie'])

    if cnf.env.CC_NAME in ['clang', 'gcc']:
        if cnf.options.debug:
            cnf.env.append_value('CFLAGS_POTFIT', ['-g', '-Wall', '-Werror', '-pedantic', '-std=c99'])
        else:
            cnf.env.append_value('CFLAGS_POTFIT', ['-O3', '-march=native', '-std=c99'])
            cnf.env.append_value('LINKFLAGS_POTFIT', ['-fPIE'])

    # potfit linker flags
    cnf.env.append_value('LINKFLAGS_POTFIT', ['-Wl,--no-undefined,--as-needed,-z,relro,-z,now'])


@conf
def _check_mpi_compiler(cnf):
    cnf.check_cfg(path='mpicc', args='--showme:compile', package='',
                  uselib_store='POTFIT', msg='Checking MPI compiler command')
    cnf.check_cfg(path='mpicc', args='--showme:link', package='',
                  uselib_store='POTFIT', msg='Checking MPI linker flags')
    cnf.check(header_name='mpi.h', features='c cprogram')
    cnf.check_cc(
        fragment='#include <mpi.h>\nint main() { MPI_Init(NULL, NULL); MPI_Finalize(); }',
        execute=True, msg='Compiling test MPI binary', okmsg='OK', errmsg='Failed', use=['POTFIT'])


@conf
def _check_math_lib_options(cnf):
    """
    Function for checking the selected math library

    The math library is checked by compiling some small code snippets
    """
    global MKLDIR

    if cnf.options.math_lib is None:
        if cnf.env.DEST_OS == 'darwin':
            cnf.options.math_lib = 'acc'
        else:
            cnf.options.math_lib = 'mkl'

    # MacOS Accelerate
    if cnf.options.math_lib == 'acc':
        cnf.check(header_name='Accelerate/Accelerate.h', features='c cprogram')
        cnf.env.append_value('DEFINES_POTFIT', ['__ACCELERATE__'])
        cnf.env.append_value('FRAMEWORK_POTFIT', ['Accelerate'])
        return

    # Intel MKL
    if cnf.options.math_lib == 'mkl':
        if cnf.options.math_lib_base_dir:
            MKLDIR = cnf.options.math_lib_base_dir
        cnf.env.append_value('INCLUDES_POTFIT', [MKLDIR + '/include'])
        cnf.env.append_value('LIBPATH_POTFIT', [MKLDIR + '/lib/intel64'])
        cnf.env.append_value('LIB_POTFIT', ['mkl_intel_lp64',
                                            'mkl_sequential', 'mkl_core', 'pthread', 'm'])
        cnf.check(header_name='mkl_vml.h', features='c cprogram', use=['POTFIT'])
        cnf.check(header_name='mkl_lapack.h', features='c cprogram', use=['POTFIT'])
        cnf.check_cc(fragment='''#include <mkl_vml.h>
          #include <mkl_lapack.h>
          int main() {
            double x = 1.0, y = 1.0, result = 0.0;
            vdPow(1, &x, &y, &result);
          }\n''', execute=True, msg='Compiling test MKL binary', okmsg='OK', errmsg='Failed', use=['POTFIT'])
        cnf.env.append_value('DEFINES_POTFIT', ['MKL'])


def _generate_target_name(cnf):
    """
    Function for generating the target name, e.g. potfit_apot_pair_mkl
    """
    name = 'potfit'
    if cnf.options.model and cnf.options.model != cnf.options.interaction:
        name += '_' + cnf.options.model
    name += '_' + cnf.options.interaction
    name += '_' + cnf.options.math_lib
    for item in dir(cnf.options):
        if item.startswith('enable_') and getattr(cnf.options, item):
            name += '_' + item[7:]

    # enable sanitizers
    if cnf.options.asan:
        cnf.env.append_value('CFLAGS_POTFIT', ['-fsanitize=address'])
        cnf.env.append_value('LINKFLAGS_POTFIT', ['-fsanitize=address'])
        name += '_asan'

    if cnf.options.debug:
        name += '_debug'

    return name


"""
Following are some small classes for handling the options and interactions.
"""


class Option:
    def __init__(self, name, description, defines):
        self.name = name
        self.description = description
        self.defines = defines

    def supports_potential(self, model):
        if model in self.supported_models:
            return True
        return False

    def __str__(self):
        desc = 'Option class:'
        for item in dir(self):
            if item[0] != '_' and item[1] != '_':
                desc += '\n\t{} = {}'.format(item, getattr(self, item))
        return desc


class OptList():
    def __init__(self):
        self.list = []
        for i in OPTIONS:
            self.list.append(Option(*i))

    def __iter__(self):
        for i in sorted(self.list, key=lambda x: x.name):
            yield i

    def get_opt(self, name):
        for pot in self.list:
            if (pot.name == name):
                return pot
        return None

    def get_opt_list(self):
        return [i.name for i in self.list]

    def get_opt_desc(self):
        return '\n'.join(sorted(['\t{:<16}{}'.format(i.name, i.description) for i in self.list]))


class Potential:
    def __init__(self, name, description, supp_models, defines, files):
        self.name = name
        self.description = description
        self.supp_models = supp_models
        self.defines = defines
        self.files = files

    def supports_model(self, model):
        if model in self.supp_models:
            return True
        if len(self.supp_models) == 0:
            # only kim supports empty model
            if self.name == 'kim':
                return True
        return False

    def __str__(self):
        desc = 'Potential class:'
        for item in dir(self):
            if item[0] != '_' and item[1] != '_':
                desc += '\n\t{} = {}'.format(item, getattr(self, item))
        return desc


class PotList():
    def __init__(self):
        self.list = []
        for i in POTENTIALS:
            self.list.append(Potential(*i))

    def __iter__(self):
        for i in sorted(self.list):
            yield i

    def get_pot(self, name):
        for pot in self.list:
            if (pot.name == name):
                return pot
        return None

    def get_pot_list(self):
        return [i.name for i in self.list]

    def get_pot_desc(self):
        return '\n'.join(sorted(['\t{:<16}{}'.format(i.name, i.description) for i in self.list]))


pl = PotList()
