import glob
import os
import pytest
import random
import string
import subprocess

class Potfit:
    def __init__(self, location, model, interaction, options = []):
        self.cwd = os.path.dirname(location)
        self.model = model
        self.interaction = interaction
        self.options = options
        for item in [x for x in glob.iglob(os.path.join(self.cwd,'asan_*.*')) if os.path.isfile(x)]:
            os.remove(item)
        self.filenames = []
        self._build()

    def run(self, param_file):
        self.run_with_args([param_file])

    def run_with_args(self, args = []):
        asan_filename = 'asan_{}'.format(''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6)))
        os.environ['ASAN_OPTIONS'] = 'log_path={},exitcode=99,strip_path_prefix={}'.format(asan_filename,os.path.abspath('..') + '/build/../')
        cmd = [os.path.join(os.path.abspath('../bin'), self.binary_name)]
        if len(args):
            cmd.extend(args)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.cwd)
        p.wait()
        self.stdout = p.stdout.read().decode('ascii')
        self.stderr = p.stderr.read().decode('ascii')
        self.returncode = p.returncode
        self._check_asan_fail(asan_filename)

    def clear(self):
        for item in list(set(self.filenames)):
            os.remove(item)

    def has_error_msg(self):
        return '[ERROR]' in self.stderr

    def create_file(self, filename, **kwargs):
        name = os.path.join(self.cwd, filename)
        self.filenames.append(name)
        f = open(name, 'w')
        if 'permission' in kwargs:
            os.chmod(name, kwargs['permission'])
        return f

    def _build(self):
        p = subprocess.Popen(['./waf', 'distclean'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd='..')
        p.wait()
        if p.returncode:
            pytest.fail('error calling "waf distclean"')
        p = subprocess.Popen(self._get_conf_cmd(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd='..')
        p.wait()
        if p.returncode:
            print(p.stderr.read().decode('ascii'))
            pytest.fail('error calling "waf configure"')
        p = subprocess.Popen(['./waf', 'build'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd='..')
        p.wait()
        if p.returncode:
            print(p.stderr.read().decode('ascii'))
            pytest.fail('error calling "waf build"')
        out = p.stdout.read().decode('ascii').split('\n')
        self.binary_name = [x for x in out if 'Linking' in x][0].split(' ')[2].split('/')[-1]
        if not len(self.binary_name):
            pytest.fail('error reading binary name')

    def _get_conf_cmd(self):
        cmd = ['./waf', 'configure', '-c', 'no', '-m', self.model, '-i', self.interaction]
        for opt in self.options:
            cmd.append('--enable-{}'.format(opt))
        cmd.append('--asan')
        cmd.append('--debug')
        cmd.append('--check-c-compiler=clang')
        return cmd

    def _check_asan_fail(self, filename):
        if self.returncode != 99:
            return
        filenames = [x[len(os.path.abspath('..'))+7:] for x in glob.iglob(os.path.join(self.cwd,'{}.*'.format(filename))) if os.path.isfile(x)]
        pytest.fail('address sanitizer detected an error, please check {}'.format('\n'.join(filenames)))
