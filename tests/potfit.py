import glob
import math
import os
import pytest
import random
import re
import string
import subprocess

from itertools import product

class Potfit:
    def __init__(self, location, model, interaction, options = [], **kwargs):
        self.cwd = os.path.dirname(location)
        self.model = model
        self.interaction = interaction
        self.options = options
        self.filenames = []
        self.mindist_pattern = re.compile('Minimal Distances Matrix:[\w\W]*?\n\n', re.MULTILINE)
        if 'git_patch' in kwargs:
          self._build(kwargs['git_patch'])
        else:
          self._build()

    def reset(self):
        self.filenames = []
        self.stdout = str()
        self.stderr = str()
        self.endpot = str()
        self.force = str()
        self.energy = str()
        self.stress = str()
        self.error = str()
        self.punish = str()
        self.rho_loc = str()
        self.returncode = 0
        self.config = None

    def create_file(self, filename, **kwargs):
        name = os.path.join(self.cwd, filename)
        self.filenames.append(name)
        f = open(name, 'w')
        if 'permission' in kwargs:
            os.chmod(name, kwargs['permission'])
        return f

    def create_param_file(self, **kwargs):
        f = self.create_file('param_file')
        if 'ntypes' in kwargs:
            f.write('ntypes {}\n'.format(int(kwargs['ntypes'])))
            del(kwargs['ntypes'])
        else:
            f.write('ntypes 1\n')
        if 'startpot' in kwargs:
            f.write('startpot {}\n'.format(kwargs['startpot']))
        else:
            f.write('startpot startpot\n')
        if 'endpot' in kwargs:
            f.write('endpot {}\n'.format(kwargs['endpot']))
        else:
            f.write('endpot endpot\n')
        if 'config' in kwargs:
            f.write('config {}\n'.format(kwargs['config']))
        else:
            f.write('config config\n')
        if 'tempfile' in kwargs:
            f.write('tempfile {}\n'.format(kwargs['tempfile']))
        else:
            f.write('tempfile tempfile\n')
        if 'output_prefix' in kwargs:
            f.write('output_prefix {}\n'.format(kwargs['output_prefix']))
        else:
            f.write('output_prefix output\n')
        for item in kwargs:
            f.write('{} {}\n'.format(item, kwargs[item]))
        f.close()
        return

    def create_potential_file(self, input):
        f = self.create_file('startpot')
        f.write(input)
        f.close()

    def call_makeapot(self, filename, args):
        os.environ['PATH'] = '{}:'.format(os.path.abspath('../util')) + os.environ['PATH']
        p = subprocess.Popen(['makeapot', '-o', filename] + args.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.cwd)
        p.wait()
        if p.returncode:
            pytest.fail('error calling "makeapot {}"'.format(args))
        self.filenames.append(filename)

    def create_config_file(self, **kwargs):
        f = self.create_file('config')
        if 'data' in kwargs:
            f.write(kwargs['data'])
            f.close()
            return
        self.config = simple_config(**kwargs)
        f.write(self.config.as_string())
        f.close()

    def run(self, param_file='param_file', args=[]):
        asan_filename = 'asan_{}'.format(''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6)))
        os.environ['ASAN_OPTIONS'] = 'log_path={},exitcode=99,strip_path_prefix={}'.format(asan_filename,os.path.abspath('..') + '/build/../')
        cmd = [os.path.join(os.path.abspath('../bin'), self.binary_name)]
        if len(args):
            cmd.extend(args)
        if param_file != None:
            cmd.append(param_file)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.cwd)
        p.wait()
        self.stdout = p.stdout.read().decode('ascii')
        self.stderr = p.stderr.read().decode('ascii')
        try:
            f = open(os.path.join(self.cwd, param_file), 'r')
            for line in f:
                for token in ['startpot', 'endpot']:
                  if token in line:
                      filename = os.path.join(self.cwd, line.split()[1])
                      if os.path.isfile(filename):
                          if not filename in self.filenames:
                              self.filenames.append(filename)
                          g = open(filename)
                          setattr(self, token, g.read())
                          g.close()
                if 'output_prefix' in line:
                    filenames = [x for x in glob.iglob(os.path.join(self.cwd, line.split()[1] + '*')) if os.path.isfile(x)]
                    for fname in filenames:
                        self.filenames.append(fname)
                        g = open(fname)
                        setattr(self, os.path.basename(fname).split('.')[1], g.read())
                        g.close()
            f.close()
        except:
            pass
        self.returncode = p.returncode
        self._check_asan_fail(asan_filename)

    def has_no_error(self):
        if self.returncode != 0:
            return False
        return not '[ERROR]' in self.stderr

    def has_no_warning(self):
        return not '[WARNING]' in self.stderr

    def has_error(self):
        if self.returncode == 0:
            return False
        return '[ERROR]' in self.stderr

    def has_correct_atom_count(self):
        return 'total of {} atoms'.format(self.config.atom_count()) in self.stdout

    def check_minimal_distance_matrix(self):
        self.config.check_matrix(self._get_minimal_distance_matrix())

    def clear(self):
        for item in list(set(self.filenames)):
            if os.path.isfile(item):
                os.remove(item)

    def _write_unit_cell(self, f, basic_size, ntypes, i, j, k):
        f.write('{} {} {} {} {} {} {}\n'.format(random.randint(0, ntypes - 1), i * basic_size, j * basic_size, k * basic_size, -random.random(), -random.random(), -random.random()))
        f.write('{} {} {} {} {} {} {}\n'.format(random.randint(0, ntypes - 1), (i + 0.5) * basic_size + random.uniform(-0.2,0.2), (j + 0.5) * basic_size + random.uniform(-0.2,0.2), (k + 0.5) * basic_size + random.uniform(-0.2,0.2), -random.random(), -random.random(), -random.random()))

    def _build(self, git_patch = None):
        p = subprocess.Popen(['./waf', 'distclean'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd='..')
        p.wait()
        if p.returncode:
            pytest.fail('error calling "waf distclean"')
        p = subprocess.Popen(self._get_conf_cmd(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd='..')
        p.wait()
        if p.returncode:
            print(p.stderr.read().decode('ascii'))
            pytest.fail('error calling "waf configure"')
        if git_patch:
          p = subprocess.Popen(['patch', '-Np1'], stdin=open(os.path.join(self.cwd, git_patch)), stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd='..')
          p.wait()
          if p.returncode:
            print(p.stderr.read().decode('ascii'))
            pytest.fail('error patching potfit source tree')
        p = subprocess.Popen(['./waf', 'build'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd='..')
        p.wait()
        if git_patch:
          q = subprocess.Popen(['git', 'checkout', '--', 'src/'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd='..')
          q.wait()
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

    def _get_minimal_distance_matrix(self):
        match = self.mindist_pattern.search(self.stdout)
        matrix = self.stdout[match.start():match.end()].strip().split('\n')[1:]
        data = []
        for i in range(len(matrix[0].split('\t')) - 2):
            data.append([float(s) for s in matrix[i + 1].strip().split('\t')[1:]])
        # check if symmetric
        if len(data) > 1:
            for i in range(len(data)):
                for j in range(i + 1, len(data)):
                    if data[i][j] != data[j][i]:
                        return None
        return data

class atom_positions():
    def __init__(self, ntypes, size, distance):
        self.iter = 0
        self.pos = []
        self.size = size
        self.distance = distance
        n = int(size / distance)
        idx = 0
        for i, j, k in product(range(n), range(n), range(n)):
            self.pos.append([idx,
                [distance * i, distance * j, distance * k],
                [random.uniform(-1,0) for _ in range(3)]])
            idx = (idx + 1) % ntypes

    def __iter__(self):
        return self

    def __next__(self):
        if self.iter >= len(self.pos):
            self.iter = 0
            raise StopIteration
        atom = self.pos[self.iter]
        self.iter += 1
        return atom

    def num_atoms(self):
        return len(self.pos)

    def mindist(self, i, j):
        min = -1
        for k, l in product(range(len(self.pos)), range(len(self.pos))):
            if self.pos[k][0] == i and self.pos[l][0] == j:
                d = self._mindist(k, l)
                if min == -1:
                    min = d
                elif d < min:
                    min = d
            if min == self.distance:
                break
        return min

    def _mindist(self, i, j):
        min = -1
        for k, l, m in product([-self.size, 0, self.size], [-self.size, 0, self.size], [-self.size, 0, self.size]):
            if (i == j) and (k == l == m == 0):
                continue
            vec = [k + self.pos[j][1][0] - self.pos[i][1][0], l + self.pos[j][1][1] - self.pos[i][1][1], m + self.pos[j][1][2] - self.pos[i][1][2]]
            if min > 0 and (vec[0] > min or vec[1] > min or vec[2] > min):
                continue
            d = math.sqrt(sum([a*a for a in vec]))
            if min == -1:
                min = d
            elif d < min:
                min = d
        return min

class simple_config():
    def __init__(self, **kwargs):
        seed = 42
        if 'seed' in kwargs:
            seed = int(kwargs['seed'])
        random.seed(seed)
        self.ntypes = 1
        if 'ntypes' in kwargs:
            self.ntypes = int(kwargs['ntypes'])
        self.size = 10
        if 'size' in kwargs:
            self.size = int(kwargs['size'])
        self.distance = 2.0
        if 'distance' in kwargs:
            self.distance = float(kwargs['distance'])
        self.useforce = 1
        if 'useforce' in kwargs:
            self.useforce = int(kwargs['useforce'])
        self.energy = -1.0
        if 'energy' in kwargs:
            self.energy = float(kwargs['energy'])
        self.contrib = 0
        if 'contrib' in kwargs:
            self.contrib = int(kwargs['contrib'])
        self.stress = None
        if 'stress' in kwargs:
            if int(kwargs['stress']):
                self.stress = ['{:2.6f}'.format(random.uniform(-1,1)) for _ in range(6)]
        self.weight = 1.0
        if 'weight' in kwargs:
            self.weight = float(kwargs['weight'])
        self.atoms = atom_positions(self.ntypes, self.size, self.distance)

    def as_string(self):
        config  = '#N {} {}\n'.format(self.atoms.num_atoms(), self.useforce)
        config += '#C {}\n'.format(' '.join([str(s) for s in range(self.ntypes)]))
        config += '## force file generated by potfit.py module\n'
        config += '#X {} 0 0\n'.format(self.size)
        config += '#Y 0 {} 0\n'.format(self.size)
        config += '#Z 0 0 {}\n'.format(self.size)
        if self.contrib:
           config += '#B_O {} {} {}\n'.format(0, 0, 0)
           config += '#B_A {} {} {}\n'.format(0.5 * self.size, 0, 0)
           config += '#B_B {} {} {}\n'.format(0, 0.5 * self.size, 0)
           config += '#B_C {} {} {}\n'.format(0, 0, 0.5 * self.size)
        if self.stress != None:
           config += '#S {}\n'.format(' '.join(self.stress))
        config += '#E 0\n'
        config += '#W {}\n'.format(self.weight)
        config += '#F\n'
        for atom in self.atoms:
            config += '{} {}'.format(atom[0], ' '.join(['{:2.6f}'.format(s) for s in atom[1]]))
            if self.useforce:
                config += ' {}\n'.format(' '.join(['{:2.6f}'.format(s) for s in atom[2]]))
            else:
                config += ' 0.0 0.0 0.0\n'
        return config

    def check_matrix(self, matrix):
        assert len(matrix) == self.ntypes
        for i, j in product(range(self.ntypes), range(self.ntypes)):
            assert math.isclose(matrix[i][j], self.atoms.mindist(i, j), rel_tol=1e-6)

    def atom_count(self):
        return self.atoms.num_atoms()
