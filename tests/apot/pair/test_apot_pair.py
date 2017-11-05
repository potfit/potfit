import pytest

@pytest.fixture(scope="module")
def potfit():
    import sys
    sys.path.insert(0, str(pytest.config.rootdir))
    import potfit
    p = potfit.Potfit(__file__, 'apot', 'pair')
    yield p
    p.clear()

def test_apot_pair_basic(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 0 1\n#T PAIR\n#I 0\n#E\n\n')
    g.write('type lj\ncutoff 6.0\nepsilon 0.1 0 1\nsigma 2.5 1 4\n')
    g.close()
    c = potfit.create_config_file('config', repeat_cell=3, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 0
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'total of 54 atoms' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    # 162 from forces (3 * (3^2 * 2))
    # 1 from energy
    # 2 from potential parameter punishment
    # 1 from potential function punishment
    # 1 from general punishment (?)
    assert 'count 167' in potfit.stdout

def test_apot_pair_wrong_potential_format(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 3 1\n#T PAIR\n#I 0\n#E\n\n')
    g.write('type lj\ncutoff 6.0\nepsilon 0.1 0 1\nsigma 2.5 1 4\n')
    g.close()
    c = potfit.create_config_file('config', repeat_cell=3, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'This potfit binary only supports analytic potentials' in potfit.stderr

def test_apot_pair_wrong_number_of_potentials(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 0 3\n#T PAIR\n#I 0\n#E\n\n')
    g.write('type lj\ncutoff 6.0\nepsilon 0.1 0 1\nsigma 2.5 1 4\n')
    g.close()
    c = potfit.create_config_file('config', repeat_cell=3, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Wrong number of data columns in PAIR potential file' in potfit.stderr
    assert 'there should be 1, but there are 3' in potfit.stderr

def test_apot_pair_wrong_potential_type(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 0 1\n#T EAM\n#I 0\n#E\n\n')
    g.write('type lj\ncutoff 6.0\nepsilon 0.1 0 1\nsigma 2.5 1 4\n')
    g.close()
    c = potfit.create_config_file('config', repeat_cell=3, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Wrong potential type found in potential file!' in potfit.stderr
    assert 'This binary only supports PAIR potentials' in potfit.stderr

def test_apot_pair_wrong_number_of_invariant_entries(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 2\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 0 3\n#T PAIR\n#I 0 1\n#E\n\n')
    g.write('type lj\ncutoff 6.0\nepsilon 0.1 0 1\nsigma 2.5 1 4\n')
    g.close()
    c = potfit.create_config_file('config', repeat_cell=3, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Not enough items in #I header line' in potfit.stderr

def test_apot_pair_no_header_end_line(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 0 1\n#T PAIR\n#I 0\n\n')
    g.write('type lj\ncutoff 6.0\nepsilon 0.1 0 1\nsigma 2.5 1 4\n')
    g.close()
    c = potfit.create_config_file('config', repeat_cell=3, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Header corrupt in file startpot' in potfit.stderr

def test_apot_pair_unknown_potential_function(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 0 1\n#T PAIR\n#I 0\n#E\n\n')
    g.write('type ljs\ncutoff 6.0\nepsilon 0.1 0 1\nsigma 2.5 1 4\n')
    g.close()
    c = potfit.create_config_file('config', repeat_cell=3, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Unknown function type in file startpot' in potfit.stderr
    assert 'please define "ljs" in functions.c' in potfit.stderr

def test_apot_pair_no_potential_function(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 0 1\n#T PAIR\n#I 0\n#E\n\n')
    g.close()
    c = potfit.create_config_file('config', repeat_cell=3, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Error while searching for analytic potentials' in potfit.stderr

def test_apot_pair_invalid_keyword(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 0 1\n#T PAIR\n#I 0\n#E\n\n')
    g.write('type lj\ncutoffs 6.0\nepsilon 0.1 0 1\nsigma 2.5 1 4\n')
    g.close()
    c = potfit.create_config_file('config', repeat_cell=3, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'No cutoff found for the 1. potential' in potfit.stderr

def test_apot_pair_not_enough_params(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 0 1\n#T PAIR\n#I 0\n#E\n\n')
    g.write('type lj\ncutoff 6.0\nepsilon 0.1 0 1\n')
    g.close()
    c = potfit.create_config_file('config', repeat_cell=3, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Error reading analytic potentials' in potfit.stderr

def test_apot_pair_not_enough_param_values(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 0 1\n#T PAIR\n#I 0\n#E\n\n')
    g.write('type lj\ncutoff 6.0\nepsilon 0.1 0\nsigma 2.5 1 4\n')
    g.close()
    c = potfit.create_config_file('config', repeat_cell=3, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Could not read parameter #1 of potential #1 in file startpot' in potfit.stderr


def test_apot_pair_not_enough_potentials(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 2\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 0 3\n#T PAIR\n#I 0 0 0\n#E\n\n')
    g.write('type lj\ncutoff 6.0\nepsilon 0.1 0 1\nsigma 2.5 1 4\n')
    g.close()
    c = potfit.create_config_file('config', repeat_cell=3, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Error while searching for analytic potentials' in potfit.stderr
