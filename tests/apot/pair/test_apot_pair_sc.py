import pytest

@pytest.fixture(scope="module")
def potfit():
    import sys
    sys.path.insert(0, str(pytest.config.rootdir))
    import potfit
    p = potfit.Potfit(__file__, 'apot', 'pair')
    yield p
    p.clear()

def test_apot_pair_sc(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    potfit.call_makeapot('startpot', '-n 1 -i pair -f eopp_sc')
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
    # 7 from potential parameter punishment
    # 1 from potential function punishment
    # 1 from general punishment (?)
    assert 'count 172' in potfit.stdout

def test_apot_pair_global_sc(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    potfit.call_makeapot('startpot', '-n 1 -i pair -f eopp_sc -g')
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
    # 7 from potential parameter punishment
    # 1 from potential function punishment
    # 1 from general punishment (?)
    assert 'count 172' in potfit.stdout

def test_apot_pair_sc_param_missing(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 0 1\n#T PAIR\n#I 0 0 0\n#E\n\n')
    g.write('type lj_sc\ncutoff 6.0\nepsilon 0.1 0 1\nsigma 2.5 1 4\n')
    g.close()
    c = potfit.create_config_file('config', repeat_cell=3, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Error reading analytic potentials' in potfit.stderr
