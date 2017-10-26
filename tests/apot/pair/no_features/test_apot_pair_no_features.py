import pytest

@pytest.fixture(scope="module")
def potfit():
    import sys
    sys.path.insert(0, str(pytest.config.rootdir))
    import potfit
    p = potfit.Potfit(__file__, 'apot', 'pair')
    yield p
    p.clear()

def test_apot_pair_no_features_simple(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
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
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert '[ERROR] Could not open file config' in potfit.stderr
