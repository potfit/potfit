import pytest

@pytest.fixture(scope="module")
def potfit():
    import sys
    sys.path.insert(0, str(pytest.config.rootdir))
    import potfit
    p = potfit.Potfit(__file__, 'tab', 'pair')
    yield p
    p.clear()

def test_tab_pair_no_features_simple(potfit):
    param_file = 'param_file'
    f = potfit.create_file(param_file)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.close()
    g = potfit.create_file('startpot')
    g.write('#F 3 1\n')
    g.write('#I 0\n')
    g.write('#E\n\n')
    g.write('1.158337 10.883950976 15\n\n')
    g.write('0.00119911169735\n7.20484717515e-06\n0.000595394899272\n0.000936445223073\n0.000484087727962\n')
    g.write('0.000200648587256\n0.0000959327293262\n0.0000609377489204\n0.000049458009818\n0.0000457044379028\n')
    g.write('0.0000444847462971\n0.0000440967362355\n0.0000439803740159\n0.0000439510638516\n0.0\n')
    g.close()
    c = potfit.create_config_file('config', 4, seed=42)
    c.close()
    potfit.run(param_file)
    assert potfit.returncode == 0
    assert 'tabulated eqdist' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'total of 4 atoms' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 3 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert 'total error sum 0.000000' in potfit.stdout
    potfit.cleanup(param_file)
