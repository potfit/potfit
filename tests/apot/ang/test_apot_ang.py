import pytest

def test_apot_ang_basic(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 3
#T ANG
#I 0 0 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file(size=6)
    potfit.run()
    assert potfit.has_no_error()
    assert 'analytic potentials' in potfit.stdout
    assert '3 ANG potential(s)' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'total of 27 atoms' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert '92 contributions' in potfit.stdout
