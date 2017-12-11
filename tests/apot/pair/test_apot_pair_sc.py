import pytest

def test_apot_pair_sc(potfit):
    potfit.create_param_file()
    potfit.call_makeapot('startpot', '-n 1 -i pair -f eopp_sc')
    potfit.create_config_file(repeat_cell=3, seed=42)
    potfit.run()
    assert potfit.has_no_error()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'total of 54 atoms' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert 'count 172' in potfit.stdout

def test_apot_pair_global_sc(potfit):
    potfit.create_param_file()
    potfit.call_makeapot('startpot', '-n 1 -i pair -f eopp_sc -g')
    potfit.create_config_file(repeat_cell=3, seed=42)
    potfit.run()
    assert potfit.has_no_error()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'total of 54 atoms' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert 'count 172' in potfit.stdout

def test_apot_pair_sc_param_missing(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0 0 0
#E

type lj_sc
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file(repeat_cell=3, seed=42)
    potfit.run()
    assert potfit.has_error()
    assert 'Error reading analytic potentials' in potfit.stderr
