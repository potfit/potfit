import pytest

def test_apot_pair_globals(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

global 1
g 1 0 10

type lj_sc
cutoff 6.0
epsilon 0.1 0.1 .1
sigma 2.5 2.5 2.5
g!
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_no_error()
    assert potfit.has_correct_atom_count()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'Read 1 global parameter(s)' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert '379 contributions' in potfit.stdout

def test_apot_pair_globals_not_used(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

global 1
g 1 0 10

type lj_sc
cutoff 6.0
epsilon 0.1 0.1 .1
sigma 2.5 2.5 2.5
g 1 0 2
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_no_error()
    assert potfit.has_correct_atom_count()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'Read 1 global parameter(s)' in potfit.stdout
    assert 'You defined global parameters but did not use them.' in potfit.stdout
    assert 'Disabling global parameters.' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert '379 contributions' in potfit.stdout

def test_apot_pair_globals_not_defined(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

type lj_sc
cutoff 6.0
epsilon 0.1 0.1 .1
sigma 2.5 2.5 2.5
g!
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Could not find global parameter g' in potfit.stderr
