import pytest

def test_apot_pair_cp(potfit):
    potfit.create_param_file(enable_cp=1)
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

cp_0 1 2 3

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_no_error()
    assert potfit.has_correct_atom_count()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Enabled 1 chemical potential(s)' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert 'count 381' in potfit.stdout

def test_apot_pair_cp_missing(potfit):
    potfit.create_param_file(enable_cp=1)
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Error while searching for chemical potentials' in potfit.stderr

def test_apot_pair_cp_parameter_missing(potfit):
    potfit.create_param_file(enable_cp=1)
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

cp_0 1 2

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Could not read chemical potential for 0. atomtype' in potfit.stderr

def test_apot_pair_cp_invalid(potfit):
    potfit.create_param_file(enable_cp=1)
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

cp0 1 0 2

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'No chemical potentials found in startpot' in potfit.stderr
