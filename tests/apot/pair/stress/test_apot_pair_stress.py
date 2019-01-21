import pytest

def test_apot_pair_stress_weight_missing(potfit):
    potfit.create_param_file()
    potfit.run()
    assert potfit.has_error()
    assert 'stress_weight' in potfit.stderr

def test_apot_pair_stress_weight_empty(potfit):
    potfit.create_param_file(stress_weight='')
    potfit.call_makeapot('startpot', '-n 1 -i pair -f eopp_sc')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Missing value in parameter file' in potfit.stderr
    assert 'stress_weight is <undefined>' in potfit.stderr

def test_apot_pair_stress_weight_invalid(potfit):
    potfit.create_param_file(stress_weight='foo')
    potfit.call_makeapot('startpot', '-n 1 -i pair -f eopp_sc')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Illegal value in parameter file' in potfit.stderr
    assert 'stress_weight is not a double' in potfit.stderr
    assert 'value = foo' in potfit.stderr

def test_apot_pair_stress_weight_out_of_bounds  (potfit):
    potfit.create_param_file(stress_weight=-1)
    potfit.call_makeapot('startpot', '-n 1 -i pair -f eopp_sc')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Illegal value in parameter file' in potfit.stderr
    assert 'stress_weight is out of bounds' in potfit.stderr

def test_apot_pair_stress_basic(potfit):
    potfit.create_param_file(stress_weight=1)
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
    potfit.create_config_file(stress=1)
    potfit.run()
    assert potfit.has_no_error()
    assert potfit.has_correct_atom_count()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Read 1 configurations (1 with forces, 1 with stresses)' in potfit.stdout
    assert 'Global stress weight: 1.000000' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert 'global stress weight w is 1.00' in potfit.stress
    assert '386 contributions' in potfit.stdout
    assert 'sum of stress-errors' in potfit.stdout
