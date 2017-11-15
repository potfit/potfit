import pytest

def test_apot_pair_stress_missing_stress_weight(potfit):
    potfit.create_param_file()
    potfit.run()
    assert potfit.has_error()
    assert 'stress_weight' in potfit.stderr

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
    potfit.create_config_file(repeat_cell=3, seed=42, stress=True)
    potfit.run()
    assert potfit.has_no_error()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Read 1 configurations (1 with forces, 1 with stresses)' in potfit.stdout
    assert 'total of 54 atoms' in potfit.stdout
    assert 'Global stress weight: 1.000000' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert 'global stress weight w is 1.00' in potfit.stress
    assert 'count 173' in potfit.stdout
    assert 'sum of stress-errors' in potfit.stdout
