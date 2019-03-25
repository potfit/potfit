import pytest

def test_apot_pair_invariant_base(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 1
#E

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
    assert '1 PAIR potential(s)' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'Optimization disabled due to 0 free parameters' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert '378 contributions' in potfit.stdout
