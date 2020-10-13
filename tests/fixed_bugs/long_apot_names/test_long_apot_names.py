import pytest

def test_fixed_bugs_long_apot_names(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

type thisisareallylongfunctionnameandhopefullyitwontcrashpotfit
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file(repeat_cell=3, seed=42)
    potfit.run()
    assert potfit.has_no_error()
