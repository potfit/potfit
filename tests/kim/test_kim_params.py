import pytest

def test_kim_params_invalid(potfit):
    potfit.create_param_file(kim_model_name='ex_model_Ar_P_LJ')
    potfit.create_potential_file('''
#F 5 1
#T ex_model_Ar_P_LJ
#C Ar
#E

KIM_CUTOFF invalid
8.150000        8.150000        8.150000

KIM_PARAM epsilon
0.010400        0.010400        0.010400

KIM_PARAM sigma
3.400000        3.400000        3.400000
''')
    potfit.create_config_file(elements=['Ar'])
    potfit.run()
    assert(potfit.has_error())
    assert('Could not associate KIM_CUTOFF name "invalid" with any parameter!' in potfit.stderr)

def test_kim_params_no_cutoff(potfit):
    potfit.create_param_file(kim_model_name='ex_model_Ar_P_LJ')
    potfit.create_potential_file('''
#F 5 1
#T ex_model_Ar_P_LJ
#C Ar
#E

KIM_PARAM invalid
8.150000        8.150000        8.150000

KIM_PARAM epsilon
0.010400        0.010400        0.010400

KIM_PARAM sigma
3.400000        3.400000        3.400000
''')
    potfit.create_config_file(elements=['Ar'])
    potfit.run()
    assert(potfit.has_error())
    assert('Parameter order mismatch, expected "cutoff" but found "invalid"' in potfit.stderr)
    assert('No cutoff parameter specified! All parameters will be optimized!' in potfit.stderr)
