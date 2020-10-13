import pytest

def test_kim_elements_not_supported(potfit):
    potfit.create_param_file(kim_model_name='ex_model_Ar_P_LJ')
    potfit.create_potential_file('''
#F 5 1
#T ex_model_Ar_P_LJ
#C He
#E

KIM_CUTOFF cutoff
8.150000        8.150000        8.150000

KIM_PARAM epsilon
0.010400        0.010400        0.010400

KIM_PARAM sigma
3.400000        3.400000        3.400000
''')
    potfit.create_config_file(elements=['He'])
    potfit.run()
    assert(potfit.has_error())
    assert('Species He not supported by the KIM model' in potfit.stderr)

def test_kim_elements_not_mismatch(potfit):
    potfit.create_param_file(kim_model_name='ex_model_Ar_P_LJ')
    potfit.create_potential_file('''
#F 5 1
#T ex_model_Ar_P_LJ
#C He
#E

KIM_CUTOFF cutoff
8.150000        8.150000        8.150000

KIM_PARAM epsilon
0.010400        0.010400        0.010400

KIM_PARAM sigma
3.400000        3.400000        3.400000
''')
    potfit.create_config_file(elements=['Fe'])
    potfit.run()
    assert(potfit.has_error())
    assert('Expected element >> He << but found element >> Fe <<.' in potfit.stderr)
