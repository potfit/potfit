import pytest

POTENTIAL_FILE_CONTENT = '''#F 5 1
#T ex_model_Ar_P_LJ
#C Ar
#E

KIM_PARAM cutoff
8.150000 8.150000 8.150000

KIM_PARAM epsilon
0.010400 0.010400 0.010400

KIM_PARAM sigma
3.400000 3.400000 3.400000
'''

def test_kim_model_params_dump_file(potfit):
    potfit.create_param_file(kim_model_name='ex_model_Ar_P_LJ', kim_model_params='dump_file')
    potfit.create_potential_file('''
#F 5 1
#T ex_model_Ar_P_LJ
#C Ar
#E
''')
    potfit.run()
    assert(POTENTIAL_FILE_CONTENT in potfit.get_file_content('ex_model_Ar_P_LJ.default'))
    assert(potfit.has_no_error())

def test_kim_model_params_dump(potfit):
    potfit.create_param_file(kim_model_name='ex_model_Ar_P_LJ', kim_model_params='dump')
    potfit.create_potential_file('''
#F 5 1
#T ex_model_Ar_P_LJ
#C Ar
#E
''')
    potfit.run()
    assert(POTENTIAL_FILE_CONTENT in potfit.stdout)
    assert(potfit.has_no_error())

def test_kim_model_params_use_default(potfit):
    potfit.create_param_file(kim_model_name='ex_model_Ar_P_LJ', kim_model_params='use_default')
    potfit.create_potential_file('''
#F 5 1
#T ex_model_Ar_P_LJ
#C Ar
#E''')
    potfit.create_config_file(elements=['Ar'])
    potfit.run()
    assert(potfit.has_no_error())
