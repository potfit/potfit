import pytest

def test_kim_basic_missing_model(potfit):
    potfit.create_param_file()
    potfit.run()
    assert(potfit.has_error())
    assert('Missing parameter: kim_model_name' in potfit.stderr)

@pytest.mark.skip(reason="Currently there is a memory leak in OpenKIM")
def test_kim_basic_invalid_model(potfit):
    potfit.create_param_file(kim_model_name='DOES_NOT_EXIST')
    potfit.run()
    assert(potfit.has_error())
    assert('KIM_Model_Create failed: 1' in potfit.stderr)

def test_kim_basic_wrong_potential_format(potfit):
    potfit.create_param_file(kim_model_name='ex_model_Ar_P_LJ')
    potfit.create_potential_file('''
#F 0 1
''')
    potfit.run()
    assert(potfit.has_error())
    assert('This potfit binary only supports KIM potentials' in potfit.stderr)

def test_kim_basic_model_mismatch(potfit):
    potfit.create_param_file(kim_model_name='ex_model_Ar_P_LJ')
    potfit.create_potential_file('''
#F 5 1
#T DOES_NOT_MATCH
''')
    potfit.run()
    assert(potfit.has_error())
    assert('Potential mismatch detected!' in potfit.stderr)
    assert('KIM potential selected is: ex_model_Ar_P_LJ' in potfit.stderr)
    assert('The potential file contains a DOES_NOT_MATCH potential.' in potfit.stderr)

def test_kim_basic_model_missing(potfit):
    potfit.create_param_file(kim_model_name='ex_model_Ar_P_LJ')
    potfit.create_potential_file('''
#F 5 1
#T ex_model_Ar_P_LJ
''')
    potfit.run()
    assert(potfit.has_error())
    assert('Unexpected end of file in startpot' in potfit.stderr)
