import pytest

#def get_installed_kim_models():
    #return ['LennardJones612_UniversalShifted__MO_959249795837_003']

#@pytest.mark.parametrize("models", get_installed_kim_models())
#def test_kim_list_model_params(potfit, models):
    #potfit.create_param_file()
    #potfit.create_potential_file('''
##F 5 1
##C Ar
##E

#model {}
#kim_list_model_params
#'''.format(models))
    #potfit.create_config_file()
    #potfit.run()
    #print(potfit.stdout)
    #print(potfit.stderr)
    #assert('kim_list_model_params" has been found in the potential file!' in potfit.stderr)
