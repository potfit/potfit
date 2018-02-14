import pytest

from itertools import product

@pytest.mark.parametrize("types,distance", [c for c in product(range(1,11), range(19,22))])
def test_minimal_distances_apot(potfit_apot, types, distance):
    potfit_apot.create_param_file(ntypes=types)
    potfit_apot.call_makeapot('startpot', '-n {} -i pair -f {}*lj'.format(types, int(types*(types+1)/2)))
    potfit_apot.create_config_file(ntypes=types, distance=distance / 10)
    potfit_apot.run()
    assert potfit_apot.has_no_error()
    potfit_apot.check_minimal_distance_matrix()

def test_minimal_distances_apot_default_value(potfit_apot):
    potfit_apot.create_param_file(ntypes=3)
    potfit_apot.call_makeapot('startpot', '-n 3 -i pair -f 6*lj')
    potfit_apot.create_config_file(data='''
#N 6 0
#C H He Li
#X 100 0 0
#Y 0 100 0
#Z 0 0 100
#E 0
#F
0 0 0 0 0 0 0
0 50 0 0 0 0 0
1 0 50 0 0 0 0
1 0 0 50 0 0 0
2 50 50 50 0 0 0
2 75 75 75 0 0 0
''')
    potfit_apot.run()
    assert potfit_apot.has_error()
    assert 'No atoms found in interaction range for potential 0!' in potfit_apot.stderr

def test_minimal_distances_apot_short(potfit_apot):
    potfit_apot.create_param_file(ntypes=1)
    potfit_apot.call_makeapot('startpot', '-n 1 -i pair -f lj')
    potfit_apot.create_config_file(data='''
#N 2 0
#C H
#X 10 0 0
#Y 0 10 0
#Z 0 0 10
#E 0
#F
0 0.000000 0.0 0.0 0 0 0
0 0.000001 0.0 0.0 0 0 0
''')
    potfit_apot.run()
    assert potfit_apot.has_error()
    assert 'Configuration 1: Distance 0.000001' in potfit_apot.stderr
    assert 'The distance 0.000001 is smaller than the beginning' in potfit_apot.stderr
    assert 'of the potential #0 (r_begin=0.000100).' in potfit_apot.stderr
    assert 'Short distance!' in potfit_apot.stderr
