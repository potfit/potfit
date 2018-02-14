import pytest

from itertools import product

POTENTIAL_BODY = '''0.00119911169735
7.20484717515e-06
0.000595394899272
0.000936445223073
0.000484087727962
0.000200648587256
0.0000959327293262
0.0000609377489204
0.000049458009818
0.0000457044379028
0.0000444847462971
0.0000440967362355
0.0000439803740159
0.0000439510638516
0.0'''

@pytest.mark.parametrize("types,distance", [c for c in product(range(1,11), range(19,22))])
def test_minimal_distances_tab(potfit_tab, types, distance):
    potfit_tab.create_param_file(ntypes=types)
    n = int(types * (types + 1) / 2)
    potfit_tab.create_potential_file('''
#F 3 {}
#I {}
#E

{}

{}
'''.format(n,
           ' '.join(['0' for _ in range(n)]),
           '\n'.join(['0.1 7 15' for _ in range(n)]),
           '\n'.join([POTENTIAL_BODY for _ in range(n)])))
    potfit_tab.create_config_file(ntypes=types, distance=distance / 10)
    potfit_tab.run()
    assert potfit_tab.has_no_error()
    potfit_tab.check_minimal_distance_matrix()

def test_minimal_distances_tab_outside(potfit_tab):
    potfit_tab.create_param_file(ntypes=3)
    n = 6
    potfit_tab.create_potential_file('''
#F 3 {}
#I {}
#E

{}

{}
'''.format(n,
           ' '.join(['0' for _ in range(n)]),
           '\n'.join(['0.1 7 15' for _ in range(n)]),
           '\n'.join([POTENTIAL_BODY for _ in range(n)])))
    potfit_tab.create_config_file(data='''
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
    potfit_tab.run()
    assert potfit_tab.has_error()
    assert 'No atoms found in interaction range for potential 0!' in potfit_tab.stderr

def test_minimal_distances_tab_short(potfit_tab):
    potfit_tab.create_param_file(ntypes=1)
    potfit_tab.create_potential_file('''
#F 3 1
#I 0
#E

2 7 15

{}
'''.format(POTENTIAL_BODY))
    potfit_tab.create_config_file(data='''
#N 2 0
#C H
#X 10 0 0
#Y 0 10 0
#Z 0 0 10
#E 0
#F
0 0 0.0 0.0 0 0 0
0 1 0.0 0.0 0 0 0
''')
    potfit_tab.run()
    assert potfit_tab.has_error()
    assert 'Configuration 1: Distance 1.000000' in potfit_tab.stderr
    assert 'The distance 1.000000 is smaller than the beginning' in potfit_tab.stderr
    assert 'of the potential #0 (r_begin=2.000000).' in potfit_tab.stderr
    assert 'Short distance!' in potfit_tab.stderr
