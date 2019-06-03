import pytest

from itertools import product

def test_overlapping_atoms(potfit_apot):
    potfit_apot.create_param_file(ntypes=1)
    potfit_apot.call_makeapot('startpot', '-n 1 -i pair -f 1*lj')
    potfit_apot.create_config_file(data='''
#N 2 0
#C H
#X 1 0 0
#Y 0 1 0
#Z 0 0 1
#E 0
#F
0 0 0 0 0 0 0
0 1 0 0 0 0 0
''')
    potfit_apot.run()
    assert potfit_apot.has_error()
    assert 'Overlapping atoms found in configuration 1!' in potfit_apot.stderr
    assert 'Atom 0 @ (0.000000, 0.000000, 0.000000)' in potfit_apot.stderr
    assert 'overlaps with atom 1 @ (1.000000, 0.000000, 0.000000)' in potfit_apot.stderr
    assert 'in this periodic copy of the unit cell: x=-1, y=0, z=0' in potfit_apot.stderr
