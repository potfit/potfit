import copy
import pytest

@pytest.mark.parametrize("ntypes,elements", [[1, ['Na']], [2, ['Na', 'Cl']], [3, ['Na', 'Cl', 'H']]])
def test_empty_element_in_config(potfit_apot, ntypes, elements):
    potfit_apot.create_param_file(ntypes=ntypes)
    potfit_apot.call_makeapot('startpot', '-n {} -i pair -f {}*lj'.format(ntypes, int(0.5 * ntypes * (ntypes + 1))))
    potfit_apot.create_config_file(ntypes=ntypes, elements=elements)
    potfit_apot.run()
    assert potfit_apot.has_no_error()
    assert potfit_apot.has_correct_atom_count()
    assert '#C {}'.format(' '.join(elements)) in potfit_apot.endpot

@pytest.mark.parametrize("ntypes,elements", [[1, ['Sodium']], [2, ['Sodium', 'Magnesium']], [3, ['Sodium', 'Magnesium', 'Kalium']]])
def test_too_long_element_name(potfit_apot, ntypes, elements):
    potfit_apot.create_param_file(ntypes=ntypes)
    potfit_apot.call_makeapot('startpot', '-n {} -i pair -f {}*lj'.format(ntypes, int(0.5 * ntypes * (ntypes + 1))))
    potfit_apot.create_config_file(ntypes=ntypes, elements=elements)
    potfit_apot.run()
    assert potfit_apot.has_no_error()
    assert potfit_apot.has_correct_atom_count()
    assert '#C {}'.format(' '.join(elements)) in potfit_apot.endpot

@pytest.mark.parametrize("ntypes", [1, 2, 3])
def test_numbered_elements(potfit_apot, ntypes):
    potfit_apot.create_param_file(ntypes=ntypes)
    potfit_apot.call_makeapot('startpot', '-n {} -i pair -f {}*lj'.format(ntypes, int(0.5 * ntypes * (ntypes + 1))))
    potfit_apot.create_config_file(ntypes=ntypes)
    potfit_apot.run()
    assert potfit_apot.has_no_error()
    assert potfit_apot.has_correct_atom_count()
    assert '#C {}'.format(' '.join([str(s) for s in range(ntypes)])) in potfit_apot.endpot

manual_config = '''
#N 2 0
#C H
#X 10 0 0
#Y 0 10 0
#Z 0 0 10
#E 0
#F
0 0 0 0 0 0 0
0 3 0 0 0 0 0
#N 4 0
#C H He
#X 10 0 0
#Y 0 10 0
#Z 0 0 10
#E 0
#F
0 0 0 0 0 0 0
0 3 0 0 0 0 0
1 0 3 0 0 0 0
1 0 0 3 0 0 0
#N 6 0
#C H {} Li
#X 10 0 0
#Y 0 10 0
#Z 0 0 10
#E 0
#F
0 0 0 0 0 0 0
0 3 0 0 0 0 0
1 0 3 0 0 0 0
1 0 0 3 0 0 0
2 3 3 0 0 0 0
2 3 0 3 0 0 0
'''

def test_element_name_config_good(potfit_apot):
    potfit_apot.create_param_file(ntypes=3)
    potfit_apot.call_makeapot('startpot', '-n 3 -i pair -f 6*lj')
    potfit_apot.create_config_file(data=manual_config.format('He'))
    potfit_apot.run()
    assert potfit_apot.has_no_error()
    assert '#C H He Li' in potfit_apot.endpot

def test_element_name_config_bad(potfit_apot):
    potfit_apot.create_param_file(ntypes=3)
    potfit_apot.call_makeapot('startpot', '-n 3 -i pair -f 6*lj')
    potfit_apot.create_config_file(data=manual_config.format('Ha'))
    potfit_apot.run()
    assert potfit_apot.has_error()
    assert 'Expected element >> He << but found element >> Ha <<' in potfit_apot.stderr

@pytest.mark.parametrize("ntypes,elements", [[1, ['Na']], [2, ['Na', 'Cl']], [3, ['Na', 'Cl', 'H']]])
def test_elements_in_config_and_potential(potfit_apot, ntypes, elements):
    potfit_apot.create_param_file(ntypes=ntypes)
    potfit_apot.call_makeapot('startpot', '-n {} -i pair -e {} -f {}*lj'.format(ntypes, ','.join(elements) ,int(0.5 * ntypes * (ntypes + 1))))
    potfit_apot.create_config_file(ntypes=ntypes, elements=elements)
    potfit_apot.run()
    assert potfit_apot.has_no_error()
    assert potfit_apot.has_correct_atom_count()
    assert '#C {}'.format(' '.join(elements)) in potfit_apot.endpot

@pytest.mark.parametrize("ntypes,elements", [[1, ['Na']], [2, ['Na', 'Cl']], [3, ['Na', 'Cl', 'He']]])
def test_elements_in_config_and_potential_mismatch(potfit_apot, ntypes, elements):
    potfit_apot.create_param_file(ntypes=ntypes)
    wrong_elements = copy.copy(elements)
    wrong_elements[-1] = elements[-1][::-1]
    potfit_apot.call_makeapot('startpot', '-n {} -i pair -e {} -f {}*lj'.format(ntypes, ','.join(wrong_elements) ,int(0.5 * ntypes * (ntypes + 1))))
    potfit_apot.create_config_file(ntypes=ntypes, elements=elements)
    potfit_apot.run()
    assert potfit_apot.has_error()
    assert 'Expected element >> {} << but found element >> {} <<'.format(wrong_elements[-1], elements[-1]) in potfit_apot.stderr

def test_elements_in_config_and_potential_mismatch_pt2(potfit_apot):
    potfit_apot.create_param_file(ntypes=2)
    potfit_apot.call_makeapot('startpot', '-n 2 -i pair -f 3*lj -e He,Se')
    potfit_apot.create_config_file(data='''
#N 6 0
#C Se He
#X 10 0 0
#Y 0 10 0
#Z 0 0 10
#E 0
#F
0 0 0 0 0 0 0
0 3 0 0 0 0 0
1 0 3 0 0 0 0
1 0 0 3 0 0 0
0 3 3 0 0 0 0
1 3 0 3 0 0 0
''')
    potfit_apot.run()
    assert potfit_apot.has_error()
    assert 'Expected element >> He << but found element >> Se <<' in potfit_apot.stderr

