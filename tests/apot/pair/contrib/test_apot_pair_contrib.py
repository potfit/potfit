import pytest

def test_apot_pair_contrib_basic_box(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file(contrib=1)
    potfit.run()
    assert potfit.has_no_error()
    assert potfit.has_correct_atom_count()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert 'count 380' in potfit.stdout

def test_apot_pair_contrib_basic_sphere(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    f = potfit.create_file('config')
    f.write('''
#N 4 1
#C Al
#X    4.04836934    0.00000000    0.00000000
#Y    0.00000000    4.04836934    0.00000000
#Z    0.00000000    0.00000000    4.04836934
#B_S 0 0 0 1
#W 50
#E -3.6884100000
#F
0           0           0           0           0           0           0
0           0     2.02418     2.02418           0           0           0
0     2.02418           0     2.02418           0           0           0
0     2.02418     2.02418           0           0           0           0''')
    f.close()
    potfit.run()
    assert potfit.has_no_error()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'total of 4 atoms' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert 'count 17' in potfit.stdout

def test_apot_pair_contrib_basic_multi_sphere(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    f = potfit.create_file('config')
    f.write('''
#N 4 1
#C Al
#X    4.04836934    0.00000000    0.00000000
#Y    0.00000000    4.04836934    0.00000000
#Z    0.00000000    0.00000000    4.04836934
#B_S 0 0 0 1
#B_S 0 1 0 1
#B_S 0 0 1 1
#B_S 1 1 0 1
#B_S 0 1 1 1
#W 50
#E -3.6884100000
#F
0           0           0           0           0           0           0
0           0     2.02418     2.02418           0           0           0
0     2.02418           0     2.02418           0           0           0
0     2.02418     2.02418           0           0           0           0''')
    f.close()
    potfit.run()
    assert potfit.has_no_error()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'total of 4 atoms' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert 'count 17' in potfit.stdout

def test_apot_pair_contrib_missing_origin(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    f = potfit.create_file('config')
    f.write('''
#N 4 1
#C Al
#X    4.04836934    0.00000000    0.00000000
#Y    0.00000000    4.04836934    0.00000000
#Z    0.00000000    0.00000000    4.04836934
#B_A 1 0 0
#B_B 0 1 0
#B_C 0 0 1
#W 50
#E -3.6884100000
#F
0           0           0           0           0           0           0
0           0     2.02418     2.02418           0           0           0
0     2.02418           0     2.02418           0           0           0
0     2.02418     2.02418           0           0           0           0''')
    f.close()
    potfit.run()
    assert potfit.has_error()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Incomplete box of contributing atoms for config 1!' in potfit.stderr

def test_apot_pair_contrib_missing_origin_value(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    f = potfit.create_file('config')
    f.write('''
#N 4 1
#C Al
#X    4.04836934    0.00000000    0.00000000
#Y    0.00000000    4.04836934    0.00000000
#Z    0.00000000    0.00000000    4.04836934
#B_O 1 1
#B_A 1 0 0
#B_B 0 1 0
#B_C 0 0 1
#W 50
#E -3.6884100000
#F
0           0           0           0           0           0           0
0           0     2.02418     2.02418           0           0           0
0     2.02418           0     2.02418           0           0           0
0     2.02418     2.02418           0           0           0           0''')
    f.close()
    potfit.run()
    assert potfit.has_error()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Error reading box vector #B_O' in potfit.stderr

def test_apot_pair_contrib_duplicate_origin(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    f = potfit.create_file('config')
    f.write('''
#N 4 1
#C Al
#X    4.04836934    0.00000000    0.00000000
#Y    0.00000000    4.04836934    0.00000000
#Z    0.00000000    0.00000000    4.04836934
#B_O 0 0 0
#B_A 1 0 0
#B_B 0 1 0
#B_C 0 0 1
#B_O 0 0 0
#W 50
#E -3.6884100000
#F
0           0           0           0           0           0           0
0           0     2.02418     2.02418           0           0           0
0     2.02418           0     2.02418           0           0           0
0     2.02418     2.02418           0           0           0           0''')
    f.close()
    potfit.run()
    assert potfit.has_error()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'There can only be one box of contributing atoms' in potfit.stderr

def test_apot_pair_contrib_missing_vector_A(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    f = potfit.create_file('config')
    f.write('''
#N 4 1
#C Al
#X    4.04836934    0.00000000    0.00000000
#Y    0.00000000    4.04836934    0.00000000
#Z    0.00000000    0.00000000    4.04836934
#B_O 0 0 0
#B_B 0 1 0
#B_C 0 0 1
#W 50
#E -3.6884100000
#F
0           0           0           0           0           0           0
0           0     2.02418     2.02418           0           0           0
0     2.02418           0     2.02418           0           0           0
0     2.02418     2.02418           0           0           0           0''')
    f.close()
    potfit.run()
    assert potfit.has_error()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Incomplete box of contributing atoms for config 1!' in potfit.stderr

def test_apot_pair_contrib_missing_vector_B(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    f = potfit.create_file('config')
    f.write('''
#N 4 1
#C Al
#X    4.04836934    0.00000000    0.00000000
#Y    0.00000000    4.04836934    0.00000000
#Z    0.00000000    0.00000000    4.04836934
#B_O 0 0 0
#B_A 0 1 0
#B_C 0 0 1
#W 50
#E -3.6884100000
#F
0           0           0           0           0           0           0
0           0     2.02418     2.02418           0           0           0
0     2.02418           0     2.02418           0           0           0
0     2.02418     2.02418           0           0           0           0''')
    f.close()
    potfit.run()
    assert potfit.has_error()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Incomplete box of contributing atoms for config 1!' in potfit.stderr

def test_apot_pair_contrib_missing_vector_C(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    f = potfit.create_file('config')
    f.write('''
#N 4 1
#C Al
#X    4.04836934    0.00000000    0.00000000
#Y    0.00000000    4.04836934    0.00000000
#Z    0.00000000    0.00000000    4.04836934
#B_O 0 0 0
#B_A 0 1 0
#B_B 0 0 1
#W 50
#E -3.6884100000
#F
0           0           0           0           0           0           0
0           0     2.02418     2.02418           0           0           0
0     2.02418           0     2.02418           0           0           0
0     2.02418     2.02418           0           0           0           0''')
    f.close()
    potfit.run()
    assert potfit.has_error()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Incomplete box of contributing atoms for config 1!' in potfit.stderr

def test_apot_pair_contrib_missing_sphere_value(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 1
#T PAIR
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    f = potfit.create_file('config')
    f.write('''
#N 4 1
#C Al
#X    4.04836934    0.00000000    0.00000000
#Y    0.00000000    4.04836934    0.00000000
#Z    0.00000000    0.00000000    4.04836934
#B_S 0 0 0
#W 50
#E -3.6884100000
#F
0           0           0           0           0           0           0
0           0     2.02418     2.02418           0           0           0
0     2.02418           0     2.02418           0           0           0
0     2.02418     2.02418           0           0           0           0''')
    f.close()
    potfit.run()
    assert potfit.has_error()
    assert 'analytic potentials' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Error reading sphere of contributing atoms' in potfit.stderr
