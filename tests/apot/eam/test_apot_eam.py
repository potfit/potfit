import pytest

def test_apot_eam_basic(potfit):
    potfit.create_param_file()
    potfit.call_makeapot('startpot', '-n 1 -i eam -f 3*lj')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_no_error()
    assert potfit.has_correct_atom_count()
    assert 'analytic potentials' in potfit.stdout
    assert '3 EAM potentials' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 0 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert 'count 389' in potfit.stdout

def test_apot_eam_wrong_potential_format(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 3 1
#T EAM
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'This potfit binary only supports analytic potentials' in potfit.stderr

def test_apot_eam_wrong_number_of_potentials(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 2
#T PAIR
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Wrong number of data columns in EAM potential file' in potfit.stderr
    assert 'there should be 3, but there are 2' in potfit.stderr

def test_apot_eam_wrong_potential_type(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 3
#T PAIR
#I 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Wrong potential type found in potential file!' in potfit.stderr
    assert 'This binary only supports EAM potentials' in potfit.stderr

def test_apot_eam_wrong_number_of_invariant_entries(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 3
#T EAM
#I 0 1
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Not enough items in #I header line' in potfit.stderr

def test_apot_eam_no_header_end_line(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 3
#T EAM
#I 0 0 0

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Header corrupt in file startpot' in potfit.stderr

def test_apot_eam_unknown_potential_function(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 3
#T EAM
#I 0 0 0
#E

type ljs
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Unknown function type in file startpot' in potfit.stderr
    assert 'please define "ljs" in functions.c' in potfit.stderr

def test_apot_eam_no_potential_function(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 3
#T EAM
#I 0 0 0
#E
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Error while searching for analytic potentials' in potfit.stderr

def test_apot_eam_invalid_keyword(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 3
#T EAM
#I 0 0 0
#E

type lj
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'No cutoff found for the 1. potential' in potfit.stderr

def test_apot_eam_not_enough_params(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 3
#T EAM
#I 0 0 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Error reading analytic potentials' in potfit.stderr
    assert 'not enough parameters for potential lj' in potfit.stderr

def test_apot_eam_not_enough_param_values(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 3
#T EAM
#I 0 0 0
#E

type lj
cutoff 6.0
epsilon 0.1 0
sigma 2.5 1 4
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Could not read parameter #1 of potential #1 in file startpot' in potfit.stderr

def test_apot_eam_not_enough_potentials(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 0 3
#T EAM
#I 0 0 0
#E

type lj
cutoff 6.0
epsilon 0.1 0 1
sigma 2.5 1 4
''')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Error while searching for analytic potentials' in potfit.stderr
