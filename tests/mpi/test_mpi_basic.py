import pytest

def test_mpi_apot_pair_single_process(potfit):
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
    potfit.create_config_file()
    potfit.run(mpi=1)
    assert potfit.has_error()
    assert 'Running potfit with a single MPI process is not supported!' in potfit.stderr
    assert 'This creates a hugh overhead and slows down the process significantly!' in potfit.stderr

def test_mpi_apot_pair_too_many_processes(potfit):
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
    potfit.create_config_file()
    potfit.run(mpi=4)
    assert potfit.has_error()
    assert 'Running potfit with more MPI processes than configurations is not supported!' in potfit.stderr
    assert 'You provided 1 configurations and requested 4 processes' in potfit.stderr
