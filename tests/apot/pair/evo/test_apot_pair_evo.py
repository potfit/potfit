import pytest

def test_apot_pair_evo_threshold_empty(potfit):
    potfit.create_param_file(evo_threshold='')
    potfit.call_makeapot('startpot', '-n 1 -i pair -f eopp_sc')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Missing value in parameter file' in potfit.stderr
    assert 'evo_threshold is <undefined>' in potfit.stderr

def test_apot_pair_evo_threshold_invalid(potfit):
    potfit.create_param_file(evo_threshold='foo')
    potfit.call_makeapot('startpot', '-n 1 -i pair -f eopp_sc')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Illegal value in parameter file' in potfit.stderr
    assert 'evo_threshold is not a double' in potfit.stderr
    assert 'value = foo' in potfit.stderr

def test_apot_pair_evo_threshold_out_of_bounds(potfit):
    potfit.create_param_file(evo_threshold=-1)
    potfit.call_makeapot('startpot', '-n 1 -i pair -f eopp_sc')
    potfit.create_config_file()
    potfit.run()
    assert potfit.has_error()
    assert 'Illegal value in parameter file' in potfit.stderr
    assert 'evo_threshold is out of bounds' in potfit.stderr
