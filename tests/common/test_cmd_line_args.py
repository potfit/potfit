import pytest

def test_cmd_line_args_none(potfit):
    potfit.run(param_file=None)
    assert potfit.has_error()
    assert 'Usage' in potfit.stderr
    assert potfit.binary_name in potfit.stderr
    assert '<paramfile>' in potfit.stderr

def test_cmd_line_args_not_exists(potfit):
    filename = 'not_exists'
    potfit.run(filename)
    assert potfit.has_error()
    assert 'Could not open parameter file' in potfit.stderr
    assert filename in potfit.stderr

def test_cmd_line_args_no_permission(potfit):
    filename = 'no_permission'
    potfit.create_file(filename, permission=0).close()
    potfit.run(filename)
    assert potfit.has_error()
    assert 'Could not open parameter file' in potfit.stderr
    assert filename in potfit.stderr
