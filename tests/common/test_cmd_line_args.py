import pytest

def test_cmd_line_args_none(potfit_apot):
    potfit_apot.run(param_file=None)
    assert potfit_apot.has_error()
    assert 'Usage' in potfit_apot.stderr
    assert potfit_apot.binary_name in potfit_apot.stderr
    assert '<paramfile>' in potfit_apot.stderr

def test_cmd_line_args_not_exists(potfit_apot):
    filename = 'not_exists'
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert 'Could not open parameter file' in potfit_apot.stderr
    assert filename in potfit_apot.stderr

def test_cmd_line_args_no_permission(potfit_apot):
    filename = 'no_permission'
    potfit_apot.create_file(filename, permission=0).close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert 'Could not open parameter file' in potfit_apot.stderr
    assert filename in potfit_apot.stderr
