import pytest

@pytest.fixture(scope="module")
def potfit():
    import sys
    sys.path.insert(0, str(pytest.config.rootdir))
    import potfit
    p = potfit.Potfit(__file__, 'apot', 'pair')
    yield p
    p.clear()

def test_cmd_line_args_none(potfit):
    potfit.run_with_args()
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Usage' in potfit.stderr
    assert potfit.binary_name in potfit.stderr
    assert '<paramfile>' in potfit.stderr

def test_cmd_line_args_not_exists(potfit):
    filename = 'not_exists'
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Could not open parameter file' in potfit.stderr
    assert filename in potfit.stderr

def test_cmd_line_args_no_permission(potfit):
    filename = 'no_permission'
    potfit.create_file(filename, permission=0).close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Could not open parameter file' in potfit.stderr
    assert filename in potfit.stderr
