import pytest

@pytest.fixture(scope="module")
def potfit():
    import sys
    sys.path.insert(0, str(pytest.config.rootdir))
    import potfit
    p = potfit.Potfit(__file__, 'apot', 'pair')
    yield p
    p.clear()

def test_invalid_param_file_empty(potfit):
    filename = 'param_file'
    potfit.create_file(filename).close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Missing parameter' in potfit.stderr
    assert 'ntypes' in potfit.stderr

def test_invalid_param_file_ntypes_empty(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Illegal value' in potfit.stderr
    assert 'out of bounds' in potfit.stderr

def test_invalid_param_file_ntypes_invalid(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes a')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Illegal value' in potfit.stderr
    assert 'out of bounds' in potfit.stderr

def test_invalid_param_file_ntypes_out_of_bounds(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes 0')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Illegal value' in potfit.stderr
    assert 'out of bounds' in potfit.stderr

def test_invalid_param_file_ntypes_negative(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes -1')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Illegal value' in potfit.stderr
    assert 'out of bounds' in potfit.stderr

def test_invalid_param_file_startpot_empty(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Missing parameter' in potfit.stderr
    assert 'startpot' in potfit.stderr

def test_invalid_param_file_startpot_valid(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Missing parameter' in potfit.stderr
    assert '[WARNING] endpot' in potfit.stderr
    assert 'config' in potfit.stderr

def test_invalid_param_file_startpot_endpot_valid(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Missing parameter' in potfit.stderr
    assert not '[WARNING] endpot' in potfit.stderr
    assert 'config' in potfit.stderr

def test_invalid_param_file_config_empty(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Missing parameter' in potfit.stderr
    assert 'config' in potfit.stderr

def test_invalid_param_file_config_valid(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Missing parameter' in potfit.stderr
    assert 'tempfile' in potfit.stderr

def test_invalid_param_file_tempfile_empty(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert 'Missing parameter' in potfit.stderr
    assert 'tempfile' in potfit.stderr

def test_invalid_param_file_tempfile_valid(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert '[ERROR] Could not open file startpot' in potfit.stderr

def test_invalid_param_file_eng_weight_invalid(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.write('eng_weight -10')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert '[ERROR] Illegal value' in potfit.stderr
    assert 'eng_weight' in potfit.stderr

def test_invalid_param_file_imdpot_empty(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.write('imdpot')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert '[ERROR] Missing value' in potfit.stderr
    assert 'imdpot' in potfit.stderr
    assert 'imdpotsteps' in potfit.stderr

def test_invalid_param_file_imdpotsteps_invalid(potfit):
    filename = 'param_file'
    f = potfit.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.write('imdpotsteps -10')
    f.close()
    potfit.run(filename)
    assert potfit.returncode == 1
    assert potfit.has_error_msg()
    assert '[ERROR] Illegal value' in potfit.stderr
    assert 'imdpotsteps is out of bounds!' in potfit.stderr
