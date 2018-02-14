import pytest

def test_invalid_param_file_empty(potfit_apot):
    filename = 'param_file'
    potfit_apot.create_file(filename).close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert 'Missing parameter' in potfit_apot.stderr
    assert 'ntypes' in potfit_apot.stderr

def test_invalid_param_file_ntypes_empty(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert 'Missing value in parameter file' in potfit_apot.stderr
    assert 'ntypes is <undefined>' in potfit_apot.stderr

def test_invalid_param_file_ntypes_invalid(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes a')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert 'Illegal value' in potfit_apot.stderr
    assert 'ntypes is not an integer!' in potfit_apot.stderr

def test_invalid_param_file_ntypes_out_of_bounds(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes 0')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert 'Illegal value' in potfit_apot.stderr
    assert 'out of bounds' in potfit_apot.stderr

def test_invalid_param_file_ntypes_negative(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes -1')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert 'Illegal value' in potfit_apot.stderr
    assert 'out of bounds' in potfit_apot.stderr

def test_invalid_param_file_startpot_empty(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert 'Missing parameter' in potfit_apot.stderr
    assert 'startpot' in potfit_apot.stderr

def test_invalid_param_file_startpot_valid(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert 'Missing parameter' in potfit_apot.stderr
    assert '[WARNING] endpot' in potfit_apot.stderr
    assert 'config' in potfit_apot.stderr

def test_invalid_param_file_startpot_endpot_valid(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert 'Missing parameter' in potfit_apot.stderr
    assert not '[WARNING] endpot' in potfit_apot.stderr
    assert 'config' in potfit_apot.stderr

def test_invalid_param_file_config_empty(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert 'Missing parameter' in potfit_apot.stderr
    assert 'config' in potfit_apot.stderr

def test_invalid_param_file_config_valid(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert 'Missing parameter' in potfit_apot.stderr
    assert 'tempfile' in potfit_apot.stderr

def test_invalid_param_file_tempfile_empty(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert 'Missing parameter' in potfit_apot.stderr
    assert 'tempfile' in potfit_apot.stderr

def test_invalid_param_file_tempfile_valid(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert '[ERROR] Could not open file startpot' in potfit_apot.stderr

def test_invalid_param_file_eng_weight_invalid(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.write('eng_weight -10')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert '[ERROR] Illegal value' in potfit_apot.stderr
    assert 'eng_weight' in potfit_apot.stderr

def test_invalid_param_file_imdpot_empty(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.write('imdpot')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert '[ERROR] Missing value' in potfit_apot.stderr
    assert 'imdpot' in potfit_apot.stderr
    assert 'imdpotsteps' in potfit_apot.stderr

def test_invalid_param_file_imdpotsteps_invalid(potfit_apot):
    filename = 'param_file'
    f = potfit_apot.create_file(filename)
    f.write('ntypes 1\n')
    f.write('startpot startpot\n')
    f.write('endpot endpot\n')
    f.write('config config\n')
    f.write('tempfile tempfile\n')
    f.write('imdpotsteps -10')
    f.close()
    potfit_apot.run(filename)
    assert potfit_apot.has_error()
    assert '[ERROR] Illegal value' in potfit_apot.stderr
    assert 'imdpotsteps is out of bounds!' in potfit_apot.stderr
