import os
import pytest
from subprocess import run

def pytest_runtest_logstart(nodeid, location):
    if not location[0].startswith('kim'):
        raise pytest.UsageError("Please run the tests from the tests/ base directory!")

potfit_obj = None

def has_kim_support():
    pkg_config_cmd = os.environ.get('PKG_CONFIG', None) or 'pkg-config'
    return run('{0} --exists libkim-api-v2'.format(pkg_config_cmd).split()).returncode == 0

def supports_models():
    pkg_config_cmd = os.environ.get('PKG_CONFIG', None) or 'pkg-config'
    libexecdir = run([pkg_config_cmd, 'libkim-api-v2', '--variable=libexecdir'], capture_output=True).stdout.decode().strip()
    models = run([os.path.join(libexecdir, 'kim-api-v2', 'kim-api-v2-collections-info'), 'models'], capture_output=True).stdout.decode().strip().split('\n')
    models = [line.split()[1] for line in models]
    return True

def get_potfit_obj(request):
    import sys
    sys.path.insert(0, str(request.config.rootdir))
    import potfit
    global potfit_obj
    if potfit_obj == None:
        potfit_obj = potfit.Potfit(__file__, model='kim')
    return potfit_obj

@pytest.fixture()
def potfit(request):
    if not has_kim_support():
        pytest.skip('No KIM installation found!')
    if not supports_models():
        pytest.skip('Not all required KIM models found!')
    p = get_potfit_obj(request)
    p.reset()
    yield p
    if os.path.isfile('kim/kim.log'):
        os.remove('kim/kim.log')
    p.clear()
