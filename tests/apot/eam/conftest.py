import pytest

def pytest_runtest_logstart(nodeid, location):
    path = location[0]
    if not path.startswith('apot/eam'):
        raise pytest.UsageError("Please run the tests from the tests/ base directory!")

potfit_obj = None

def get_potfit_obj(request):
    import sys
    sys.path.insert(0, str(request.config.rootdir))
    import potfit
    global potfit_obj
    if potfit_obj == None:
        potfit_obj = potfit.Potfit(__file__, interaction='eam', model='apot')
    return potfit_obj

@pytest.fixture()
def potfit(request):
    p = get_potfit_obj(request)
    p.reset()
    yield p
    p.clear()
