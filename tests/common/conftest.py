import pytest

def pytest_runtest_logstart(nodeid, location):
    path = location[0]
    if not path.startswith('common'):
        raise pytest.UsageError("Please run the tests from the tests/ base directory!")

potfit_apot_obj = None
potfit_tab_obj = None

def get_potfit_apot_obj():
    import sys
    sys.path.insert(0, str(pytest.config.rootdir))
    import potfit
    global potfit_apot_obj
    if potfit_apot_obj is None:
        potfit_apot_obj = potfit.Potfit(__file__, 'apot', 'pair')
    return potfit_apot_obj

def get_potfit_tab_obj():
    import sys
    sys.path.insert(0, str(pytest.config.rootdir))
    import potfit
    global potfit_tab_obj
    if potfit_tab_obj is None:
        potfit_tab_obj = potfit.Potfit(__file__, 'tab', 'pair')
    return potfit_tab_obj

@pytest.fixture()
def potfit_apot():
    p = get_potfit_apot_obj()
    p.reset()
    yield p
    p.clear()

@pytest.fixture()
def potfit_tab():
    p = get_potfit_tab_obj()
    p.reset()
    yield p
    p.clear()
