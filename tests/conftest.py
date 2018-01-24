import os
import pytest

def pytest_sessionstart(session):
    """ Remove all lingering asan files before running """
    _clean_dir(os.path.curdir)

def _clean_dir(dirname):
    for root, dir, files in os.walk(dirname):
        for item in files:
            if item.startswith('asan'):
                os.remove(os.path.join(root, item))
        for item in dir:
            _clean_dir(item)

@pytest.fixture(scope='session', autouse=True)
def configure_html_report_env(request):
    keep_items = ['BUILD_', 'GIT_', 'JENKINS', 'JOB_NAME', 'NODE_NAME']
    for key in list(request.config._metadata.keys()):
      remove = True
      for allowed_string in keep_items:
        if allowed_string in key:
          remove = False
          break
      if remove:
        del request.config._metadata[key]
