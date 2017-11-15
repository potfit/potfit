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
