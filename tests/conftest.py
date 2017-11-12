import os
import pytest

def pytest_sessionstart(session):
    """ Remove all lingering asan files """
    clean_dir(os.path.curdir)

def clean_dir(dirname):
    for root, dir, files in os.walk(dirname):
        for item in files:
            if item.startswith('asan'):
                os.remove(os.path.join(root, item))
        for item in dir:
            clean_dir(item)
