import pytest
import numpy

# Set numpy errors (not sure why this isn't the default...)
numpy.seterr(invalid = 'raise', divide = 'raise')

def pytest_addoption(parser):
    parser.addoption(
        "--run-slow", action="store_true", default=False, help="run slow tests")
    parser.addoption(
        "--run-parallel-only", action="store_true", default=False, help="only run tests that can be ran in parallel")
    parser.addoption(
        "--run-series-only", action="store_true", default=False, help="only run tests that can be ran in series"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as being slow")


def pytest_collection_modifyitems(config, items):
#     if config.getoption("--run-slow"):
#         # Asked to run slow tests, do nothing.
#         return
    

    skip_submit = pytest.mark.skip(reason="need --run-slow option to run")
    for item in items:
        if     ("slow" in item.keywords and not config.getoption("--run-slow") and not config.getoption("--run-series-only")) \
            or ("slow" in item.keywords and config.getoption("--run-parallel-only")) \
            or ("slow" not in item.keywords and config.getoption("--run-series-only")):
            item.add_marker(skip_submit)
