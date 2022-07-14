import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--run-slow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as being slow")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-slow"):
        # Asked to run slow tests, do nothing.
        return
    
    skip_submit = pytest.mark.skip(reason="need --run-slow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_submit)
