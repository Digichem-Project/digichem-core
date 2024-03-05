import atexit
from contextlib import ExitStack

# TODO: Can switch to non-backport importlib.resources once >= 3.9
import importlib_resources

def get_resource(name):
    """
    Get a pathlib path object to a package resource.
    """
    file_manager = ExitStack()
    atexit.register(file_manager.close)
    ref = importlib_resources.files('digichem') / name
    return file_manager.enter_context(importlib_resources.as_file(ref))