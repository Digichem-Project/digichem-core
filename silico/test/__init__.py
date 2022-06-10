"""Unit tests for silico"""

from pathlib import Path
import pkg_resources

def data_directory():
    return Path(pkg_resources.resource_filename('silico', 'test/data'))

import silico.logging
silico.logging.set_logging_level("DEBUG")