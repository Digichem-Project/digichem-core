from logging import getLogger
import logging
import sys
import os
from pathlib import Path

import silico.logging
import warnings

def init_logger(logger_name):
    """
    Init the program wide logger.
    """
    logging.captureWarnings(True)
    
    logger = getLogger(logger_name)
    warnings_logger = logging.getLogger("py.warnings")
    
    # The console handler, where we'll print most messages.
    consoleHandler = logging.StreamHandler(sys.stderr)
    # Handle everything.
    consoleHandler.setLevel(logging.DEBUG)
    # Set its formatter.
    var_formatter = silico.logging.Variable_formatter(logger, default_warning_formatter = warnings.formatwarning)
    consoleHandler.setFormatter(var_formatter)
    # Add the handler.
    logger.addHandler(consoleHandler)
    warnings_logger.addHandler(consoleHandler)
    warnings.formatwarning = var_formatter.formatWarning
    
    
    
def init_obabel():
    """
    Set-up openbabel.
    
    When frozen with pyinstaller we take a version of the openbabel C library with us (along with the relevant python bindings of course).
    This library is split into several .so files corresponding to the various formats obabel supports, and while the main libopenbabel.so file is found automatically, these supplementary library files are not.
    So, when we are frozen, we manually set the location of these library files so openbabel will work.
    If we are not frozen we do not do this as we expect openbabel to be correctly configured.
    """
#     # The sys attribute 'frozen' is our flag, '_MEIPASS' is the dir location.
#     # https://pyinstaller.readthedocs.io/en/stable/runtime-information.html#run-time-information
    if silico.frozen:
        # We need to tell openbabel where its library components are.
        os.environ['BABEL_LIBDIR'] = str(Path(sys._MEIPASS, "openbabel", "lib", "3.1.1"))
        
        # And also data.
        os.environ['BABEL_DATADIR'] = str(Path(sys._MEIPASS, "openbabel", "data", "3.1.1"))