import warnings
import logging
import sys
from multiprocessing import Lock

from silico.logging.format import Variable_formatter

# The handler object silico uses for logging.
LOGGING_HANDLER = None

# The name of the silico logger, can be passed to logging.get_logger() to get the silico logging object.
LOGGER_NAME = "silico"

# This probably isn't necessary.
LOGGER_LOCK = Lock()


def init_logger():
    """
    Init the program wide logger.
    """
    global LOGGING_HANDLER
    
    logging.captureWarnings(True)
    
    logger = logging.getLogger(LOGGER_NAME)
    warnings_logger = logging.getLogger("py.warnings")
    
    # The console handler, where we'll print most messages.
    #LOGGING_HANDLER = logging.StreamHandler(sys.stderr)
    LOGGING_HANDLER = Handler(sys.stderr)
    # Handle everything.
    LOGGING_HANDLER.setLevel(logging.DEBUG)
    # Set its formatter.
    var_formatter = Variable_formatter(logger, default_warning_formatter = warnings.formatwarning)
    LOGGING_HANDLER.setFormatter(var_formatter)
    # Add the handler.
    logger.addHandler(LOGGING_HANDLER)
    warnings_logger.addHandler(LOGGING_HANDLER)
    warnings.formatwarning = var_formatter.formatWarning


def get_logger():
    """
    Get the logger used by all parts of silico.
    """
    return logging.getLogger(LOGGER_NAME)


class Handler(logging.StreamHandler):
    """
    """
    
    def emit(self, *args, **kwargs):
        """
        """
        # Acquire lock.
        LOGGER_LOCK.acquire()
        try:
            # Do the actual work.
            return super().emit(*args, **kwargs)
            
        finally:
            # Release lock.
            LOGGER_LOCK.release()