import warnings
import logging
import sys
import textwrap
from subprocess import CalledProcessError

# The handler object silico uses for logging.
LOGGING_HANDLER = None

# The name of the silico logger, can be passed to logging.get_logger() to get the silico logging object.
LOGGER_NAME = "silico"


def init_logger():
    """
    Init the package wide logger.
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
    
    
def set_logging_level(log_level, verbose = None):
    """
    Set the logging level of the silico logger object.
    
    :param log_level: The base logging level as a string.
    :param verbose: An integer which specifies how much to increase the logging level by. If verbose is zero (or None) then the logging level is determined only by log_level.
    """
    logger = get_logger()
    
    # Set from log_level first.
    if log_level == "OFF":
        logger.setLevel(60)
    else:
        logger.setLevel(log_level)
    
    # Now adjust with verbosity.
    if verbose is not None:
        # Set from verbosity.
        new_level = logger.level - verbose * 10
        
        # Don't allow us to reach 0 (because this is actually 'UNSET').
        if new_level <= 0:
            new_level = 10
        
        # And set.
        logger.setLevel(new_level)
        
#     # Adjust warnings if necessary.
#     if logger.level == 10:
#         warnings.simplefilter('always', DeprecationWarning)


def get_logger():
    """
    Get the logger used by all parts of silico.
    """
    return logging.getLogger(LOGGER_NAME)


class Handler(logging.StreamHandler):
    """
    """
    
class Variable_formatter(logging.Formatter):
    """
    The logging formatter used by silico, the format changes depending on a logger's log_level.
    """
    
    default_format = 'silico: %(levelname)s: %(message)s'
    
    def __init__(self, logger, fmt=None, datefmt=None, style='%', *, default_warning_formatter):
        if fmt is None:
            fmt = self.default_format
        super().__init__(fmt, datefmt, style)
        # Save our logger.
        self.logger = logger
        self.default_warning_formatter = default_warning_formatter
        
    def formatWarning(self, message, category, filename, lineno, line=None):
        """
        Format a warning (from the warnings module).
        
        :return: The formatted warning message.
        """
        if self.logger.getEffectiveLevel() > logging.DEBUG:
            return "{}: {}".format(category.__name__, message)
        else:
            return self.default_warning_formatter(message, category, filename, lineno, line=None)
                
    def formatException(self, exc_info):
        """
        Format an exception.
        
        In Variable_formatter, how we do this depends on our log level.
        
        :return: The formatted exception.
        """
        if self.logger.getEffectiveLevel() > logging.DEBUG:
            return self.exception_to_str(exc_info[1])
        else:
            return textwrap.indent(super().formatException(exc_info), "  ")
        
    @classmethod
    def process_output_to_str(self, process_exception):
        """
        Get additional descriptive text from a CalledProcessError exception.
        
        
        """
        output = ""
        if process_exception.stdout is not None and process_exception.stdout != "\n":
            output += "\n" + process_exception.stdout.strip()
        if process_exception.stderr is not None and process_exception.stderr != "\n":
            output += "\n" + process_exception.stderr.strip()
        return output
        
    @classmethod
    def exception_to_str(self, exception):
        """
        Recursively get a string representation of an exception.
        
        :return: A string representation of exception and any prior exceptions.
        """
        excstr = textwrap.indent("{}: {}".format(type(exception).__name__, exception), "\t")
        
        # If the exception is a CalledProcessError, we'll also append the stdout/stderr.
        if isinstance(exception, CalledProcessError):
            excstr += textwrap.indent(self.process_output_to_str(exception), "\t\t")
        
        if exception.__cause__ is not None:
            excstr += "\n{}".format(self.exception_to_str(exception.__cause__))
        elif exception.__context__ is not None and not exception.__suppress_context__:
            excstr += "\n{}".format(self.exception_to_str(exception.__context__))
        
        return excstr