import cclib.parser.gaussianparser
import cclib.parser.turbomoleparser

from silico.parser.base import Parser
from silico.parser.gaussian import Gaussian_parser
from silico.parser.turbomole import Turbomole_parser

def get_parser(log_file):
    """
    Get an output file parser of appropriate type.
    
    This is a convenience function.
    """
    # First get child files if we are a dir.
    log_files, parent, kwfiles = Parser.find_logs(log_file)
    
    # We'll use cclib to guess the file type for us.
    log_file_type = type(cclib.io.ccopen([str(found_log_file) for found_log_file in log_files]))
    
    # Either get a more complex parser if we have one, or just return the base parser.
    if log_file_type == cclib.parser.gaussianparser.Gaussian:
        return Gaussian_parser.from_log(log_file)
    elif log_file_type == cclib.parser.turbomoleparser.Turbomole:
        return Turbomole_parser.from_log(log_file)
    else:
        return Parser.from_log(log_file)
    