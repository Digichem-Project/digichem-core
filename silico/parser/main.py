from silico.parser.base import Parser
from silico.parser.gaussian import Gaussian_parser
from silico.parser.turbomole import Turbomole_parser

def get_parser(log_file):
    """
    Get an output file parser of appropriate type.
    
    This is a convenience function.
    """  
    # First get a general parser.
    base_parser = Parser(log_file)
    
    # Either get a more complex parser if we have one, or just return the base parser.
    if base_parser.data.metadata['package'] == "Gaussian":
        return Gaussian_parser.from_parser(base_parser)
    elif base_parser.data.metadata['package'] == "Turbomole":
        return Turbomole_parser.from_parser(base_parser)