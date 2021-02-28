# General imports.
import cclib.parser.gaussianparser
import cclib.parser.turbomoleparser

# Silico imports.
from silico.parser.main import Parser
from silico.parser.gaussian import Gaussian_parser
from silico.parser.turbomole import Turbomole_parser
from silico.misc.base import is_iter
from silico.result.multi.base import Merged

def get_parser(*log_files, **aux_files):
    """
    Get an output file parser of appropriate type.
    
    This is a convenience function.
    """
    # First get child files if we are a dir.
    found_log_files = Parser.find_log_files(log_files[0])
    
    # We'll use cclib to guess the file type for us.
    log_file_type = type(cclib.io.ccopen([str(found_log_file) for found_log_file in found_log_files]))
    
    # Either get a more complex parser if we have one, or just return the base parser.
    if log_file_type == cclib.parser.gaussianparser.Gaussian:
        return Gaussian_parser.from_logs(*log_files, **aux_files)
    elif log_file_type == cclib.parser.turbomoleparser.Turbomole:
        return Turbomole_parser.from_logs(*log_files, **aux_files)
    else:
        return Parser.from_logs(*log_files, **aux_files)
    
def parse_log_files(*log_files, alignment_class = None, **aux_files):
    """
    Get a single result object by parsing a number of computational log files.
    
    Multiple log files can be given both from the same calculation, or from multiple different calculations.
    If multiple different calculations are given, the individually parsed results will be merged together (which may give bizarre results if the calculations are unrelated, eg if they are of different molecules).
    
    Example:
        parse_log_files(['calc1/primary.log', 'calc1/secondary.log'], 'calc2/calc.log', 'calc3/calc.log')
    Would parse three separate calculations (calc1, calc2 and calc3), of which the first is contained in two output files (primary.log and secondary.log), merging the result sets together.
    
    :param log_files: A list of paths to computational chemistry log files to parse. If more than one file is given, each is assumed to correspond to a separate calculation in which case the parsed results will be merged together. In addition, each given 'log file' can be an iterable of log file paths, which will be considered to correspond to an individual calculation.
    :param alignment_class: An optional alignment class to use to reorientate each molecule.
    :return: A single result_set.
    """
    # Process our given log files.
    results = []   
    for log_file_or_list in log_files:
        # First decide if this is a real log file path or is actually a list.
        if not isinstance(log_file_or_list, str) and is_iter(log_file_or_list):
            # Multiple paths.
            logs = log_file_or_list
        else:
            # Single path.
            logs = (log_file_or_list,)
            
        # Parse this calculation output.
        result = get_parser(logs, **aux_files).process(alignment_class)
        
        # Add to our list.
        results.append(result)
    
    # If we have more than one result, merge them together.
    if len(results) > 1:
        return Merged.from_results(results, alignment_class = alignment_class)
    else:
        return results[0]
            
        
        
        
        
    