# General imports.
import cclib.parser.gaussianparser
import cclib.parser.turbomoleparser
import multiprocessing.pool
from functools import partial
from itertools import filterfalse

# Silico imports.
from silico.parser.main import Parser
from silico.parser.gaussian import Gaussian_parser
from silico.parser.turbomole import Turbomole_parser
from silico.misc.base import is_iter
from silico.result.multi.base import Merged
from silico.result.alignment import Minimal
from silico.result.result import Result_set
from silico.exception.base import Silico_exception
import silico.logging

def class_from_log_files(*log_files):
    """
    Get a parser class based on some calculation log files.
    """
    # First get child files if we are a dir.
    found_log_files = Parser.find_log_files(log_files[0])
    
    # We'll use cclib to guess the file type for us.
    try:
        log_file_type = type(cclib.io.ccopen([str(found_log_file) for found_log_file in found_log_files]))
        
    except Exception as e:
        # cclib couldn't figure out the file type, it probably wasn't a .log file.
        raise Silico_exception("Could not determine file type of file(s): '{}'; are your sure these are computational log files?".format(", ".join((str(log_file) for log_file in log_files)))) from e
    
    # Either get a more complex parser if we have one, or just return the base parser.
    if log_file_type == cclib.parser.gaussianparser.Gaussian:
        return Gaussian_parser
    elif log_file_type == cclib.parser.turbomoleparser.Turbomole:
        return Turbomole_parser
    else:
        return Parser

def from_log_files(*log_files, **aux_files):
    """
    Get an output file parser of appropriate type.
    
    This is a convenience function.
    """
    return class_from_log_files(*log_files).from_logs(*log_files, **aux_files)
    
def parse_calculation(*log_files, alignment_class = None, **aux_files):
    """
    Parse a single calculation result.
    
    :param log_files: A number of calculation result files corresponding to the same calculation.
    :param alignment_class: An optional alignment class to use to reorientate the molecule.
    :param aux_files: Optional auxiliary calculation files corresponding to the calculation.
    :return: A single Result_set object.
    """
    if alignment_class is None:
        alignment_class = Minimal
            
    return from_log_files(*log_files, **aux_files).process(alignment_class)

def multi_parser(log_file, alignment_class):
        """
        The inner function which will be called in parallel to parse files.
        """
        try:
            return parse_calculation(log_file, alignment_class = alignment_class)
            
        except Exception:
            silico.logging.get_logger().warning("Unable to parse calculation result file '{}'; skipping".format(log_file), exc_info = True)

def parse_multiple_calculations(*log_files, alignment_class = None, pool = None, init_func = None, init_args = None, processes = 1):
    """
    Parse a number of separate calculation results in parallel.
    
    If the argument 'pool' is not given, a multiprocessing.pool object will first be created using the arguments init_func and num_CPUs.
    
    :param log_files: A number of calculation result files corresponding to different calculations.
    :param alignment_class: An optional alignment class to use to reorientate the molecule.
    :param pool: An optional subprocessing.pool object to use for parallel parsing.
    :param init_func: An optional function to call to init each newly created process.
    :param processes: The max number of processes to create the new pool object with; if the number of given log_files is less than processes, then len(log_files) will be used instead.
    """
    if len(log_files) == 0:
        # Give up now.
        return []
    
    if len(log_files) < processes:
        processes = len(log_files)
    
    # Sort out our pool if we need to.
    own_pool = False
    if pool is None:
        own_pool = True
        pool = multiprocessing.Pool(processes, initializer = init_func, initargs = init_args if init_args is not None else [])
        #pool = multiprocessing.pool.ThreadPool(processes)
    
    # Do some parsing.
    try:
        results = list(
            filterfalse(lambda x: x is None,
                pool.map(partial(multi_parser, alignment_class = alignment_class), log_files)
            )
        )
        
        return results
    
    finally:
        # Do some cleanup if we need to.
        if own_pool:
            pool.close()
    
def parse_calculations(*results, alignment_class = None, aux_files = None):
    """
    Get a single result object by parsing a number of computational log files.
    
    Multiple log files can be given both from the same calculation, or from multiple different calculations.
    If multiple different calculations are given, the individually parsed results will be merged together (which may give bizarre results if the calculations are unrelated, eg if they are of different molecules).
    
    Example:
        parse_calculations(['calc1/primary.log', 'calc1/secondary.log'], 'calc2/calc.log', 'calc3/calc.log')
    Would parse three separate calculations (calc1, calc2 and calc3), of which the first is contained in two output files (primary.log and secondary.log), merging the result sets together.
    
    :param log_files: A list of paths to computational chemistry log files to parse. If more than one file is given, each is assumed to correspond to a separate calculation in which case the parsed results will be merged together. In addition, each given 'log file' can be an iterable of log file paths, which will be considered to correspond to an individual calculation.
    :param alignment_class: An alignment class to use to reorientate each molecule.
    :param aux_files: A list of dictionaries of auxiliary files. The ordering of aux_files should match that of log_files.
    :return: A single Result_set object (or child thereof).
    """
    if alignment_class is None:
            alignment_class = Minimal
    
    if aux_files is None:
        aux_files = []
    
    # Process our given log files.
    parsed_results = []   
    for index, log_result_or_list in enumerate(results):
        if isinstance(log_result_or_list, Result_set):
            # Nothing we need to do.
            result = log_result_or_list
        else:
            # First check if this is an already parsed result.
            # Next, decide if this is a real log file path or is actually a list.
            if not isinstance(log_result_or_list, str) and is_iter(log_result_or_list):
                # Multiple paths.
                logs = log_result_or_list
            else:
                # Single path.
                logs = (log_result_or_list,)
                
            # See if we have some aux files we can use.
            try:
                aux = aux_files[index]
            except IndexError:
                aux = {}
                
            # Parse this calculation output.
            result = parse_calculation(*logs, alignment_class = alignment_class, **aux)
        
        # Add to our list.
        parsed_results.append(result)
    
    # If we have more than one result, merge them together.
    if len(parsed_results) > 1:
        return Merged.from_results(*parsed_results, alignment_class = alignment_class)
    else:
        return parsed_results[0]
            
        
        
        
        
    