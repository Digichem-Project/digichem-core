"""Utilities and handy functions for reading and parsing result files."""

# General imports.
import cclib.parser.gaussianparser
import cclib.parser.turbomoleparser
import multiprocessing
from functools import partial
from itertools import filterfalse
from pathlib import Path
import itertools
import shutil
import tempfile
import collections
import warnings

# Silico imports.
from silico.parse.base import Cclib_parser
from silico.parse.gaussian import Gaussian_parser
from silico.parse.turbomole import Turbomole_parser
from silico.parse.orca import Orca_parser
from silico.misc.base import is_iter
from silico.result.multi import Merged
from silico.result.alignment import Minimal
from silico.result.result import Result_set
from silico.exception.base import Silico_exception
import silico.log
from silico.parse.dump import Yaml_multi_parser, Json_multi_parser


custom_parsing_formats = {
    "sir": Yaml_multi_parser,
    "yaml": Yaml_multi_parser,
    "sij": Json_multi_parser,
    "json": Json_multi_parser,
}

def find_log_files_from_hint(hint):
    """
    Find output (log) files from a given hint.
    
    :param hint: A path to a file to use as a hint to find additional log files. hint can optionally be a directory, in which case files inside this directory will be found.
    :returns: A list of found log files.
    """
    hint = Path(hint)
    
    # First, find our parent dir.
    # hint may actually be a dir.
    if hint.is_dir():
        # Look for all .log files.
        # File extensions that we recognise.
        log_types = itertools.chain(["*." + custom_format for custom_format in custom_parsing_formats], ["*.log", ]) # "*.out" disabled for now...
        parent = hint
        log_files = [found_log_file for found_log_file in itertools.chain(*[parent.glob(log_type) for log_type in log_types])]
        
        #log_files = [Path(found_log_file) for found_log_file in iglob(str(Path(parent, "*.log")))]
        # Remove any 'silico.log' files as we know these are not calc log files.
        # We don't actually write 'silico.log' files anymore either (we use silico.out instead),
        # but older versions did...
        log_files = [log_file for log_file in log_files if log_file.name not in ["silico.log", "silico.out"]]
    else:
        parent = hint.parent
        log_files = [hint]
    
    # If we have a computational style log file, look for others.
    if hint.suffix not in ["." + custom_format for custom_format in custom_parsing_formats]:
        # Try and find job files.
        # These files have names like 'job.0', 'job.1' etc, ending in 'job.last'.
        for number in itertools.count():
            # Get the theoretical file name.
            job_file_path = Path(parent, "job.{}".format(number))
            
            # See if it exists (and isn't the log_file given to us).
            if job_file_path.exists():
                # Add to list.
                log_files.append(job_file_path)
            else:
                # We've found all the numbered files.
                break
                    
        # Look for other files.
        for maybe_file_name in ("basis", "control", "mos", "alpha", "beta", "coord", "gradient", "aoforce", "job.last", "numforce/aoforce.out"):
            maybe_file_path = Path(parent, maybe_file_name)
            
            if maybe_file_path.exists():
                # Found it.
                log_files.append(maybe_file_path)
            
    # Make sure we only have unique log files.
    # We also now reverse our ordering, so that files given earlier by the user have priority.
    if len(log_files) > 0:
        return log_files
            
    # If we have no log files, and there's a directory called Output or Result that we can use, try again using that as the hint.
    elif hint.is_dir():
        if Path(hint, "Results").is_dir():
            log_files = find_log_files_from_hint(Path(hint, "Results"))
        
        # If we still have nothing, try the output dir.
        if len(log_files) == 0 and Path(hint, "Output").is_dir():
            log_files = find_log_files_from_hint(Path(hint, "Output"))
            
    return log_files


def find_log_files(*hints):
    """
    Find log files from a number of given hints.
    
    Each hint should point to an existing log file to parse, or a directory containing such log files.
    Each log file given (and found) should refer to the same calculation. 
    """
    # Get a list of found log files.
    log_files = [found_log for hint in hints for found_log in find_log_files_from_hint(hint)]
    
    # Make sure we only have unique log files.
    # We also now reverse our ordering, so that files given earlier by the user have priority.
    log_files = list(reversed(list(dict.fromkeys([log_file.resolve() for log_file in log_files]))))

    return log_files
    

def class_from_log_files(*log_files, format_hint = "auto"):
    """
    Get a parser class based on some calculation log files.
    
    :param format_hint: A hint as to the format of the given log files. Either 'auto' (to guess), 'log' (calc log file), 'sir' (silico result file) or 'sid' (silico database file).
    """
    if format_hint in custom_parsing_formats:
        return custom_parsing_formats[format_hint]
    
    elif format_hint == "auto" and len(log_files) > 0 and log_files[0].suffix[1:] in custom_parsing_formats:
        return custom_parsing_formats[log_files[0].suffix[1:]]
    
    if format_hint not in ["cclib", "auto"]:
        raise ValueError("Unrecognised format hint '{}'".format(format_hint))
    
    # We'll use cclib to guess the file type for us.
    try:
        log_file_type = type(cclib.io.ccopen([str(found_log_file) for found_log_file in log_files]))
        
    except Exception as e:
        # cclib couldn't figure out the file type, it probably wasn't a .log file.
        raise Silico_exception("Could not determine file type of file(s): '{}'; are you sure these are computational log files?".format(", ".join((str(log_file) for log_file in log_files)))) from e
    
    # Either get a more complex parser if we have one, or just return the base parser.
    if log_file_type == cclib.parser.gaussianparser.Gaussian:
        return Gaussian_parser
    elif log_file_type == cclib.parser.turbomoleparser.Turbomole:
        return Turbomole_parser
    elif log_file_type == cclib.parser.orcaparser.ORCA:
        return Orca_parser
    else:
        return Cclib_parser

def from_log_files(*log_files, format_hint = "auto", **aux_files):
    """
    Get an output file parser of appropriate type.
    
    This is a convenience function.
    
    :param format_hint: A hint as to the format of the given log files. Either 'auto' (to guess), 'log' (calc log file), 'sir' (silico result file) or 'sid' (silico database file).
    """
    log_files = find_log_files(*log_files)
    return class_from_log_files(*log_files, format_hint = format_hint).from_logs(*log_files, **aux_files)
    
def parse_calculation(*log_files, alignment_class = None, parse_all = False, format_hint = "auto", **aux_files):
    """
    Parse a single calculation result.
    
    :param log_files: A number of calculation result files corresponding to the same calculation.
    :param alignment_class: An optional alignment class to use to reorientate the molecule.
    :param parse_all: Whether to parse all results in the given log file. If True, a list of result sets will be returned, if False, only the first result will be returned if there are multiple.
    :param format_hint: A hint as to the format of the given log files. Either 'auto' (to guess), 'log' (calc log file), 'sir' (silico result file) or 'sid' (silico database file).
    :param aux_files: Optional auxiliary calculation files corresponding to the calculation.
    :return: A single Result_set object.
    """
    if alignment_class is None:
        alignment_class = Minimal
        
    # Open files for reading (handles archives for us).
    with open_for_parsing(*log_files) as open_log_files:
        
        if parse_all:
            return from_log_files(*open_log_files, format_hint = format_hint, **aux_files).process_all(alignment_class)
        
        else:       
            return from_log_files(*open_log_files, format_hint = format_hint, **aux_files).process(alignment_class)

def multi_parser(log_file, alignment_class, format_hint = "auto"):
        """
        The inner function which will be called in parallel to parse files.
        
        :param format_hint: A hint as to the format of the given log files. Either 'auto' (to guess), 'log' (calc log file), 'sir' (silico result file) or 'sid' (silico database file).
        """
        try:
            return parse_calculation(log_file, alignment_class = alignment_class, parse_all = True, format_hint = format_hint)
            
        except Exception:
            silico.log.get_logger().warning("Unable to parse calculation result file '{}'; skipping".format(log_file), exc_info = True)

def parse_multiple_calculations(*log_files, alignment_class = None, pool = None, init_func = None, init_args = None, format_hint = "auto", processes = 1):
    """
    Parse a number of separate calculation results in parallel.
    
    If the argument 'pool' is not given, a multiprocessing.pool object will first be created using the arguments init_func and num_cpu.
    
    :param log_files: A number of calculation result files corresponding to different calculations.
    :param alignment_class: An optional alignment class to use to reorientate the molecule.
    :param pool: An optional subprocessing.pool object to use for parallel parsing.
    :param init_func: An optional function to call to init each newly created process.
    :param format_hint: A hint as to the format of the given log files. Either 'auto' (to guess), 'log' (calc log file), 'sir' (silico result file) or 'sid' (silico database file).
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
        result_lists = list(
            filterfalse(lambda x: x is None,
                pool.map(partial(multi_parser, alignment_class = alignment_class, format_hint = format_hint), log_files)
            )
        )
        
        return [result for result_list in result_lists for result in result_list]
    
    finally:
        # Do some cleanup if we need to.
        if own_pool:
            pool.close()
    
def parse_calculations(*results, alignment_class = None, format_hint = "auto", aux_files = None):
    """
    Get a single result object by parsing a number of computational log files.
    
    Multiple log files can be given both from the same calculation, or from multiple different calculations.
    If multiple different calculations are given, the individually parsed results will be merged together (which may give bizarre results if the calculations are unrelated, eg if they are of different molecules).
    
    Example:
        parse_calculations(['calc1/primary.log', 'calc1/secondary.log'], 'calc2/calc.log', 'calc3/calc.log')
    Would parse three separate calculations (calc1, calc2 and calc3), of which the first is contained in two output files (primary.log and secondary.log), merging the result sets together.
    
    :param log_files: A list of paths to computational chemistry log files to parse. If more than one file is given, each is assumed to correspond to a separate calculation in which case the parsed results will be merged together. In addition, each given 'log file' can be an iterable of log file paths, which will be considered to correspond to an individual calculation.
    :param alignment_class: An alignment class to use to reorientate each molecule.
    :param format_hint: A hint as to the format of the given log files. Either 'auto' (to guess), 'log' (calc log file), 'sir' (silico result file) or 'sid' (silico database file).
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
            result_list = log_result_or_list
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
            result_list = parse_calculation(*logs, alignment_class = alignment_class, parse_all = True, format_hint = format_hint, **aux)
        
        # Add to our list.
        parsed_results.extend(result_list)
    
    # If we have more than one result, merge them together.
    if len(parsed_results) > 1:
        return Merged.from_results(*parsed_results, alignment_class = alignment_class)
    else:
        return parsed_results[0]


class open_for_parsing():
    """
    A context manager for opening log files for parsing. Returns a list of pathlib.Path objects suitable for parsing.
    
    Currently, the main purpose of this context manager is to intelligently handle unpacking of archives (.zip, .tar etc) for parsing.
    """
    
    def __init__(self, *log_files):
        log_files = [Path(log_file).resolve() for log_file in log_files]
        
        # Check we haven't been given any duplicate log files.
        # This is just for convenience, if duplicates have been given the user has probably made a mistake.               
        duplicates = [path for path, number in collections.Counter(log_files).items() if number > 1]
        
        for duplicate in duplicates:
            warnings.warn("Ignoring duplicate log file: ".format(duplicate))
        
        # Remove duplicates but retain order.
        self.log_files = list(dict.fromkeys(log_files).keys())
        
        # A list of tempfile.TemporaryDirectory objects that should be closed when we are finished.
        self.temp_dirs = []
    
    @property
    def archive_formats(self):
        """
        Get a list of supported archive formats.
        """
        return list(itertools.chain(*[extensions for name, extensions, desc in shutil.get_unpack_formats()]))
        
    def __enter__(self):
        """
        'Open' files for reading.
        """
        new_log_files = []
        
        formats = self.archive_formats
        
        for log_file in self.log_files:
            
            found_child_archive = None
            
            # If 'log_file' is a directory, check for an archive inside called 'Output.xxx'.
            for archive_format in formats:
                child_archive = Path(log_file, "Output" + archive_format)
                if child_archive.exists():
                    if not found_child_archive:
                        # Found an Output dir archive, use this instead.
                        new_log_files.extend(self.extract(child_archive))
                        found_child_archive = child_archive
                    
                    else:
                        # For now, only care about the first.
                        warnings.warn("Ignoring subsequent Output archive '{}'; already found '{}'".format(child_archive, found_child_archive))
            
            # No need to check 'found_child_archive' here; a file cannot simultaneously be a directory containing an archive and also an archive itself.
            if "".join(log_file.suffixes) in formats:
                # This is an archive format.                
                # Add any files/directories that were unpacked.
                new_log_files.extend(self.extract(log_file))
                
            elif not found_child_archive:
                # Non-archive file, add normally.
                new_log_files.append(log_file)
            
        return new_log_files
    
    def extract(self, file_name):
        """
        Extract an archive and return the contained log files.
        """
        # Get a temp dir to extact to.
        tempdir = tempfile.TemporaryDirectory()
        self.temp_dirs.append(tempdir)
        
        # Extract to it.
        silico.log.get_logger().info("Extracting archive '{}'...".format(file_name))
        shutil.unpack_archive(file_name, tempdir.name)
        
        # Add any files/directories that were unpacked.
        return Path(tempdir.name).glob("*")
    
    def cleanup(self):
        """
        Perform cleanup of any open files.
        
        Alias for __exit__().
        """
        return self.__exit__()
    
    def __exit__(self, exc_type = None, exc_value = None, traceback = None):
        """
        'Close' any open files.
        """
        for tempdir in self.temp_dirs:
            tempdir.cleanup()
        
        