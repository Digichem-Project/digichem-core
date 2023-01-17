# General imports.
from pathlib import Path
import cclib.io
import itertools
import pwd
import os

# Silico imports
import silico.log
from silico.result.orbital import Molecular_orbital_list,\
    Beta_orbital
from silico.result.metadata import Metadata
from silico.result.result import Result_set
from silico.result.atom import Atom_list
from silico.result.tdm import Transition_dipole_moment
from silico.result.excited_state import Excited_state_list
from silico.result.energy import Energies
from silico.result.ground_state import Ground_state
from silico.result.soc import SOC_list
from silico.result.dipole_moment import Dipole_moment
from silico.result.vibration import Vibrations_list
from silico.exception.base import Silico_exception
from silico.result.emission import Relaxed_excited_state


# NOTE: This is a repeat of the list in util to avoid circular import nonsense.
custom_parsing_formats = [
    "sir",
    "sij",
    "yaml",
    "json"
]

class Parser_abc():
    """ABC for all parsers."""
    
    def __init__(self, *log_files, raw_data = None, **kwargs):
        """
        Top level constructor for calculation parsers.
        
        :param log_files: A list of output file to analyse/parse. The first log_file given will be used for naming purposes.
        """
        # Set our name.
        self.log_file_paths = self.sort_log_files([Path(log_file) for log_file in log_files if log_file is not None])
        
        # Panic if we have no logs.
        if len(self.log_file_paths) == 0:
            raise Silico_exception("Cannot parse calculation output; no available log files. Are you sure the given path is a log file or directory containing log files?")
        
        # An object that we will populate with raw results.
        self.data = raw_data
        
        # A result set object that we'll populate with results.
        self.results = None
        
        # Parse (if we haven't already).
        try:
            if self.data is None:
                self.parse()
        except Exception:
            raise Silico_exception("Error parsing calculation result '{}'".format(self.description))

    @classmethod
    def from_logs(self, *given_log_files, **kwargs):
        """
        Intelligent constructor that will attempt to guess the location of files from a given log file(s).
        
        :param given_log_files: Output file(s) to parse or a directory of output files to parse.
        """
        found_log_files = [found_log.resolve() for given_log_file in given_log_files for found_log in self.find_log_files(given_log_file)]
        
        return self(*found_log_files, **kwargs)
    
    @classmethod
    def find_log_files(self, hint):
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
            log_types = itertools.chain(["*." + custom_format for custom_format in custom_parsing_formats], ["*.log", "*.out"])
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
            log_files = list(reversed(list(dict.fromkeys(log_files))))
            return log_files
                
        # If we have no log files, and there's a directory called Output or Result that we can use, try again using that as the hint.
        elif hint.is_dir():
            if Path(hint, "Results").is_dir():
                log_files = self.find_log_files(Path(hint, "Results"))
            
            # If we still have nothing, try the output dir.
            if len(log_files) == 0 and Path(hint, "Output").is_dir():
                log_files = self.find_log_files(Path(hint, "Output"))
                
        return log_files
    
    @classmethod
    def sort_log_files(self, log_files):
        """
        Sort a list of log files into a particular order, if required for this parser.
        """
        return log_files
    
    @property
    def log_file_path(self):
        """
        The main log file.
        """
        for log_file in self.log_file_paths:
            if log_file.suffix.lower() == ".log":
                return log_file
            
        return self.log_file_paths[0]
            
    @property
    def name(self):
        """
        Short name to describe this calculation result.
        """
        return self.log_file_path.with_suffix("").name
    
    @property
    def description(self):
        """
        A name/path that describes the file(s) being parsed, used for error messages etc.
        """
        return self.log_file_path
    
    def parse(self):
        """
        Extract results from our output files.
        """
        raise NotImplementedError("Implement in subclass")
            
    def process_all(self, alignment_class):
        """
        Get all the Result set objects produced by this parser.
        
        :param: alignment_class: An alignment class object to use to reorientate atoms. If not specified the Minimal alignment method will be used by default.
        :return: A list of the populated result sets.
        """
        self.process(alignment_class)
        return [self.results]

    def process(self, alignment_class):
        """
        Get a Result set object from this parser.
        
        :param: alignment_class: An alignment class object to use to reorientate atoms. If not specified the Minimal alignment method will be used by default.
        :return: The populated result set.
        """
        
        # Get our result set.
        self.results = Result_set(self.load_result_part(Metadata))
        
        # First get our list of MOs (because we need them for excited states too.)
        self.results.orbitals = self.load_result_part(Molecular_orbital_list)
        self.results.beta_orbitals = self.load_result_part(Molecular_orbital_list, cls = Beta_orbital)
        
        # Assign total levels.
        self.results.orbitals.assign_total_level(self.results.beta_orbitals)
        
        # Our alignment orientation data.
        self.results.atoms = self.load_result_part(alignment_class)
        self.results.raw_atoms = self.load_result_part(Atom_list)
        
        # TEDM and TMDM.
        self.results.transition_dipole_moments = self.load_result_part(Transition_dipole_moment)
        
        # Excited states.
        self.results.excited_states = self.load_result_part(Excited_state_list)
        
        # Energies.
        self.results.energies = self.load_result_part(Energies)
        
        # Our ground state.
        self.results.ground_state = self.load_result_part(Ground_state)
        
        # And a similar list but also including the ground.
        self.results.energy_states = Excited_state_list()
        self.results.energy_states.append(self.results.ground_state)
        self.results.energy_states.extend(self.results.excited_states)
        
        # SOC.
        self.results.soc = self.load_result_part(SOC_list)
        
        # PDM
        self.results.pdm = self.load_result_part(Dipole_moment)
        
        # Frequencies.
        self.results.vibrations = self.load_result_part(Vibrations_list)
        
        # Finally, try and set emission.
        self.results.emission.vertical, self.results.emission.adiabatic = Relaxed_excited_state.guess_from_results(self.results)
        
        # Return the populated result set for convenience.
        return self.results
            
    def load_result_part(self, result_cls, *, data = None, **kwargs):
        """
        Get part of a result file.
        
        For most parsers, this will simply call from_parser() of the given class, but some parsers do something more interesting.
        Any arguments other than cls will be parsed to the underlying function.
        """
        # TODO: This mechanism is no-longer needed and should be removed.
        data = self if data is None else data
        return result_cls.from_parser(data, **kwargs)


class Cclib_parser(Parser_abc):
    """
    ABC for parsers that use cclib to do most of their work for them.
    """
    
    # A dictionary of recognised auxiliary file types.
    INPUT_FILE_TYPES = {}
    
    def __init__(self, *log_files, **aux_files):
        """
        Top level constructor for calculation parsers.
        
        :param log_files: A list of output file to analyse/parse. The first log_file given will be used for naming purposes.
        :param aux_files: A dictionary of auxiliary input files related to the calculation.
        """
        # Also save our aux files, stripping None.
        self.aux_files = {name: aux_file for name,aux_file in aux_files.items() if aux_file is not None}
        
        super().__init__(*log_files)
        
    @classmethod
    def from_logs(self, *given_log_files, **kwargs):
        """
        Intelligent constructor that will attempt to guess the location of files from a given log file(s).
        
        :param given_log_files: Output file(s) to parse or a directory of output files to parse.
        """
        found_log_files = [found_log.resolve() for given_log_file in given_log_files for found_log in self.find_log_files(given_log_file)]
        
#         # Make sure we only have unique log files.
#         # We also now reverse our ordering, so that files given earlier by the user have priority.
#         log_files = list(reversed(list(dict.fromkeys(found_log_files))))
        
        # Next, have a look for aux. files.
        aux_files = {}
        
        for found_log_file in found_log_files:
            aux_files.update(self.find_aux_files(found_log_file))
            
        # Finally, update our aux_files with kwargs, so any user specified aux files take precedence.
        aux_files.update(kwargs)
        
        return self(*found_log_files, **aux_files)
    
    @classmethod
    def find_aux_files(self, hint):
        """
        Find auxiliary files from a given hint.
        
        :param hint: A path to a file to use as a hint to find additional files.
        :returns: A dictionary of found aux files.
        """
        hint = Path(hint)
        
        # Now have a look for aux. input files, which are defined by each parser's INPUT_FILE_TYPES
        aux_files = {}
        for file_type in self.INPUT_FILE_TYPES:
            for extension in file_type.extensions:
                if hint.with_suffix(extension).exists():
                    aux_files[self.INPUT_FILE_TYPES[file_type]] = hint.with_suffix(extension)
        
        return aux_files
    
    def parse(self):
        """
        Extract results from our output files.
        
        This default implementation will parse the log file with cclib and save the result to self.data.
        """
        # We start by using cclib to get most of the data we need.
        
        # Output a message (because this is slow).
        silico.log.get_logger().info("Parsing calculation result '{}'".format(self.description))
        
        # Use cclib to open our log files.
        # ccread will accept a list of log files to read, but will sometimes choke if the list contains only one entry,
        # in which case we give it only one file name.
        # ccread will also choke if we give it pathlib Paths.
        file_paths = [str(log_file) for log_file in self.log_file_paths]
        
        # Get data from cclib.
        self.data = cclib.io.ccread(file_paths if len(file_paths) > 1 else file_paths[0])
        
        # Do some setup.
        self.pre_parse()
        
        # Do our own parsing (if any).
        for log_file_path in self.log_file_paths:
            with open(log_file_path, "rt") as log_file:
                for line in log_file:
                    self.parse_output_line(log_file, line)
        
        # Make any final adjustments.
        self.post_parse()
        
    def parse_output_line(self, log_file, line):
        """
        Perform custom line-by-line parsing of an output file.
        
        This method will be called for each line of each log-file given to the parser (although be aware that some implementations may skip some lines during parsing),
        and it allows for data not supported by cclib to be extracted. It is program specific.
        """
        # Do nothing.
        
    def pre_parse(self):
        """
        Perform any setup before line-by-line parsing.
        """
        # Do nothing.
        
    def post_parse(self):
        """
        Perform any required operations after line-by-line parsing.
        """
        # Add our name.
        if self.name is not None:
            self.data.metadata['name'] = self.name
            
        # Add current username.
        # TODO: It would probably be better if we used the name of the user who owns the output file, rather than the current user...
        self.data.metadata['user'] = self.get_current_username()
        
        # Set our file paths.
        self.data.metadata['log_files'] = self.log_file_paths
        self.data.metadata['aux_files'] = self.aux_files
        
    @classmethod
    def get_current_username(self):
        """
        Get the username of the currently logged in user.
        
        :returns: The username as a string, or None if the username could not be determined.
        """
        # Based on https://stackoverflow.com/questions/842059/is-there-a-portable-way-to-get-the-current-username-in-python
        try:
            return pwd.getpwuid(os.getuid()).pw_name
        
        except Exception as e:
            return None
        
    @classmethod
    def au_to_debye(self, au):
        """
        Convert a dipole moment in au to debye.
        """
        return au * 2.541746473
    
    @classmethod
    def bohr_to_angstrom(self, bohr_distance):
        """
        Convert a length in bohr to angstrom.
        """
        return bohr_distance * 0.529177210903
    