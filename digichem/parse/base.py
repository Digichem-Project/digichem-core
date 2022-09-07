# General imports.
from pathlib import Path
import cclib.io
from glob import iglob
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

# TODO: Why is this a Result_set?!
class Parser(Result_set):
    """
    Top-level abstract class for calculation result parsers.
    """
    
    # A dictionary of recognised auxiliary file types.
    INPUT_FILE_TYPES = {}
    
    def __init__(self, *log_files, **aux_files):
        """
        Top level constructor for calculation parsers.
        
        :param log_files: A list of output file to analyse/parse. The first log_file given will be used for naming purposes.
        :param aux_files: A dictionary of auxiliary input files related to the calculation.
        """
        # Set our name.
        self.log_file_paths = [Path(log_file) for log_file in log_files if log_file is not None]
        
        # Panic if we have no logs.
        if len(self.log_file_paths) == 0:
            raise Silico_exception("Cannot parse calculation output; no available log files. Are you sure the given path is a log file or directory containing log files?")
        
        # Also save our aux files, stripping None.
        self.aux_files = {name: aux_file for name,aux_file in aux_files.items() if aux_file is not None}
        
        # An object that we will populate with raw results.
        self.data = None
        
        # A result set object that we'll populate with results.
        self.results = None
        
        # Parse.
        try:
            self.parse()
        except Exception:
            raise Silico_exception("Error parsing calculation result '{}'".format(self.description))
        
    @classmethod
    def from_logs(self, *given_log_files, **kwargs):
        """
        Intelligent constructor that will attempt to guess the location of files from a given log file(s).
        
        :param log_file: Output file to parse or a directory of output files to parse.
        """
        found_log_files = [found_log.resolve() for given_log_file in given_log_files for found_log in self.find_log_files(given_log_file)]
        
        # Make sure we only have unique log files.
        # We also now reverse our ordering, so that files given earlier by the user have priority.
        log_files = list(reversed(list(dict.fromkeys(found_log_files))))
        
        # Next, have a look for aux. files.
        aux_files = {}
        
        for found_log_file in found_log_files:
            aux_files.update(self.find_aux_files(found_log_file))
            
        # Finally, update our aux_files with kwargs, so any user specified aux files take precedence.
        aux_files.update(kwargs)
        
        return self(*log_files, **aux_files)
    
    @classmethod
    def find_log_files(self, hint):
        """
        Find output (log) files from a given hint.
        
        :param hint: A path to a file to use as a hint to find additional log files. hint can optionally be a directory, in which case files inside this directory will be found.
        :returns: A list of found log files.
        """
        hint = Path(hint)
        
        attempt = 1
        while attempt < 2:
            attempt += 1
            # First, find our parent dir.
            # hint may actually be a dir.
            if hint.is_dir():
                # Look for all .log files.
                parent = hint
                log_files = [Path(found_log_file) for found_log_file in iglob(str(Path(parent, "*.log")))]
                # Remove any 'silico.log' files as we know these are not calc log files.
                # We don't actually write 'silico.log' files anymore either (we use silico.out instead),
                # but older versions did...
                log_files = [log_file for log_file in log_files if log_file.name != "silico.log"]
            else:
                parent = hint.parent
                log_files = [hint]
                                
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
            for maybe_file_name in ("basis", "control", "mos", "alpha", "beta", "coord", "gradient", "aoforce", "job.last"):
                maybe_file_path = Path(parent, maybe_file_name)
                
                if maybe_file_path.exists():
                    # Found it.
                    log_files.append(maybe_file_path)
                    
            # If we have no log files, and there's a directory called Output that we can use, try again using that as the hint.
            if len(log_files) == 0 and hint.is_dir() and Path(hint, "Output").is_dir():
                hint = Path(hint, "Output")
                attempt -= 1
                
        
        return log_files
    
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
        
        self.data = cclib.io.ccread(file_paths if len(file_paths) > 1 else file_paths[0])
        
        self.parse_metadata()
        
    def parse_metadata(self):
        """
        Parse additional calculation metadata.
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
        
    def process(self, alignment_class):
        """
        Get a Result set object from this parser.
        
        :param: alignment_class: An alignment class object to use to reorientate atoms. If not specified the Minimal alignment method will be used by default.
        :return: The populated result set.
        """
        
        # Get our result set.
        self.results = Result_set(metadata = Metadata.from_parser(self))
        
        # First get our list of MOs (because we need them for excited states too.)
        self.results.orbitals = Molecular_orbital_list.from_parser(self)
        self.results.beta_orbitals = Molecular_orbital_list.from_parser(self, cls = Beta_orbital)
        
        # Assign total levels.
        self.results.orbitals.assign_total_level(self.results.beta_orbitals)
        
        # Metadata.
#         self.results.metadata = Metadata.from_parser(self)
        
        # Our alignment orientation data.
        self.results.alignment = alignment_class.from_parser(self)
        self.results.atoms = Atom_list.from_parser(self)
        
        # TEDM and TMDM.
        #self.results.transition_dipole_moments = Electric_transition_dipole_moment.list_from_parser(self)
        #self.results.magnetic_transition_dipole_moments = Magnetic_transition_dipole_moment.list_from_parser(self)
        self.results.transition_dipole_moments = Transition_dipole_moment.list_from_parser(self)
        
        # Excited states.
        self.results.excited_states = Excited_state_list.from_parser(self)
        
        # Energies.
        self.results.energies = Energies.from_parser(self)
        
        # Our ground state.
        self.results.ground_state = Ground_state.from_parser(self)
        
        # And a similar list but also including the ground.
        self.results.energy_states = Excited_state_list()
        self.results.energy_states.append(self.results.ground_state)
        self.results.energy_states.extend(self.results.excited_states)
        
        # SOC.
        #self.results.soc = Spin_orbit_coupling.list_from_parser(self)
        self.results.soc = SOC_list.from_parser(self)
        
        # PDM
        self.results.dipole_moment = Dipole_moment.from_parser(self)
        
        # Frequencies.
        self.results.vibrations = Vibrations_list.from_parser(self)
        
        # Finally, try and set emission.
        self.results.vertical_emission, self.results.adiabatic_emission = Relaxed_excited_state.guess_from_results(self.results)
        
        # Return the populated result set for convenience.
        return self.results
        
        
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
    