# General imports.
from pathlib import Path
import pwd
import os

from digichem.exception.base import Digichem_exception
from digichem.result.orbital import Molecular_orbital_list,\
    Beta_orbital
from digichem.result.metadata import Metadata
from digichem.result.result import Result_set
from digichem.result.atom import Atom_list
from digichem.result.tdm import Transition_dipole_moment
from digichem.result.excited_state import Excited_state_list
from digichem.result.energy import Energies
from digichem.result.ground_state import Ground_state
from digichem.result.soc import SOC_list
from digichem.result.dipole_moment import Dipole_moment
from digichem.result.vibration import Vibrations_list
from digichem.result.emission import Relaxed_excited_state
from digichem.result.nmr import NMR_shielding, NMR_spin_couplings_list, NMR_list
from digichem.result.alignment.base import Minimal, Alignment


# NOTE: This is a repeat of the list in util to avoid circular import nonsense.
custom_parsing_formats = [
    "sir",
    "sij",
    "yaml",
    "json"
]

class Parser_abc():
    """ABC for all parsers."""

    def __init__(self, *, raw_data = None, options, **kwargs):
        """
        Top level constructor for calculation parsers.
        """
        # An object that we will populate with raw results.
        self.data = raw_data
        
        # A result set object that we'll populate with results.
        self.results = None

        # Config options.
        self.options = options
        
        # Parse (if we haven't already).
        try:
            if self.data is None:
                self.parse()
        except Exception:
            raise Digichem_exception("Error parsing calculation result '{}'".format(self.description))
    
    @property
    def name(self):
        """
        Short name to describe this calculation result.
        """
        return "Parser"
    
    @property
    def description(self):
        """
        A name/path that describes the file(s) being parsed, used for error messages etc.
        """
        return "Parser"

    def parse(self):
        """
        Extract results from our output files.
        """
        self._parse()
    
        # Make any final adjustments.
        self.post_parse()
    
    def _parse(self):
        """
        Extract results from our output files.
        """
        raise NotImplementedError("Implement in subclass")
        
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
    
    def post_parse(self):
        """
        Perform any required operations after line-by-line parsing.
        """ 
        # Add current username.
        # TODO: It would probably be better if we used the name of the user who owns the output file, rather than the current user...
        self.data.metadata['user'] = self.get_current_username()
            
    def process_all(self):
        """
        Get all the Result set objects produced by this parser.
        
        :param options: A Digichem options nested dictionary containing options to control parsing.
        :return: A list of the populated result sets.
        """
        self.process()
        return [self.results]
    
    def process(self):
        """
        Get a Result set object from this parser.
        
        :param options: A Digichem options nested dictionary containing options to control parsing.
        :return: The populated result set.
        """
        # Get our result set.
        self.results = Result_set(
            _id = self.data._id,
            metadata = Metadata.from_parser(self),
            aux = self.data._aux if hasattr(self.data, '_aux') else None
            )
        
        alignment_class = Alignment.from_class_handle(self.options['alignment']) if self.options['alignment'] is not None else Minimal
        
        # First get our list of MOs (because we need them for excited states too.)
        self.results.orbitals = Molecular_orbital_list.from_parser(self)
        self.results.beta_orbitals = Molecular_orbital_list.from_parser(self, cls = Beta_orbital)
        
        # Assign total levels.
        self.results.orbitals.assign_total_level(self.results.beta_orbitals)
        
        # Our alignment orientation data.
        self.results.atoms = alignment_class.from_parser(self)
        self.results.raw_atoms = Atom_list.from_parser(self)
        
        # TEDM and TMDM.
        self.results.transition_dipole_moments = Transition_dipole_moment.from_parser(self)
        
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
        self.results.soc = SOC_list.from_parser(self)
        
        # PDM
        self.results.pdm = Dipole_moment.from_parser(self)
        
        # Frequencies.
        self.results.vibrations = Vibrations_list.from_parser(self)
        
        # NMR.
        self.results.nmr_shieldings = NMR_shielding.dict_from_parser(self)
        self.results.nmr_couplings = NMR_spin_couplings_list.from_parser(self)
        self.results.nmr = NMR_list.from_parser(self)
        
        # Finally, try and set emission.
        self.results.emission.vertical, self.results.emission.adiabatic = Relaxed_excited_state.guess_from_results(self.results)
        
        # Return the populated result set for convenience.
        return self.results


class File_parser_abc(Parser_abc):
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
            raise Digichem_exception("Cannot parse calculation output; no available log files. Are you sure the given path is a log file or directory containing log files?")
        
        super().__init__(raw_data=raw_data, **kwargs)

    @classmethod
    def from_logs(self, *log_files, **kwargs):
        """
        Intelligent constructor that will attempt to guess the location of aux files from a given log file(s).
        
        :param log_files: Output file(s) to parse or a directory of output files to parse.
        """
        # This default implementation does nothing smart.
        return self(*log_files, **kwargs)
    
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
        
    def post_parse(self):
        """
        Perform any required operations after line-by-line parsing.
        """
        super().post_parse()

        # Add our name.
        if self.name is not None:
            self.data.metadata['name'] = self.name
        
        # Set our file paths.
        self.data.metadata['log_files'] = self.log_file_paths
        self.data.metadata['auxiliary_files'] = self.auxiliary_files
