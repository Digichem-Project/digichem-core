from logging import getLogger
from pathlib import Path
import cclib.io

from silico.result.molecular_orbitals import Molecular_orbital_list,\
    Beta_orbital
from silico.result.result import Metadata, Result_set
from silico.result.atoms import Atom_list
from silico.result.transition_dipole_moment import Transition_dipole_moment
from silico.result.excited_states import Excited_state_list
from silico.result.energy import CC_energy_list, MP_energy_list, SCF_energy_list
from silico.result.ground_state import Ground_state
from silico.result.spin_orbit_coupling import SOC_list
from silico.result.dipole_moment import Dipole_moment
from silico.result.vibrations import Vibration_list
from silico.exception.base import Silico_exception
import silico

class Parsed_data():
    """
    An empty class used to store attributes.
    """

class Parser(Result_set):
    """
    Top-level abstract class for calculation result parsers.
    """
    
    def __init__(self, *log_files, _parsed_data = None):
        """
        Top level constructor for calculation parsers.
        
        :param log_files: A list of output file to analyse/parse. The first log_file given will be used for naming purposes.
        """
        # Set our name.
        self.log_file_paths = [Path(log_file) for log_file in log_files if log_file is not None]
        
        # An object that we will populate with raw results.
        self.data = _parsed_data
        
        # A result set object that we'll populate with results.
        self.results = None
        
        # Parse.
        try:
            self.parse()
        except Exception:
            raise Silico_exception("Error parsing calculation result '{}'".format(self.description))
    
    @classmethod
    def from_parser(self, parser):
        """
        Get a new parser object from an existing parser object.
        
        :param parser: An existing parser object.
        :return: A new parser, which will be of the current type (which is normally a different type to parser).
        """
        return self.from_log(*parser.log_file_paths, _parsed_data = parser.data)
    
    @classmethod
    def from_log(self, log_file, **kwargs):
        """
        Intelligent constructor that will attempt to guess the location of files from a given log file.
        
        :param log_file: Output file to parse.
        """
        return self(log_file, **kwargs)
    
    @property
    def name(self):
        """
        Short name to describe this calculation result.
        """
        return self.log_file_paths[0].with_suffix("").name
    
    @property
    def description(self):
        """
        A name/path that describes the file(s) being parsed, used for error messages etc.
        """
        return self.log_file_paths[0]
    
    def parse(self):
        """
        Extract results from our output files.
        
        This default implementation will parse the log file with cclib and save the result to self.data, but only if self.data is already None.
        """
        # We start by using cclib to get most of the data we need.
        
        # Read and parse the file with cclib if we don't already have some data.
        if self.data is None:
            # Output a message (because this is slow).
            getLogger(silico.logger_name).info("Parsing calculation result '{}'".format(self.description))
            
            # Use cclib to open our log files.
            # ccread will accept a list of log files to read, but will sometimes choke if the list contains only one entry,
            # in which case we give it only one file name.
            # ccread will also choke if we give it pathlib Paths.
            file_paths = [str(log_file) for log_file in self.log_file_paths]
            
            self.data = cclib.io.ccread(file_paths if len(file_paths) > 1 else file_paths[0])
        
    def process(self, alignment_class):
        """
        Get a Result set object from this parser.
        
        :param: alignment_class: An alignment class object to use to reorientate atoms.
        :return: The populated result set.
        """
        # Get our result set.
        self.results = Result_set()
        
        # First get our list of MOs (because we need them for excited states too.
        self.results.molecular_orbitals = Molecular_orbital_list.from_parser(self)
        self.results.beta_orbitals = Molecular_orbital_list.from_parser(self, cls = Beta_orbital)
        
        # Metadata.
        self.results.metadata = Metadata.from_parser(self)
        
        # Our alignment orientation data.
        self.results.alignment = alignment_class.from_parser(self)
        self.results.atoms = Atom_list.from_parser(self)
        
        # TDM.
        self.results.transition_dipole_moments = Transition_dipole_moment.list_from_parser(self)
        
        # Excited states.
        self.results.excited_states = Excited_state_list.from_parser(self)
        
        # Energies.
        self.results.CC_energies = CC_energy_list.from_parser(self)
        self.results.MP_energies = MP_energy_list.from_parser(self)
        self.results.SCF_energies = SCF_energy_list.from_parser(self)
        
        # Our ground state.
        self.results.ground_state = Ground_state.from_parser(self)
        
        # And a similar list but also including the ground.
        self.results.energy_states = Excited_state_list()
        self.results.energy_states.append(self.results.ground_state)
        self.results.energy_states.extend(self.results.excited_states)
        
        # SOC.
        #self.results.spin_orbit_coupling = Spin_orbit_coupling.list_from_parser(self)
        self.results.spin_orbit_coupling = SOC_list.from_parser(self)
        
        # PDM
        self.results.dipole_moment = Dipole_moment.from_parser(self)
        
        # Finally, frequencies.
        self.results.vibrations = Vibration_list.from_parser(self)
        
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
    
                
class Sub_parser():
    """
    Top-level class for sub-parsers.
    """
    
    def __init__(self, parser):
        """
        :param parser: The top-level parser.
        """
        self.parser = parser