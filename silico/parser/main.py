from silico.result.molecular_orbitals import Molecular_orbital_list,\
    Beta_orbital
from silico.result.result import Metadata, Result_set
from silico.result.atoms import Atom_list
from silico.result.transition_dipole_moment import Transition_dipole_moment
from silico.result.excited_states import Excited_state_list
from silico.result.energy import CC_energy_list, MP_energy_list, SCF_energy_list
from silico.result.ground_state import Ground_state
from silico.result.spin_orbit_coupling import Spin_orbit_coupling
from silico.result.dipole_moment import Dipole_moment
from silico.result.vibrations import Vibration_list
from silico.exception.base import Silico_exception
from logging import getLogger
import silico

class Parser(Result_set):
    """
    Top-level abstract class for calculation result parsers.
    """
    
    def __init__(self, name):
        """
        Top level constructor for calculation parsers.
        
        :param name: The name of the calculation.
        """
        # Set our name.
        self.name = name
        
        # An object that we will populate with raw results.
        self.data = None
        
        # A result set object that we'll populate with results.
        self.results = None
        
        # Parse.
        # Output a message (because this is slow).
        getLogger(silico.logger_name).info("Parsing calculation result '{}'".format(self.description))
        try:
            self.parse()
        except Exception:
            raise Silico_exception("Error parsing calculation result '{}'".format(self.description))
        
    @property
    def description(self):
        """
        A name/path that describes the file(s) being parsed, used for error messages etc.
        """
        return self.name
    
    def parse(self):
        """
        Extract results from our output files.
        """
        raise NotImplementedError()
        
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
        self.results.CC_energies = CC_energy_list.from_parser(self),
        self.results.MP_energies = MP_energy_list.from_parser(self),
        self.results.SCF_energies = SCF_energy_list.from_parser(self),
        
        # Our ground state.
        self.results.ground_state = Ground_state.from_parser(self)
        
        # And a similar list but also including the ground.
        self.results.energy_states = Excited_state_list()
        self.results.energy_states.append(self.ground_state)
        self.results.energy_states.extend(self.excited_states)
        
        # SOC.
        self.results.spin_orbit_coupling = Spin_orbit_coupling.list_from_parser(self)
        
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