# Code for processing results from calculations and generating reports.

# General imports.
from logging import getLogger

# Silico imports.
from silico.exception import Result_unavailable_error
from silico.result.alignment.AAA import Adjusted_average_angle
from silico.result.alignment.AA import Average_angle
from silico.result.alignment.FAP import Furthest_atom_pair
from silico.result.alignment import Minimal
from silico.result import Result_object
from silico.exception.base import Silico_exception
import silico.result.excited_states
from silico.result.emission import Relaxed_excited_state

        
class Result_set(Result_object):
    """
    Class that represents a collection of results from a calculation.
    
    This class is a bit heavy and might not be around for long...
    """
    
    def __init__(
            self,
            metadata = None,
            CC_energies = None,
            MP_energies = None,
            SCF_energies = None,
            atoms = None,
            alignment = None,
            dipole_moment = None,
            molecular_orbitals = None,
            beta_orbitals = None,
            ground_state = None,
            excited_states = None,
            energy_states = None,
            vertical_emission = None,
            adiabatic_emission = None,
            vibrations = None,
            spin_orbit_coupling = None):
        """
        Constructor for Result_set objects.
        
        :param gaussian_log_file: The Gaussian log file from which these results were read.
        :param metadata: Optional Metadata result object.
        :param CC_energies: Optional Energy_list object of coupled-cluster energies.
        :param MP_energies: Optional Energy_list object of Moller-Plesset energies.
        :param SCF_energies: Optional Energy_list object of self-consistent field energies (SCF is the type of energy printed for normal HF and DFT).
        :param atoms: Optional Atom_list object of atom positions.
        :param dipole_moment: Optional dipole_moment object.
        :param molecular_orbitals: Optional Molecular_orbital_list object.
        :param beta_orbitals: Optional Beta MOs. If this is not None, then molecular_orbitals is assumed to refer to the Alpha MOs.
        :param excited_states: Optional Excited_state_list object.
        :param vertical_emission: A Relaxed_excited_state object representing the vertical emission energy.
        :param adiabatic_emission: A Relaxed_excited_state object representing the adiabatic emission energy.
        :param vibrations: Optional molecular Vibrations object.
        :param spin_orbit_coupling: A list of spin_orbit_coupling.
        """
        super().__init__()
        self.metadata = metadata
        self.CC_energies = CC_energies
        self.MP_energies = MP_energies
        self.SCF_energies = SCF_energies
        self.dipole_moment = dipole_moment
        self.atoms = atoms
        self.alignment = alignment
        self.molecular_orbitals = molecular_orbitals
        self.beta_orbitals = beta_orbitals
        self.ground_state = ground_state
        self.excited_states = excited_states
        self.energy_states = energy_states
        self.vibrations = vibrations
        self.vertical_emission = vertical_emission
        self.adiabatic_emission = adiabatic_emission
        self.spin_orbit_coupling = spin_orbit_coupling
    
    
    
    @property
    def title(self):
        """
        A string Title describing this result.
        """
        title = ", ".join(self.metadata.calculations)
        if "Excited States" in self.metadata.calculations:
            # Add multiplicity based on ES.
            mult = self.excited_states.group()
            mult_strings = [silico.result.excited_states.Energy_state.multiplicity_number_to_string(multiplicity).capitalize() for multiplicity in mult]
            if len(mult_strings) <= 2:
                # Add the strings.
                title += " ({})".format(", ".join(mult_strings))
            else:
                # Add something non-specific.
                title += " (Various Multiplicities)"
        else:
            title += " ({})".format(silico.result.excited_states.Energy_state.multiplicity_number_to_string(self.metadata.system_multiplicity).capitalize())
        return title
    
    @property
    def energy(self):
        """
        The total energy of this calculation.
        
        This convenience property is the energy at the highest level of theory available (CC > MP > SCF).
        
        :raises Result_unavailable_error: If no total energy is available.
        """
        # Try CC first.
        if len(self.CC_energies) > 0:
            return self.CC_energies.final
        elif len(self.MP_energies) > 0:
            return self.MP_energies.final
        else:
            return self.SCF_energies.final
    
        
    @property
    def transition_dipole_moment(self):
        """
        The S1 dipole moment, commonly referred to as THE transition dipole moment (although this name is ambiguous).
        
        None is returned if the S1 dipole moment is not available.        
        """
        try:
            S1 = self.excited_states.get_state("S(1)")
            return S1.transition_dipole_moment
        except Result_unavailable_error:
            # No S1 available.
            return None

        
    def add_emission(self,
            vertical_emission_ground_result = None,
            adiabatic_emission_ground_result = None,
            emission_excited_result = None,
            emission_excited_state = None
            ):
        """
        Add additional result sets containing emission energies to this result set.
        
        # TODO: Rethink the way we handle emission energy.
        
        :param vertical_emission_ground_result: An optional additional Result_set object which will be used as the ground state for calculation of vertical emission energy.
        :param adiabatic_emission_ground_result: An optional additional Result_set object which will be used as the ground state for calculation of adiabatic emision energy.
        :param emission_excited_result: An optional additional Result_set object which will be used as the excited state for calculation of vertical and/or adiabatic emission energy. If emission_ground_result is give, then emission_excited_result is not optional.
        :param emission_excited_state: Optionally either an Excited_state object or a string describing one ('S(1)' etc) for use in calculating vertical and/or adiabatic emission energy.
        """
        # Set our emission energies.
        # For vertical, we can also use the adiabatic excited energy
        try:
            self._set_emission('vertical', vertical_emission_ground_result, emission_excited_result, emission_excited_state)
        except Exception:
            getLogger(silico.logger_name).warning("Could not load vertical emission energy", exc_info = True)
        try:
            self._set_emission('adiabatic', adiabatic_emission_ground_result, emission_excited_result, emission_excited_state)
        except Exception:
            getLogger(silico.logger_name).warning("Could not load adiabatic emission energy", exc_info = True)
        
    def _set_emission(self, transition_type, emission_ground_result = None, emission_excited_result = None, emission_excited_state = None):
        """
        Helper function.
        """
        if transition_type != "vertical" and transition_type != "adiabatic":
            raise ValueError("Unknown transition_type '{}'".format(transition_type))
        
        # Set our emission energy.
        if emission_excited_result is not None:
            setattr(self, transition_type + "_emission", Relaxed_excited_state.from_results(
                self,
                ground_state_result = emission_ground_result,
                excited_state_result = emission_excited_result,
                transition_type = transition_type,
                excited_state = emission_excited_state))
        elif emission_ground_result is not None:
            raise Silico_exception("Cannot calculate emission energy; no excited state given for ground state '{}'".format(emission_ground_result.metadata.name))
        
    
                