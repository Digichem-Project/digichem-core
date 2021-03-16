# Contains classes for representing the ground state of our system.
# This file is closely related to the excited_states.py file.
# These definitions may change or move.

from silico.result.excited_states import Energy_state
from silico.exception.base import Result_unavailable_error
from builtins import classmethod

class Ground_state(Energy_state):
    """
    Class for representing the ground state of a system/molecule.
    """
    
    def __init__(self, charge, multiplicity, energy):
        """
        Constructor for ground state objects.
        
        :param charge: The charge of the ground states as a signed integer.
        :param multiplicity: The multiplicity as an integer of float.
        :param energy: The absolute energy of the GS in eV.
        """
        # The 'level' refers to the ordering of energy states, the ground state is always the lowest so have level 0.
        # Similarly, the multiplicity_level refers to the ordering of energy states that share the same multiplicity. Regardless of what multiplicity our GS has, we will be the lowest by definition.
        # Energy is our energy above the GS, so that will always be 0.
        super().__init__(level = 0, multiplicity = multiplicity, multiplicity_level = 0, energy = 0.0)
        self.charge = charge
        # Also save our absolute energy.
        self.absolute_energy = energy
        
    def __eq__(self, other):
        """
        Is this ground state object equal to another?
        """
        return self.charge == other.charge and self.multiplicity == other.multiplicity and self.energy == other.energy
    
    @classmethod
    def from_parser(self, parser):
        """
        Create a Ground_state object from an output file parser.
        
        As different energy types can be present in a single calculation, the following (arbitrary?) precedence is used:
        1) Coupled-cluster energy.
        2) The highest level Moller-Plesset energy.
        3) SCF energy.
        
        :param parser: An output file parser.
        """
        return self.from_energies(parser.results.metadata.charge, parser.results.metadata.multiplicity, parser.results.CC_energies, parser.results.MP_energies, parser.results.SCF_energies)
    
    @classmethod
    def from_energies(self, charge, multiplicity, CC_energies, MP_energies, SCF_energies):
        """
        Create a ground state object selecting from various lists of energies.
        
        At least one of CC_energies, MP_energies or SCF_energies should not be empty.
        
        :param charge: Charge of the ground state.
        :param multiplicity: Multiplicity of the ground state.
        :param CC_energies: List of coupled-cluster energies.
        :param MP_energies: List of Moller-Plesset energies.
        :param SCF_energies: List of self-consistent field energies.
        """
        if len(CC_energies) > 0:
            energy = CC_energies.final
        elif len(MP_energies) > 0:
            energy = MP_energies.final
        elif len(SCF_energies) > 0:
            energy = SCF_energies.final
        else:
            # We have no energies that we can understand. Give up.
            raise Result_unavailable_error("Ground_state", "No ground state energies that we understand are available")
        # Return our constructed object.
        return self(charge, multiplicity, energy)
    
    
    
    