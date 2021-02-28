# General imports.
from logging import getLogger

# Silico imports.
import silico
from silico.result.result import Result_set
from silico.result.metadata import Metadata
from silico.result.atoms import Atom_list
from silico.result.ground_state import Ground_state
from silico.result.excited_states import Excited_state_list
from silico.result.emission import Relaxed_excited_state
from silico.exception.base import Silico_exception


class Merged(Result_set):
    """
    A type of result set that represents one or more separate calculations merged into one result set.
    """
    
    def __init__(self, metadatas, *args, **kwargs):
        """
        Constructor for Merged result sets.
        """
        super().__init__(*args, **kwargs)
        self.metadatas = metadatas
    
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
        
    @classmethod
    def from_results(self, *results, alignment_class):
        """
        Create a Merged result set object from a number of result sets.
        
        :param *results: List of result sets to merge.
        :param alignment_class: An alignment class to use. 
        """
        # First, get a merged metadata object.
        metadatas = [result.metadata for result in results]
        merged_metadata = Metadata.merge(metadatas)
        
        # Get a new ground state.
        ground_state = Ground_state.from_energies(merged_metadata.system_charge, merged_metadata.system_multiplicity, merged_metadata.CC_energies, merged_metadata.MP_energies, merged_metadata.SCF_energies)

        # 'Merge' our atoms.
        atoms = Atom_list.merge([result.atoms for result in results], charge = merged_metadata.system_charge)
        # And alignment.
        alignment = alignment_class(atoms, charge = merged_metadata.system_charge)
        
        # Merge remaining attributes.
        attrs = {}
        for attr in ["CC_energies", "MP_energies", "SCF_energies", "dipole_moment", "molecular_orbitals", "beta_orbitals", "excited_states", "vibrations", "spin_orbit_coupling"]:
            attrs[attr] = type(getattr(results[0], attr)).merge([getattr(result, attr) for result in results])
        
        # Build new list of energy states.
        energy_states = Excited_state_list()
        energy_states.append(ground_state)
        energy_states.extend(attrs['excited_states'])
        
        return self(
            metadata = merged_metadata,
            metadatas = metadatas,
            ground_state = ground_state,
            atoms = atoms,
            alignment = alignment,
            energy_states = energy_states,
            **attrs
            )
