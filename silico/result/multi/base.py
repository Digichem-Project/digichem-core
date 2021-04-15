# General imports.

# Silico imports.
from silico.result.result import Result_set
from silico.result.metadata import Metadata
from silico.result.atoms import Atom_list
from silico.result.ground_state import Ground_state
from silico.result.excited_states import Excited_state_list
from silico.result.emission import Relaxed_excited_state


class Merged(Result_set):
    """
    A type of result set that represents one or more separate calculations merged into one result set.
    """
    
    def __init__(self, metadatas, *args, vertical_emission = None, adiabatic_emission = None, **kwargs):
        """
        Constructor for Merged result sets.
        :param vertical_emission: An optional dictionary of Relaxed_excited_state objects representing vertical emission energies (one for each multiplicity).
        :param adiabatic_emission: An optional dictionary of Relaxed_excited_state objects representing vertical adiabatic energies (one for each multiplicity).
        """
        super().__init__(*args, **kwargs)
        self.metadatas = metadatas
        self.vertical_emission = vertical_emission if vertical_emission is not None else {}
        self.adiabatic_emission = adiabatic_emission if vertical_emission is not None else {}
            
    @classmethod
    def from_results(self, *results, alignment_class):
        """
        Create a Merged result set object from a number of result sets.
        
        :param *results: List of result sets to merge.
        :param alignment_class: An alignment class to use. 
        """
        # First, get a merged metadata object.
        metadatas = [result.metadata for result in results]
        merged_metadata = Metadata.merge(*metadatas)

        # 'Merge' our atoms.
        atoms = Atom_list.merge(*[result.atoms for result in results], charge = merged_metadata.charge)
        # And alignment.
        alignment = alignment_class(atoms, charge = merged_metadata.charge)
        
        # Merge remaining attributes.
        attrs = {}
        for attr in ["CC_energies", "MP_energies", "SCF_energies", "dipole_moment", "molecular_orbitals", "beta_orbitals", "excited_states", "vibrations", "spin_orbit_coupling"]:
            attrs[attr] = type(getattr(results[0], attr)).merge(*[getattr(result, attr) for result in results])
        
        # Get a new ground state.
        ground_state = Ground_state.from_energies(merged_metadata.charge, merged_metadata.multiplicity, attrs['CC_energies'], attrs['MP_energies'], attrs['SCF_energies'])
        
        # Build new list of energy states.
        energy_states = Excited_state_list()
        energy_states.append(ground_state)
        energy_states.extend(attrs['excited_states'])
        
        merged_results =  self(
            metadata = merged_metadata,
            metadatas = metadatas,
            ground_state = ground_state,
            atoms = atoms,
            alignment = alignment,
            energy_states = energy_states,
            **attrs
            )
        
        # Try and guess emission.
        merged_results.vertical_emission, merged_results.adiabatic_emission = Relaxed_excited_state.guess_from_results(*results)
        
        # Done.
        return merged_results
