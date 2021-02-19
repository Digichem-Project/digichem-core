from silico.result.result import Result_set


class Merged(Result_set):
    """
    A type of result set that represents one or more separate calculations merged into one result set.
    """
    
    def __init__(self, *results):
        """
        """
        
        super().__init__(
            metadata,
            CC_energies,
            MP_energies,
            SCF_energies,
            atoms,
            alignment,
            dipole_moment,
            molecular_orbitals,
            beta_orbitals,
            ground_state,
            excited_states,
            energy_states,
            vertical_emission,
            adiabatic_emission,
            vibrations,
            spin_orbit_coupling):