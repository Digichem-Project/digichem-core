# Code for processing results from calculations and generating reports.

# General imports.
import itertools

# Silico imports.
from silico.exception import Result_unavailable_error
from silico.result.alignment.AAA import Adjusted_average_angle
from silico.result.alignment.AA import Average_angle
from silico.result.alignment.FAP import Furthest_atom_pair
from silico.result.alignment import Minimal
from silico.result import Result_object
import silico.result.excited_state

        
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
            transition_dipole_moments = None,
            molecular_orbitals = None,
            beta_orbitals = None,
            ground_state = None,
            excited_states = None,
            energy_states = None,
            vibrations = None,
            spin_orbit_coupling = None):
        """
        Constructor for Result_set objects.
        
        :param metadata: Optional Metadata result object.
        :param CC_energies: Optional Energy_list object of coupled-cluster energies.
        :param MP_energies: Optional Energy_list object of Moller-Plesset energies.
        :param SCF_energies: Optional Energy_list object of self-consistent field energies (SCF is the type of energy printed for normal HF and DFT).
        :param atoms: Optional Atom_list object of atom positions.
        :param dipole_moment: Optional dipole_moment object.
        :param transition_dipole_moments: Optional list of TDMs.
        :param molecular_orbitals: Optional Molecular_orbital_list object.
        :param beta_orbitals: Optional Beta MOs. If this is not None, then molecular_orbitals is assumed to refer to the Alpha MOs.
        :param excited_states: Optional Excited_state_list object.
        :param vibrations: Optional molecular Vibrations object.
        :param spin_orbit_coupling: A list of spin_orbit_coupling.
        """
        super().__init__()
        self.metadata = metadata
        self.results = (self,)
        self.CC_energies = CC_energies
        self.MP_energies = MP_energies
        self.SCF_energies = SCF_energies
        self.dipole_moment = dipole_moment
        self.transition_dipole_moments = transition_dipole_moments
        self.atoms = atoms
        self.alignment = alignment
        self.molecular_orbitals = molecular_orbitals
        self.beta_orbitals = beta_orbitals
        self.ground_state = ground_state
        self.excited_states = excited_states
        self.energy_states = energy_states
        self.vibrations = vibrations
        self.vertical_emission = {}
        self.adiabatic_emission = {}
        self.spin_orbit_coupling = spin_orbit_coupling
        
    @property
    def metadatas(self):
        """
        Property providing access to the list of metadatas of the calculations that were merged together.
        """
        return [result.metadata for result in self.results]
    
    @property
    def level_of_theory(self):
        """
        A short-hand summary of the methods and basis sets used.
        """
        methods = []
        basis_sets = []
        for metadata in self.metadatas:
            if len(metadata.converted_methods) > 0:
                method = metadata.converted_methods[-1]
#             for method in metadata.converted_methods:
                if method not in methods:
                    methods.append(method)
                
            if metadata.basis_set is not None and metadata.basis_set not in basis_sets:
                basis_sets.append(metadata.basis_set)
                        
        return "/".join(itertools.chain(methods, basis_sets))
    
    @property
    def title(self):
        """
        A string Title describing this result.
        """
        title = ", ".join(self.metadata.calculations)
        if "Excited States" in self.metadata.calculations:
            # Add multiplicity based on ES.
            mult = self.excited_states.group()
            mult_strings = [silico.result.excited_state.Energy_state.multiplicity_number_to_string(multiplicity).capitalize() for multiplicity in mult]
            if len(mult_strings) <= 2:
                # Add the strings.
                title += " ({})".format(", ".join(mult_strings))
            else:
                # Add something non-specific.
                title += " (Various Multiplicities)"
        else:
            title += " ({})".format(silico.result.excited_state.Energy_state.multiplicity_number_to_string(self.metadata.multiplicity).capitalize())
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
    def electric_transition_dipole_moment(self):
        """
        The S1 electric dipole moment, commonly referred to as THE electric transition dipole moment (although this name is ambiguous).
        
        None is returned if the S1 dipole moment is not available.
        """
        tdm = self.transition_dipole_moment
        if tdm is not None:
            return tdm.electric
        
    @property
    def magnetic_transition_dipole_moment(self):
        """
        The S1 magnetic dipole moment, commonly referred to as THE magnetic transition dipole moment (although this name is ambiguous).
        
        None is returned if the S1 dipole moment is not available.
        """
        tdm = self.transition_dipole_moment
        if tdm is not None:
            return tdm.magnetic

    @property
    def transition_dipole_moment(self):
        """
        The S1 dipole moment (both electric and magnetic), commonly referred to as THE transition dipole moment (although this name is ambiguous).
        
        None is returned if the S1 dipole moment is not available.
        """
        try:
            S1 = self.excited_states.get_state("S(1)")
            return S1.transition_dipole_moment
        except Result_unavailable_error:
            # No S1 available.
            return None
        