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
import warnings

        
class Result_set(Result_object):
    """
    Class that represents a collection of results from a calculation.
    
    This class is a bit heavy and might not be around for long...
    """
    
    def __init__(
            self,
            metadata = None,
            energies = None,
            atoms = None,
            alignment = None,
            pdm = None,
            transition_dipole_moments = None,
            orbitals = None,
            beta_orbitals = None,
            ground_state = None,
            excited_states = None,
            energy_states = None,
            vibrations = None,
            soc = None):
        """
        Constructor for Result_set objects.
        
        :param metadata: Optional Metadata result object.
        :param energies: Energies result object.
        :param atoms: Optional Atom_list object of atom positions.
        :param pdm: Optional dipole_moment object.
        :param transition_dipole_moments: Optional list of TDMs.
        :param orbitals: Optional Molecular_orbital_list object.
        :param beta_orbitals: Optional Beta MOs. If this is not None, then molecular_orbitals is assumed to refer to the Alpha MOs.
        :param excited_states: Optional Excited_state_list object.
        :param vibrations: Optional molecular Vibrations object.
        :param soc: A list of spin_orbit_coupling.
        """
        super().__init__()
        self.metadata = metadata
        self.results = (self,)
        self.energies = energies
        self.pdm = pdm
        self.transition_dipole_moments = transition_dipole_moments
        self.atoms = atoms
        self.alignment = alignment
        self.orbitals = orbitals
        self.beta_orbitals = beta_orbitals
        self.ground_state = ground_state
        self.excited_states = excited_states
        self.energy_states = energy_states
        self.vibrations = vibrations
        self.vertical_emission = {}
        self.adiabatic_emission = {}
        self.soc = soc
        
    @property
    def CC_energies(self):
        warnings.warn("CC_energies is deprecated, use energies.cc instead", DeprecationWarning)
        return self.energies.cc
    
    @CC_energies.setter
    def CC_energies(self, value):
        warnings.warn("CC_energies is deprecated, use energies.cc instead", DeprecationWarning)
        self.energies.cc = value
    
    @property
    def MP_energies(self):
        warnings.warn("MP_energies is deprecated, use energies.mp instead", DeprecationWarning)
        return self.energies.mp
    
    @MP_energies.setter
    def MP_energies(self, value):
        warnings.warn("MP_energies is deprecated, use energies.mp instead", DeprecationWarning)
        self.energies.mp = value
    
    @property
    def SCF_energies(self):
        warnings.warn("SCF_energies is deprecated, use energies.scf instead", DeprecationWarning)
        return self.energies.scf
    
    @SCF_energies.setter
    def SCF_energies(self, value):
        warnings.warn("SCF_energies is deprecated, use energies.scf instead", DeprecationWarning)
        self.energies.scf = value
        
    @property
    def molecular_orbitals(self):
        warnings.warn("molecular_orbitals is deprecated, use orbitals instead", DeprecationWarning)
        return self.orbitals
    
    @molecular_orbitals.setter
    def molecular_orbitals(self, value):
        warnings.warn("molecular_orbitals is deprecated, use orbitals instead", DeprecationWarning)
        self.orbitals = value
        
    @property
    def spin_orbit_coupling(self):
        warnings.warn("spin_orbit_coupling is deprecated, use soc instead", DeprecationWarning)
        return self.soc
    
    @spin_orbit_coupling.setter
    def spin_orbit_coupling(self, value):
        warnings.warn("spin_orbit_coupling is deprecated, use soc instead", DeprecationWarning)
        self.soc = value
        
    @property
    def dipole_moment(self):
        warnings.warn("dipole_moment is deprecated, use pdm instead", DeprecationWarning)
        return self.pdm
    
    @dipole_moment.setter
    def dipole_moment(self, value):
        warnings.warn("dipole_moment is deprecated, use pdm instead", DeprecationWarning)
        self.pdm = value
        
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
        warnings.warn("energy is deprecated, use energies.final instead", DeprecationWarning)
        return self.energies.final
    
        
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
        
    def dump(self):
        return {
            "metadata": self.metadata.dump(),
            "atoms": self.atoms.dump(),
            "pdm": self.dipole_moment.dump(),
            "energies": self.energies.dump(),
            "ground_state": self.ground_state.dump(),
            "orbitals": self.orbitals.dump(),
            "beta_orbitals": self.beta_orbitals.dump(),
            "excited_states": self.excited_states.dump(),
            "soc": self.spin_orbit_coupling.dump(),
            "vibrations": self.vibrations.dump(),
            "adiabatic_emission": {key:value.dump() for key,value in self.adiabatic_emission.items()},
            "vertical_emission": {key:value.dump() for key,value in self.vertical_emission.items()}
        }
        
        