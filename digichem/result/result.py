# Code for processing results from calculations and generating reports.

# General imports.
import itertools
import warnings

# Silico imports.
from silico.exception import Result_unavailable_error
from silico.result.alignment.AAA import Adjusted_average_angle
from silico.result.alignment.AA import Average_angle
from silico.result.alignment.FAP import Furthest_atom_pair
from silico.result.alignment import Minimal
from silico.result import Result_object
import silico.result.excited_state
from silico.result.emission import Emissions


class Result_set(Result_object):
    """
    Class that represents a collection of results from a calculation.
    
    This class is a bit heavy and might not be around for long...
    """
    
    def __init__(self, **attributes):
        """Constructor for Result_set objects."""
        super().__init__()
        self._id = attributes.pop('database_id', None)
        self.results = (self,)
        self.emission = attributes.pop('emission', Emissions())
        
        for attr_name, attribute in attributes.items():
            setattr(self, attr_name, attribute)
        
    @property
    def alignment(self):
        warnings.warn("alignment is deprecated, use atoms instead", DeprecationWarning)
        return self.atoms

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
    def vertical_emission(self):
        warnings.warn("vertical_emission is deprecated, use emission.vertical instead", DeprecationWarning)
        return self.emission.vertical
    
    @vertical_emission.setter
    def vertical_emission(self, value):
        warnings.warn("vertical_emission is deprecated, use emission.vertical instead", DeprecationWarning)
        self.emission.vertical = value
    
    @property
    def adiabatic_emission(self):
        warnings.warn("adiabatic_emission is deprecated, use emission.adiabatic instead", DeprecationWarning)
        return self.emission.adiabatic
    
    @adiabatic_emission.setter
    def adiabatic_emission(self, value):
        warnings.warn("adiabatic_emission is deprecated, use emission.adiabatic instead", DeprecationWarning)
        self.emission.adiabatic = value
        
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
        
    def dump(self, silico_options):
        "Dump the data contained in this result set, serialising it to a hierarchy of dicts that can be saved in various formats."
        # Start with our DB ID if we have one.
        dump_dic = {}
        if self._id is not None:
            dump_dic['_id'] = self._id
            
        dump_dic.update({
            "metadata": self.metadata.dump(silico_options),
            "ground_state": self.ground_state.dump(silico_options) if self.ground_state is not None else None,
            "energies": self.energies.dump(silico_options),
            "atoms": self.atoms.dump(silico_options),
            "raw_atoms": self.raw_atoms.dump(silico_options),
            "orbitals": self.orbitals.dump(silico_options),
            "beta_orbitals": self.beta_orbitals.dump(silico_options),
            "pdm": self.pdm.dump(silico_options) if self.pdm is not None else None,
            "excited_states": self.excited_states.dump(silico_options),
            "soc": self.soc.dump(silico_options),
            "vibrations": self.vibrations.dump(silico_options),
            "emission": self.emission.dump(silico_options)
        })
        
        return dump_dic