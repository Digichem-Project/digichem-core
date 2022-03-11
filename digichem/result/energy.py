# Classes for total system energies.
from silico.result import Result_container
import scipy.constants
from silico.exception.base import Result_unavailable_error
from silico.result import Unmergeable_container_mixin

class Energy_list(Result_container, Unmergeable_container_mixin):
    """
    Class that represents a list of calculated energies. Storing energies as a list is useful because optimisations will print energies at each step, and it can be useful to look back and see how the opt has progressed.
    """
    
    # The type of energy; corresponds to the name of the attribute provided by cclib.
    cclib_energy_type = ""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
    def __float__(self):
        """
        Our float representation, which is the last energy in our list.
        """
        return self.final
    
    def __str__(self):
        """
        Our string representation, same as the float but in text form.
        """
        return str(self.__float__())
    
    def assign_levels(self):
        """
        (Re)assign total levels of the objects in this list.
        
        The contents of this list will be modified in place.
        """
        # Energies don't have levels.
        pass
    
    def human_energy_type(self):
        """
        Human readable name of this type of energy (for example, self-consistent field).
        """
        if self.energy_type == "SCF":
            return "self-consistent field"
        elif self.energy_type == "MP":
            return "Møller–Plesset"
        elif self.energy_type == "CC":
            return "coupled-cluster"
        else:
            raise ValueError("Unrecognised energy type {}".format(self.energy_type))

    @property
    def final(self):
        """
        The 'final' energy in this list, useful for getting the final optimised energy.
        
        :raises IndexError: If this list is empty.
        """
        try:
            return self[-1]
        except IndexError:
            raise Result_unavailable_error(self.energy_type + ' energy', 'there is no {} energy'.format(self.energy_type))
    
    @classmethod
    def eV_to_kJ(self, eV):
        """
        Convert a value in eV to kJ.
        
        :param eV: The value in electron volts.
        :return: The value in kilojoules
        """
        return eV * 1.602176634e-22
    
    @classmethod
    def kJ_to_kJmol(self, kJ):
        """
        Convert a value in kJ to kJmol-1.
        
        :param kJ: The value in kilojoules.
        :return: The value in kilojoules per mol
        """
        return kJ * scipy.constants.Avogadro
    
    @classmethod
    def eV_to_kJmol(self, eV):
        """
        Convert a value in eV to kJmol-1.
        
        :param eV: The value in electron volts.
        :return: The value in kilojoules per mol
        """
        return self.kJ_to_kJmol(self.eV_to_kJ(eV))

    @classmethod
    def from_parser(self, parser):
        """
        Get an Energy_list from an output file parser.
        
        :param parser: An output file parser.
        :return: The populated Energy_list object. The object will be empty if the energy is not available.
        """
        try:
            return self(getattr(parser.data, self.cclib_energy_type))
        except AttributeError:
            return self()
    
class SCF_energy_list(Energy_list):
    """
    List of Self-consistent field energies.
    """
    cclib_energy_type = "scfenergies"
    energy_type = "SCF"
    
    # A warning issued when attempting to merge non-equivalent objects
    MERGE_WARNING = "Attempting to merge list of SCF energies that are not identical; non-equivalent energies will be ignored"
    
class CC_energy_list(Energy_list):
    """
    List of coupled-cluster energies.
    """
    cclib_energy_type = "ccenergies"
    energy_type = "CC"
    
    # A warning issued when attempting to merge non-equivalent objects
    MERGE_WARNING = "Attempting to merge list of CC energies that are not identical; non-equivalent energies will be ignored"
    
class MP_energy_list(Energy_list):
    """
    List of Moller-Plesset energies.
    """
    cclib_energy_type = "mpenergies"
    energy_type = "MP"
    
    # A warning issued when attempting to merge non-equivalent objects
    MERGE_WARNING = "Attempting to merge list of MP energies that are not identical; non-equivalent energies will be ignored"
    
    @classmethod
    def from_parser(self, parser, order = -1):
        """
        Get an MP_energy_list from an output file parser.
        
        Note that unlike other calculated energies, MP energies are calculated sequentially up to the requested order. So for example, MP4 first calculates the MP1 energy (which is the same as the uncorrected energy), then MP2, MP3 and finally MP4.
        
        :param parser: An output file parser.
        :param order: The order of the MP correction to get (ie, what is n in MPn). A value of -1 will get the highest order MP energy.
        :return: The populated Energy_list object. The object will be empty if the energy is not available.
        """
        try:
            return self([energy[order] for energy in super().from_parser(parser)])
        except IndexError:
            # No energy.
            return self()
            
    
    
    