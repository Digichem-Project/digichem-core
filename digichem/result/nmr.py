import numpy

from silico.result.base import Result_object, Result_container


class NMR():
    """
    A result object containing all the NMR related data for a single atom.
    """
    
    def __init__(self, atom, shielding, couplings):
        """
        :param atom: The atom these NMR parameters relate to.
        :param shielding: The chemical shielding object for this atom.
        :param couplings: A dictionary of all the coupling interactions calculated for this atom.
        """
        self.atom = atom
        self.shielding = shielding
        self.couplings = couplings
    
    @classmethod
    def list_from_parser(self, parser):
        """
        """
        

class NMR_tensor_ABC(Result_object):
    """ABC for classes that contain dicts of NMR tensors."""
    
    tensor_names = ()
    
    def __init__(self, tensors, total_isotropic):
        self.tensors = tensors
        self.total_isotropic = total_isotropic
    
    def eigenvalues(self, tensor = "total"):
        """
        Calculate the eigenvalues for a given tensor.
        
        :param tensor: The name of a tensor to calculate for (see tensor_names). Use 'total' for the total tensor.
        """
        try:
            return numpy.linalg.eigvals(self.tensors[tensor])
        
        except KeyError:
            if tensor not in self.tensor_names:
                raise ValueError("The tensor '{}' is not recognised") from None
            
            elif tensor not in self.tensors:
                raise ValueError("The tensor '{}' is not available") from None
        
    def isotropic(self, tensor = "total"):
        """
        Calculate the isotropic value for a given tensor.
        
        :param tensor: The name of a tensor to calculate for (see tensor_names). Use 'total' for the total tensor.
        """
        eigenvalues = self.eigenvalues(tensor)
        return sum(eigenvalues) / len(eigenvalues)


class NMR_shielding(NMR_tensor_ABC):
    """
    A result object to represent the chemical shielding of an atom.
    """
    
    tensor_names = ("paramagnetic", "diamagnetic", "total")
    
    def __init__(self, atom, tensors, total_isotropic):
        """
        :param atom: The atom this shielding applies to.
        :param tensors: A dictionary of tensors.
        :param total_isotropic: The isotropic value for the total tensor.
        """
        super().__init__(tensors, total_isotropic)
        self.atom = atom
        
    @classmethod
    def list_from_parser(self, parser):
        """
        Create a list of NMR shielding objects from an output file parser.
        
        :param parser: An output file parser.
        :return: A list of NMR_shielding objects. The list will be empty if no NMR data is available.
        """
        shieldings = []
        try:
            for atom_index, tensors in parser.data.nmrtensors.items():
                total_isotropic = tensors.pop("isotropic")
                shieldings.append(self(parser.atoms[atom_index], tensors, total_isotropic))
        
        except AttributeError:
            return []
        
        return shieldings


class NMR_spin_couplings_list(Result_container):
    """A collection of NMR spin-spin couplings."""
    
    @classmethod
    def from_parser(self, parser):
        """
        Get an NMR_spin_couplings_list object from an output file parser.
        
        :param parser: An output file parser.
        :return: A NMR_spin_couplings_list object. The list will be empty if no NMR data is available.
        """
        return self(NMR_spin_coupling.list_from_parser(parser))


class NMR_spin_coupling(Result_object):
    """
    A result object to represent spin-spin NMR couplings.
    """
    
    tensor_names = ("paramagnetic", "diamagnetic", "fermi", "spin-dipolar", "spin-dipolar-fermi", "total")
    
    def __init__(self, atoms, isotopes, tensors, total_isotropic):
        """
        :param atoms: Tuple of atoms that this coupling is between.
        :param isotopes: Tuple of the specific isotopes of atoms.
        :param tensors: A dictionary of tensors.
        :param total_isotropic: The isotropic value for the total tensor.
        """
        super().__init__(tensors, total_isotropic)
        self.atoms = atoms
        self.isotopes = isotopes
        
    @classmethod
    def list_from_parser(self, parser):
        """
        Create a list of NMR coupling objects from an output file parser.
        
        :param parser: An output file parser.
        :return: A list of NMR_spin_coupling objects. The list will be empty if no NMR data is available.
        """
        couplings = []
        try:
            for atom_tuple, isotopes in parser.data.nmrcouplingtensors.items():
                for isotope_tuple, tensors in isotopes.items():
                    total_isotropic = tensors.pop("isotropic")
                    couplings.append(self((parser.atoms[atom_tuple[0]], parser.atoms[atom_tuple[1]]), isotope_tuple, tensors, total_isotropic))
        
        except AttributeError:
            return []
        
        return couplings