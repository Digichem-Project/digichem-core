import numpy
from itertools import filterfalse

from silico.result.base import Result_object, Result_container
from silico.exception.base import Result_unavailable_error


class NMR_list(Result_container):
    
    @classmethod
    def from_parser(self, parser):
        return self(NMR.list_from_parser(parser))
    
    def find(self, criteria = None, *, label = None, index = None):
        """
        """
        if label is None and index is None and criteria is None:
            raise ValueError("One of 'criteria', 'label' or 'index' must be given")
        
        if criteria is not None:
            if str(criteria).isdigit():
                index = int(criteria)
            
            else:
                label = criteria
        
        # Now get our filter func.
        if label is not None:
            filter_func = lambda nmr: nmr.atom.label != label
        
        elif index is not None:
            filter_func = lambda nmr: nmr.atom.index != index
        
        # Now search.
        found = type(self)(filterfalse(filter_func, self))
        
        try:
            return found[0]
        
        except IndexError:
            if label is not None:
                criteria_string = "label = '{}'".format(label)
            
            elif index is not None:
                criteria_string = "index = '{}'".format(index)
            
            raise Result_unavailable_error("NMR", "could not NMR data for atom '{}'".format(criteria_string))


class NMR(Result_object):
    """
    A result object containing all the NMR related data for a single atom.
    
    For a given atom, this class will containg:
        - The chemical shielding of this atom (including a breakdown by tensors, if available).
        - The spin-spin coupling constants between this atom and all other atoms for which couplings were calculated (also with a breakdown by tensor, if available).
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
        return list([self(atom, parser.results.nmr_shieldings[atom], parser.results.nmr_couplings.between(atom)) for atom in parser.results.atoms if atom in parser.results.nmr_shieldings])
        
    def dump(self, silico_options):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "atom": self.atom.label,
            "shielding": self.shielding.dump(silico_options)
        }
        

class NMR_tensor_ABC(Result_object):
    """ABC for classes that contain dicts of NMR tensors."""
    
    tensor_names = ()
    
    def __init__(self, tensors, total_isotropic):
        self.tensors = tensors
        self.total_isotropic = total_isotropic
    
    def eigenvalues(self, tensor = "total", real_only = True):
        """
        Calculate the eigenvalues for a given tensor.
        
        :param tensor: The name of a tensor to calculate for (see tensor_names). Use 'total' for the total tensor.
        """
        try:
            return numpy.array([val.real for val in numpy.linalg.eigvals(self.tensors[tensor])])
        
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
    
    def dump(self, silico_options):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "tensors": {t_type: {"value": list(list(map(float, dim)) for dim in tensor), "units": "ppm"} for t_type, tensor in self.tensors.items()},
            "eigenvalues": {t_type: {"value": list(list(map(float, self.eigenvalues(t_type)))), "units": "ppm"} for t_type in self.tensors},
            "isotropic": {t_type: {"value": float(self.isotropic(t_type)), "units": "ppm"} for t_type in self.tensors}
        }


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
    def dict_from_parser(self, parser):
        """
        Create a dict of NMR shielding objects from an output file parser.
        Each key is the atom being shielded.
        
        :param parser: An output file parser.
        :return: A list of NMR_shielding objects. The list will be empty if no NMR data is available.
        """
        shieldings = {}
        try:
            for atom_index, tensors in parser.data.nmrtensors.items():
                total_isotropic = tensors.pop("isotropic")
                shieldings[parser.results.atoms[atom_index]] = self(parser.results.atoms[atom_index], tensors, total_isotropic)
        
        except AttributeError:
            return {}
        
        return shieldings
    
    def dump(self, silico_options):
        """
        Get a representation of this result object in primitive format.
        """
        dump_dic = {
            "atom": self.atom.label,
        }
        
        dump_dic.update(super().dump(silico_options))
        return dump_dic


# We could look at some more advanced type of container for couplings.
# A dictionary might make sense, but it's difficult to choose how exactly the keys should be arranged and ordered.
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
    
    def between(self, atom1, atom1_isotopes = None, atom2 = None, atom2_isotopes = None):
        """
        Return a list containing all couplings involving either one or two atoms.
        """
#         new_list = self()
#         
#         for coupling in self:
#             if atom1 in coupling.atoms and
#             (atom2 is None or atom2 in coupling.atoms) and
#             (atom1_isotopes is None or coupling.isotopes[coupling.atoms.index(atom1)] in atom1_isotopes) and
#             (atom2_isotopes is None or atom2 is None or coupling.isotopes[coupling.atoms.index(atom2)] in atom2_isotopes)
#             
#             :
#                 # Decide if any of our isotopes match.
#                 atom1_index
        return type(self)([coupling for coupling in self
            if atom1 in coupling.atoms and
            (atom2 is None or atom2 in coupling.atoms) and
            (atom1_isotopes is None or coupling.isotopes[coupling.atoms.index(atom1)] in atom1_isotopes) and
            (atom2_isotopes is None or atom2 is None or coupling.isotopes[coupling.atoms.index(atom2)] in atom2_isotopes)
        ])

class NMR_spin_coupling(NMR_tensor_ABC):
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
                    couplings.append(self((parser.results.atoms[atom_tuple[0]], parser.results.atoms[atom_tuple[1]]), isotope_tuple, tensors, total_isotropic))
        
        except AttributeError:
            return []
        
        return couplings