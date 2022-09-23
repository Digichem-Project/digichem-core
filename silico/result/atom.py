# General imports
import math
import periodictable
from itertools import zip_longest
from openbabel import pybel

# Silico imports
from silico.result import Result_container
from silico.result import Result_object
from silico.result import Unmergeable_container_mixin
from silico.exception.base import Result_unavailable_error


class Atom_list(Result_container, Unmergeable_container_mixin):
    """
    Class for representing a group of atoms.
    """
    
    # A warning issued when attempting to merge non-equivalent atom lists.
    MERGE_WARNING = "Attempting to merge lists of atoms that are not identical; non-equivalent atoms will be ignored"
    
    def __init__(self, *args, charge = None, **kwargs):
        super().__init__(*args, **kwargs)
        self.charge = charge if charge is not None else 0
        
    @property
    def mass(self):
        """
        The total mass of all the atoms in this set.
        
        :return: The mass (in Daltons).
        """
        try:
            return sum([atom.mass for atom in self])
        except TypeError:
            # Exact mass not available.
            raise Result_unavailable_error("Exact mass") from None
        
    @property
    def molar_mass(self):
        """
        The molar mass of the molecule (takes into account different isotopes and relative isotope abundances, unlike the mass attribute).
        
        :return: The mass (in Daltons).
        """
        return self.formula.mass
        
    @property
    def element_dict(self):
        """
        Get a dictionary where each key is one of the elements in the atom list (C, H, N etc) and the value is the number of that element that appears in the atom list.
        
        :return: The element dictionary.
        """
        atoms = {}
        for atom in self:
            # Try and increment the count of the atom.
            try:
                atoms[atom.element.symbol] += 1
            except KeyError:
                # Add the new atom.
                atoms[atom.element.symbol] = 1
        return atoms;
    
    @property
    def formula(self):
        """
        Get a formula representation of this atom list.
        
        :return: The formula as a periodictable.formula object (which can be safely cast to string).
        """
        # A dictionary where each key is a type of atom (N, C, H etc) and the value is the number of that atom.
        atoms = self.element_dict
        # Build a string rep.
        form_string = ""
        for atom in atoms:
            form_string += "{}{}".format(atom, atoms[atom])
        # Add our charge.
        #if self.charge
        #form_string += "{{{0}}}".format(self.charge)
        # Get and return the formula object.
        return periodictable.formula(form_string)
        #return periodictable.formula([atoms[key] for key in atoms])
        
    @property
    def formula_string(self):
        """
        Get a formula representation of this atom list as a string, including optional charge.
        """
        # Get the base formula
        formula_string = str(self.formula)
        
        # Add charge, if we have one.
        if self.charge == 1:
            formula_string += " +"
        elif self.charge == -1:
            formula_string += " -"
        elif self.charge != 0:
            formula_string += " {}{}".format(abs(self.charge), "-" if self.charge < 0 else "+")
            
        return formula_string
    
    @property
    def smiles(self):
        """
        Get this geometry in (canonical) SMILES format.
        
        This property uses the pybel smiles algorithm.
        """
        molecule = pybel.readstring("xyz", self.to_xyz())
        return molecule.write("can").strip()
        
    @property
    def X_length(self):
        return self.get_axis_length(0)
    
    @property
    def Y_length(self):
        return self.get_axis_length(1)
    
    @property
    def Z_length(self):
        return self.get_axis_length(2)
    
    def get_axis_length(self, axis):
        """
        Calculate the length of an axis, defined as the distance required in that axis to contain all the atoms of the set.
        
        :param axis: The axis to calculate for as an integer (0: X-axis, 1: Y-axis, 2: Z-axis).
        :return The length (in angstroms).
        """
        if not 0 <= axis <= 2:
            # Axis is invalid.
            raise ValueError("Axis '{}' is out of bounds. Possible  values are 0 (X), 1 (Y) or 2 (Z)")
        
        # First sort our list of atoms in terms of x, y or z coord.
        sorted_atoms = sorted(self, key = lambda atom: atom.coords[axis])
                
        # Now the axis length is simply the difference between the greatest and the smallest.
        return sorted_atoms[-1].coords[axis] - sorted_atoms[0].coords[axis]
    
    def get_linear_ratio(self):
        """
        Get the linear ratio of the molecule.
        
        The linear ratio is defined as 1 - (Y_length / X_length).
        
        :return: The ratio, from 0 (non-linear) to 1 (linear).
        """
        return 1- (self.Y_length / self.X_length)
    
    def get_planar_ratio(self):
        """
        Get the planar ratio of the molecule.
        
        The planar ratio is defined as 1 - (Z_length / Y_length).
        
        :return: The ratio, from 0 (non-planar) to 1 (planar).
        """
        
        return 1- (self.Z_length / self.Y_length)
    
    def get_X_axis_angle(self, start_coord, end_coord):
        """
        Get the angle between a line and the X axis.
        
        :param start_coord: A (X, Y, Z) tuple of coordinates of the start of the line.
        :param end_coord: A (X, Y, Z) tuple of coordinates of the end of the line.
        :return: The angle (in radians). 
        """
        return self.get_theta(math.sqrt( (end_coord[2] - start_coord[2])**2 + (end_coord[1] - start_coord[1])**2 ), end_coord[0] - start_coord[0])
    
    def get_XY_plane_angle(self, start_coord, end_coord):
        """
        Get the angle between a line and the XY plane.
        
        :param start_coord: A (X, Y, Z) tuple of coordinates of the start of the line.
        :param end_coord: A (X, Y, Z) tuple of coordinates of the end of the line.
        :return: The angle (in radians). 
        """
        # The 'secondary' axis is the opposite side of our triangle.
        secondary_axis = end_coord[2] - start_coord[2]
        # The 'primary' axis is the adjacent side of our triangle, which we can get with pythagoras.
        primary_axis = math.sqrt( (end_coord[0] - start_coord[0])**2 + (end_coord[1] - start_coord[1])**2 )
        return self.get_theta(secondary_axis, primary_axis)
    
    def dump(self, silico_options):
        dump_dict = {
            "formula": self.formula_string,
            "smiles": self.smiles,
            "exact_mass": {
                "value": float(self.mass) if self.safe_get("mass") is not None else None,
                "units": "g mol^-1" 
            },
            "molar_mass": {
                "value": self.molar_mass,
                "units": "g mol^-1",
                },
            "num_atoms": len(self),
            "x-extension": {
                "value": float(self.X_length),
                "units": "Å"
            },
            "y-extension": {
                "value": float(self.Y_length),
                "units": "Å"
            },
            "z-extension": {
                "value": float(self.Z_length),
                "units": "Å"
            },
            "linearity_ratio": float(self.get_linear_ratio()),
            "planarity_ratio": float(self.get_planar_ratio()),
            "values": super().dump(silico_options),
        }
        return dump_dict
    
    @classmethod
    def from_parser(self, parser):
        """
        Get an Atom_list object from an output file parser.
        
        :param parser: An output file parser.
        :param charge: Charge of the system.
        :return: A list of TDM objects.
        """
        return self(Atom.list_from_parser(parser), charge = parser.results.metadata.charge)
    
    @classmethod
    def merge(self, *multiple_lists, charge):
        """
        Merge multiple lists of atoms into a single list.
         
        Note that it does not make logical sense to combine different list of atoms into one; hence the method only ensures that all given lists (which are not empty) are the same and then returns the first (non empty) given.
        If the atom lists are not equivalent, a warning will be issued.
        """
        return super().merge(*multiple_lists, charge = charge)
    
    def to_xyz(self):
        """
        Convert this list of atoms to xyz format.
        """
        # First, the number of atoms.
        xyz = "{}\n\n".format(len(self))
        
        # Then coordinates.
        # No effort is made here to truncate coordinates to a certain precision.
        for atom in self:
            xyz += "{}    {}    {}    {}\n".format(atom.element.symbol, atom.coords[0], atom.coords[1], atom.coords[2])
            
        return xyz


class Atom(Result_object):
    """
    Class that represents an atom.
    """
    
    def __init__(self, atomic_number, coords, mass = None):
        """
        Construct an Atom class.
        
        :param atomic_number: The atomic/proton number of the atom. Conventional wisdom suggests this has to be an integer, but we make no check here.
        :param mass: The mass of the atom (in Daltons). This isn't always available for some reason (eg, Gaussian Freq calculations), but when it is it identifies the isotope of the given atom.
        :param coords: The coords of the atom, as an (x, y, z) tuple.        
        """
        # Just save each of our attributes.
        self.mass = mass
        self.coords = coords
        # Get our element class.
        self.element = periodictable.elements[atomic_number]
        
    def __str__(self):
        """
        Stringify this atom.
        """
        return "'{}' at '{}'".format(self.element, self.coords)
        
    def __eq__(self, other):
        """
        Is this atom equal to another?
        """
        # Atoms are considered equal if they are the same element in the same position.
        return self.element == other.element and self.coords == other.coords
            
    def distance(self, foreign_atom):
        """
        Get the distance between this atom and another atom.
        
        :return: The distance. The units depend on the units of the atoms' coordinates. If the two atoms have coordinates of different units, then you will obviously get bizarre results.
        """
        
        return math.sqrt( (self.coords[0] - foreign_atom.coords[0])**2 + (self.coords[1] - foreign_atom.coords[1])**2 + (self.coords[2] - foreign_atom.coords[2])**2)
    
    def dump(self, silico_options):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "element": self.element.number,
            "coords": {
                "x": {
                    "value": float(self.coords[0]),
                    "units": "Å", 
                },
                "y": {
                    "value": float(self.coords[1]),
                    "units": "Å", 
                },
                "z": {
                    "value": float(self.coords[2]),
                    "units": "Å", 
                }
            },
            "mass": {
                "value": float(self.mass) if self.mass is not None else None,
                "units": "g mol^-1"
            }
        }
        
    @classmethod
    def list_from_dump(self, data, result_set):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return [self(atom_dict['element'], atom_dict['coords']['value'], atom_dict['mass']['value']) for atom_dict in data['atoms']]
    
    @classmethod
    def list_from_parser(self, parser):
        """
        Get a list of Atom objects from an output file parser.
        
        :param parser: An output file parser.
        :result: A list of Atom objects. An empty list is returned if no atom data is available.
        """
        # First pack our data together to make is easier to loop through.
        try:
            atomnos = parser.data.atomnos
            # Atom coords contains a list for each iteration, we only want the last.
            atomcoords = parser.data.atomcoords[-1]
            atommasses = getattr(parser.data, 'atommasses', [])
            
            # Atommasses sometimes is longer than atomcoords or atomnos.
            # This might be a cclib bug, it is not clear.
            if len(atommasses) > len(atomnos):
                # Take the last values only.
                atommasses = atommasses[-len(atomnos):]
            
        except AttributeError:
            # No atom data available.
            return []
        
        # Zip.
        zip_data = zip_longest(atomnos, atomcoords, atommasses, fillvalue = None)
        
        # Loop through and rebuild our objects.
        return [self(atomic_number, tuple(coords), mass) for atomic_number, coords, mass in zip_data]
        
        