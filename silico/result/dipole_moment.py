# General imports
import math

# Silico imports
from silico.result import Result_object
from silico.result.angle import Angle

class Dipole_moment(Result_object):
    """
    Class that represents a dipole moment.
    """
    
    # A warning issued when attempting to merge non-equivalent objects
    MERGE_WARNING = "Attempting to merge list of DPMs that are not identical; non-equivalent DPMs will be ignored"
    
    @property
    def magnitude(self):
        """
        The dipole moment magnitude in debye.
        
        :deprecated: 'magnitude' is an inaccurate name for this attribute, use 'total' instead.
        """
        return self.total
        
    
    @property
    def total(self):
        """
        The dipole moment in debye.
        """
        #return math.sqrt( (self.vector_coords[0] - self.origin_coords[0]) ** 2 + (self.vector_coords[1] - self.origin_coords[1]) ** 2 + (self.vector_coords[2] - self.origin_coords[2]) ** 2)
        return math.sqrt( self.vector_coords[0] ** 2 + self.vector_coords[1] ** 2 + self.vector_coords[2] ** 2)
    
    
    def __init__(self, origin_coords, vector_coords, atoms = None):
        """
        Constructor for Dipole_moment objects.
        
        :param origin_coords: Tuple of (x, y, z) coordinates of where this dipole moment starts.
        :param vector_coords: Tuple of (x, y, z) coordinates of where this dipole moment ends.
        :param atoms: Optional atoms atom list object which is used to calculate addition data about this dipole moment. 
        """
        super().__init__()
        # The start of our vector, normally (0,0,0).
        self.origin_coords = origin_coords
        # The end of our vector.
        self.vector_coords = vector_coords
        
        # Save our atoms object.
        self.atoms = atoms if atoms is not None else []
        
        # Save a name describing which dipole we are (permanent vs transition etc).
        self.dipole_type = "permanent"
    
    @property
    def name(self):
        """
        Name that describes this dipole moment.
        """
        return "PDM"
    
    @property
    def origin_coords(self):
        """
        The origin coords of this vector as a tuple of (x, y, z). origin_coords is automatically realigned by the atoms alignment object of this dipole, use _origin_coords if you do not want this behaviour.
        """
        if len(self.atoms) > 0:
            return self.atoms.apply_transformation(self._origin_coords)
        else:
            return self._origin_coords
        
    @origin_coords.setter
    def origin_coords(self, value):
        """
        Set the origin coords of this dipole.
        """
        self._origin_coords = value
        
    @property
    def vector_coords(self):
        """
        The ending coords of this vector as a tuple of (x, y, z). vector_coords is automatically realigned by the atoms alignment object of this dipole, use _vector_coords if you do not want this behaviour.
        """
        if len(self.atoms) > 0:
            return self.atoms.apply_transformation(self._vector_coords)
        else:
            return self._vector_coords
        
    @vector_coords.setter
    def vector_coords(self, value):
        """
        Set the ending coords of this dipole.
        """
        self._vector_coords = value
        
    def __float__(self):
        """
        Floatify this Dipole_moment.
        """
        return self.magnitude
    
    def __eq__(self, other):
        """
        Is this dipole moment equal to another?
        
        A dipole moment is considered equivalent only if all coordinates match.
        """
        return self.origin_coords == other.origin_coords and self.vector_coords == other.vector_coords
    
    @property
    def X_axis_angle(self):
        """
        The angle between this dipole moment and the X axis of the atom set of the molecule.
        """
        if self.total == 0:
            return Angle(0)
        elif len(self.atoms) > 0:
            return Angle(math.fabs(self.atoms.get_X_axis_angle(self.origin_coords, self.vector_coords)), "rad")
        else:
            return None
        
    @property
    def XY_plane_angle(self):
        """
        The angle between this dipole moment and the XY plane of the atom set of the molecule.
        """
        if self.total == 0:
            return Angle(0)
        elif len(self.atoms) > 0:
            return Angle(math.fabs(self.atoms.get_XY_plane_angle(self.origin_coords, self.vector_coords)), "rad")
        else:
            return None
            
    @classmethod
    def from_parser(self, parser):
        """
        Construct a Dipole_moment object from an output file parser.
        
        :param parser: An output file parser.
        :result: A single Dipole_moment object, or None if no dipole information is available.
        """
        try:
            return self(parser.data.moments[0], parser.data.moments[1], parser.results.alignment)
        except AttributeError:
            return None
    
    