import math
from silico.result.base import Result_object
from silico.result.angle import Angle
from silico.file.cube import Cube_maker
from pathlib import Path
from silico.image.vmd import Dipole_image_maker

class Dipole_moment(Result_object):
	"""
	Class that represents a dipole moment.
	"""	
	
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
		return math.sqrt( (self.vector_coords[0] - self.origin_coords[0]) ** 2 + (self.vector_coords[1] - self.origin_coords[1]) ** 2 + (self.vector_coords[2] - self.origin_coords[2]) ** 2)
	
	
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
		self.atoms = atoms
		
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
		if self.atoms is not None:
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
		if self.atoms is not None:
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
	
	@property
	def X_axis_angle(self):
		"""
		The angle between this dipole moment and the X axis of the atom set of the molecule.
		"""
		if self.total == 0:
			return Angle(0)
		elif self.atoms is not None:
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
		elif self.atoms is not None:
			return Angle(math.fabs(self.atoms.get_XY_plane_angle(self.origin_coords, self.vector_coords)), "rad")
		else:
			return None
			
	def set_file_options(self, output_dir, output_name, cube_file = None, **kwargs):
		"""
		Set the options that will be used to create images from this object.
		
		:param output_dir: A pathlib Path object to the directory within which our files should be created.
		:param output_name: A string that will be used as the start of the file name of the files we create.
		"""
		# Get ourselves a cube file maker if we need one.
		if cube_file is None:
			cube_file = Cube_maker.from_image_options(Path(output_dir, "Dipole Moment", output_name + ".dipole.cube"), cubegen_type = "MO", orbital = "HOMO", **kwargs)
			
		# Get our image.
		self._files['dipole_image'] = Dipole_image_maker.from_image_options(Path(output_dir, "Dipole Moment", output_name + ".dipole.jpg"), cube_file = cube_file, dipole_moment = self, **kwargs)
	
	def cleanup_intermediate_files(self):
		"""
		Remove any intermediate files that may have been created by this object.
		"""
		# Remove our cube file.
		super().cleanup_intermediate_files('cube_file')
	
	@property
	def dipole_image(self):
		return self.get_file('dipole_image')
	
# 	def realign(self, atoms):
# 		"""
# 		Realign this dipole moment to a new coordinate system.
# 		
# 		:param atoms: An atom Alignment object.
# 		"""
# 		# Realign both our coordinates.
# 		self.origin_coords = atoms.apply_transformation(self.origin_coords)
# 		self.vector_coords = atoms.apply_transformation(self.vector_coords)
# 		# Done.
		
	@classmethod
	def from_cclib(self, ccdata, atoms = None):
		"""
		Construct a Dipole_moment object from the data provided by cclib.
		
		:param ccdata: Result object as provided by cclib.
		:param atoms: Optional Atom Alignment object that will be used to calculate addition data and to optionally realign the coordinates of our dipole.
		:result: A single Dipole_moment object, or None if no dipole information is available.
		"""
		try:
			dipole = self(ccdata.moments[0], ccdata.moments[1], atoms)
			return dipole
		except AttributeError:
			return None

	
	
	