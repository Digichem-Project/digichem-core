import yaml
from silico.file.convert.gaussian import Gaussian_input_parser
from silico.file.convert.babel import Openbabel_converter
from pathlib import Path

# Custom formats to allow literal strings in yaml output.
# Adapted from https://stackoverflow.com/questions/6432605/any-yaml-libraries-in-python-that-support-dumping-of-long-strings-as-block-liter

# Custom 'str' class that will be dumped literally.
class literal_str(str): pass
# The dumper (which recognises the literal_str class).
def literal_str_representer(dumper, data):
	return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
# Add to yaml.
yaml.add_representer(literal_str, literal_str_representer)

class Silico_input():
	"""
	Class that represents an input file in the silico input (.si) format.
 	
	The .si format is yaml based, and stores the input geometry in a modified xyz format along with charge and multiplicity.
	"""
	
	def __init__(self,  geometry, charge = None, multiplicity = None, *, name = None):
		"""
		Constructor for .si files.
 		
		:param geometry: The molecular geometry in .si format. Use from_xyz() instead if your format is in xyz.
		:param charge: The molecular charge.
		:param multiplicity: The molecular multiplicity (as an integer).
		:param name: Name of the system/molecule
		"""
		self.geometry = geometry
		self.charge = charge
		self.multiplicity = multiplicity
		self.name = name

	@classmethod
	def from_xyz(self, geometry, *args, **kwargs):
		"""
		Create a Silico_input object from a molecule in xyz format.
		
		:param geometry: The input geometry in xyz format.
		:param charge: The molecular charge.
		:param multiplicity: The molecular multiplicity (as an integer).
		:param name: Name of the system/molecule
		"""
		# Split the xyz format on newlines.
		split_geom = geometry.split("\n")
		
		# Remove the first line (the number of atoms) and the second (optional comment line)
		split_geom = split_geom[2:]
		
		# Call our main constructor.
		return self("\n".join(split_geom), *args, **kwargs)
	
	@classmethod
	def from_com(self, geometry, charge = None, multiplicity = None, *args, **kwargs):
		"""
		Create a Silico_input object from a molecule in gaussian input format (.com, .gjf etc).
		
		:param geometry: The input geometry in gaussian format.
		:param charge: The molecular charge.
		:param multiplicity: The molecular multiplicity (as an integer).
		:param name: Name of the system/molecule
		"""
		# Get a parser for our gaussian input file.
		parser = Gaussian_input_parser(geometry)
		
		# If we've not been given a charge and/or multi, use the values from the input file.
		if charge is None:
			charge = parser.charge
		if multiplicity is None:
			multiplicity = parser.multiplicity
			
		# Continue with other constructors.
		return self.from_xyz(parser.xyz, charge = charge, multiplicity = multiplicity, *args, **kwargs)
		
		
	@classmethod
	def from_file(self, file_name, file_type = None, *args, gen3D = None, name = None, **kwargs):
		"""
		Create a Silico_input object from a file in arbitrary format.
		
		:param file_name: Name/path of the input file to read from.
		:param file_type: The format of the file; a string recognised by openbabel. If not given, an attempt will be made to guess from the file name (see Openbabel_converter.type_from_file_name()).
		:param charge: The molecular charge.
		:param multiplicity: The molecular multiplicity (as an integer).
		:param name: Name of the system/molecule
		"""
		file_name = Path(file_name)
		
		# We convert all formats to gaussian input formats (because this format contains charge and multiplicity, which we can extract).
		com_file = Openbabel_converter.from_file(file_name, file_type, gen3D = gen3D).convert("com")
		
		# Get a name from the file if none were given.
		if name is None:
			name = file_name.with_suffix("").name
		
		# Continue with other constructors.
		return self.from_com(com_file, *args, name = name, **kwargs)
		
	def to_file(self, file):
		"""
		Write this silico input file to an open file object.
		
		:param file: A file opened with open().
		"""
		file.write(self.yaml)
		
	def to_format(self, file_type):
		"""
		Get this input file in an arbitrary format.
		
		:param file_type: The format of the file; a string recognised by openbabel.
		"""
		return Openbabel_converter.get_cls("xyz")(input_file = self.xyz, input_file_type = "xyz", gen3D = False).convert(file_type)

	@property
	def implicit_charge(self):
		"""
		The charge of the molecule/system, accounting for cases where no explicit charge is set.
		""" 
		if self.charge is None:
			return 0
		else:
			return self.charge
		
	@property
	def implicit_multiplicity(self):
		"""
		The multiplicity (as an integer) of the molecule/system, accounting for cases where no explicit multiplicity is set.
		""" 
		if self.multiplicity is None:
			return 1
		else:
			return self.multiplicity
		
	@property
	def dict(self):
		"""
		Get this input file as a dict.
		"""
		return {
			'name': self.name,
			'charge': self.charge,
			'multiplicity': self.multiplicity,
			'geometry': self.geometry,
		}
		
	@property
	def yaml(self):
		"""
		Get this input file in yaml format.
		"""
		# Get in dict format.
		intermediate = self.dict
		
		# Wrap geometry in literal format (so it will appear line by line).
		intermediate['geometry'] = literal_str(intermediate['geometry'])
		
		# Convert.
		return yaml.dump(intermediate, sort_keys=False)
		
		
	@property
	def xyz(self):
		"""
		Get the geometry of this input file in XYZ format.
		"""
		return "{}\n\n{}".format(len(self.geometry.split("\n")), self.geometry) 
	