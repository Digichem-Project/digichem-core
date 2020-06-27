from pathlib import Path
from datetime import datetime

class Silico_Directory():
	"""
	Top level class for Directory helper classes.
	
	Sadly, we can't subclass pathlib Paths yet.
	"""
	
	# Characters that are considered unsafe and are replaced with underscores.
	UNSAFE_CHARS = ["/", "\\"]

	def __init__(self, *args, **kwargs):
		"""
 		Constructor for Calculation_directory objects.
 		"""
		
		self.path = Path(*args, **kwargs)
		
	def __str__(self):
		return str(self.path)
	
	@classmethod
	def safe_name(self, file_name, unsafe_chars = None):
		"""
		Get a safe version of a file name.
		
		Note that both '/' and '\' are considered unsafe, so supplying paths to this function is probably not what you want to do.
		
		:param file_name: The unsafe file_name.
		:param unsafe_chars: A list of characters to replace. If None is supplied, the default list is used.
		"""
		if unsafe_chars is None:
			unsafe_chars = self.UNSAFE_CHARS
		
		safe_str = file_name
		for unsafe_char in self.UNSAFE_CHARS:
			safe_str = safe_str.replace(unsafe_char, "_")
		return safe_str

class Molecule_directory(Silico_Directory):
	"""
	Class that represents directory holding several calculations for a single molecule.
	"""		
		
	@classmethod
	def from_calculation(self, calculation):
		"""
		Create a Molecule_directory object from a Calculation_target.
		"""
		return self(calculation.output, calculation.name)
	
	def create_structure(self):
		"""
		Create the directory structure of this Molecule_directory object.
		
		:return: True normally, False if the structure already exists.
		"""
		try:
			self.path.mkdir(parents = True)
		except FileExistsError:
			return False
		return True
		
	def get_calculation_dirs(self):
		"""
		Get a list of Calculation_directory objects that are in this molecule dir.
		"""
		return [Calculation_directory(calc_path) for calc_path in self.path.iterdir()]

class Calculation_directory(Silico_Directory):
	"""
	Class that represents the directory hierarchy where we submit calculations to.
	"""
	
	def __init__(self, molecule_directory, name, directory_time = None, create = False):
		"""
		Constructor for 
		"""
		self.molecule_directory = molecule_directory
		self.name = name
		self.directory_time = directory_time if directory_time is not None else datetime.now()
		if create:
			self.molecule_directory.create_structure()
			self.create_structure(True)
			
	@classmethod
	def from_calculation(self, calculation, create = False):
		"""
		Create a Calculation_directory object from a Calculation_target object.
		"""
		return self(Molecule_directory.from_calculation(calculation), self.safe_name(calculation._CONFIG_NAME), create = create)		
	
	@property
	def path(self):
		"""
		Pathlib Path object to this calculation dir.
		"""
		#return Path(str(self.molecule_directory), self.name, self.directory_time.strftime("%d-%m-%Y %M-%H-%S"))
		return Path(str(self.molecule_directory), self.name)
	
	@property
	def input_directory(self):
		"""
		Full path to the calculation input directory.
		"""
		return Path(str(self) + "/Input")
	
	@property
	def output_directory(self):
		"""
		Full path to the calculation output directory.
		"""
		return Path(str(self) + "/Output")
	
	@property
	def result_directory(self):
		"""
		Full path to the results directory.
		"""
		return Path(str(self) + "/Results")
	
	@property
	def scratch_directory(self):
		"""
		Full path to the scratch directory.
		
		Note that this directory is not the scratch directory as written to by quantum chemistry programs, but rather where we attempt to save the scratch in cases of error etc.
		"""
		return Path(str(self) + "/Scratch")
	
	@property
	def report_directory(self):
		"""
		Full path to the calculation report directory.
		"""
		return Path(str(self) + "/Report")
	
	def create_structure(self, adjust = False):
		"""
		Create the directory structure of this Calculation_directory object.
		
		:return True normally, False if the directory structure has already been made.
		"""
		retval = 0
		
		# First try and make our base directory.
		counter = 1
		dir_name = self.name
		while True:
			try:
				self.name = dir_name + " {}".format(str(counter).zfill(2)) if counter != 1 else dir_name
				self.path.mkdir(parents = True)
				break
			except FileExistsError:
				if adjust:
					counter +=1
				else:
					raise
		
		
		for sub_dir in [self.input_directory, self.output_directory]:
			try:
				sub_dir.mkdir(parents = True)
			except FileExistsError:
				retval += 1
		return retval != 3
		
		
		