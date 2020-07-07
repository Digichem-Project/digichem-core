import pkg_resources
from pathlib import Path
import yaml
import glob

from silico.config.base import Config
from silico.config.main import Silico_options
from silico.config.configurable import Configurable

class Config_loader(list):
	"""
	Class for reading config options from a YAML config file.
	
	Each item in this list is a dict, each key in these dicts is a key from a YAML document (note that each YAML file can contain several YAML documents).
	"""
	
	# The directory that contains the top-level config file (which contains final default values if all else fails).
	@classmethod
	def DEFAULT_CONFIG_LOCATION(self):
		return Path(pkg_resources.resource_filename('silico', "data/config"))
	DEFAULT_CONFIG_NAME = "default.yaml"
	
	# The directory that contains the master config.
	@classmethod
	def MASTER_CONFIG_LOCATION(self):
		return Path(pkg_resources.resource_filename('silico', "data/config"))
	MASTER_CONFIG_NAME = "silico.yaml"
	
	# The name of the directory in which all found yaml files will be loaded.
	MULTI_DIRECTORY_NAME = "conf.d"
	
	# The location of the system-wide config file (which has precedence over the master).
	@classmethod
	def SYSTEM_CONFIG_LOCATION(self):
		return Path("/etc/silico")
	
	# The location of the user specific config file (which has precedence over the system).
	@classmethod
	def USER_CONFIG_LOCATION(self):
		return Path(Path.home(), ".config/silico") 
	
	def __init__(self, config_file, *, file_path = None):
		"""
		Constructor for config file loaders.
		
		:param config_file: A loaded config file as a string (or an open file object to read from).
		"""
		super().__init__()
		
		self.file_path = file_path
		
		# Load our options if we've been asked to.
		if config_file is not None:
			self.load(config_file)
			
	@classmethod
	def from_file(self, config_file_path):
		"""
		Alternative constructor that first opens a config file path before reading it.
		"""
		with open(config_file_path, "rt") as config_file:
			return self(config_file, file_path = config_file_path)
		
	@classmethod
	def silico(self):
		"""
		Load all Silico config files.
		
		Config files are searched for in the following order (increasing precedence) 1) the silico install directory, 2) /etc/silico, 3) the user's ~/.config/silico directory
		
		Note that configurables will be added to the returned Silico_options, but they will not be resolved.
		
		:return: A Silico_options object (a fancy dict).
		"""
		# First, we always load the default config object.
		configs = [self.default_config()]
		
		# Next the 'master' (located in the silico install directory).
		configs.extend(self.list_from_location(self.MASTER_CONFIG_LOCATION()))
		# Next the system (located at /etc/silico)
		configs.extend(self.list_from_location(self.SYSTEM_CONFIG_LOCATION()))
		# Finally the user specific (located at ~/.config/silico)
		configs.extend(self.list_from_location(self.USER_CONFIG_LOCATION()))
		
		# And return a merged version.
		return Silico_options(configs)
	
		
	def load(self, config_file):
		"""
		Load the config options from the config file that this class represents.
		
		This list will first be emptied before being updated with the config file's contents.
		
		:return: Nothing.
		"""
		# Clear our self first (necessary?)
		self.clear()
		
		# Read our YAML file.
		data = yaml.safe_load_all(config_file)
		# Convert empty/None configs to empty dicts.
		self.extend([self.pre_process(datum) if datum is not None else {} for datum in data])
	
	def pre_process(self, config):
		"""
		Pre-process a config file immediately after it has been read.
		"""
		# All done.
		return Config(config, FILE_NAME = self.file_path) if 'TYPE' not in config else Configurable(config, FILE_NAME = self.file_path)
				
	@classmethod
	def default_config(self):
		"""
		Get a Config_file_loader object that contains configs in the default config file.
		
		:return: Config_file_loader object.
		"""
		# Get a full path to our resource.
		# And get the object
		return self.from_file(Path(self.DEFAULT_CONFIG_LOCATION(), self.DEFAULT_CONFIG_NAME))
		
	@classmethod
	def list_from_location(self, location):
		"""
		Load all configs found under a given location.
		
		:param location: The path to a directory to load from. In this directory, the files 'silico.yaml' and 'conf.d/**/*.yaml' will be loaded.
		"""
		# The names of files to load.
		config_files = [Path(location, self.MASTER_CONFIG_NAME)]
		
		# And add the names of any files found under conf.d
		confd_files = glob.glob(str(Path(location, self.MULTI_DIRECTORY_NAME, "**", "*.yaml")), recursive = True)
		# Now sort them by filename so we load in the correct order.
		confd_files.sort(key= lambda file_path: Path(file_path).name)
		# Add add after the master file (so we have greater precedence).
		config_files.extend(confd_files)
		
		# Now load.
		config_objects = []
		for config_file in config_files:
			try:
				config_objects.append(self.from_file(config_file))
			except FileNotFoundError:
				# This is ok.
				pass

				
		# All done.
		return config_objects
	
	