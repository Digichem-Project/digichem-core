from pathlib import Path
import yaml
import pkg_resources
from logging import getLogger
import logging
from silico.exception.base import Config_loader_exception, Bad_config_exception
import glob
from itertools import filterfalse
from copy import deepcopy
import silico
import itertools

class Config(dict):
	"""
	Class for holding config options from various sources.
	"""
	
	def __init__(self, configs = None, load_defaults = True):
		"""
		Constructor for Config objects.
		"""
		super().__init__()
		
		# Set our default properly.
		if configs is None:
			configs = []
		
		# We keep a list of all our config sources as well as merging them into ourself. Might be useful if you wanted to compare different options.
		self.configs = []
		
		# First, we always load the default config object.
		if load_defaults:
			self.add_config(Config_file_loader.default_config())
		
		# Add our list.
		for config in configs:
			self.add_config(config)
			
	def __getitem__(self, key):
		# TODO: This only works for the top-level dict object for now, we need to replace other dicts with Config objects.
		try:
			return super().__getitem__(key)
		except KeyError:
			raise KeyError("Unknown config option '{}'".format(key))

	def add_config(self, config):
		"""
		Add a set of options to this config object.
		
		The functioning of this method is controlled by the 'CONFIG_TYPE' key in config. If this key is None (or is missing entirely), then config is taken to contain normal config options which are added to this object so that later options will override earlier ones.
		If 'CONFIG_TYPE' is not 'None' or missing, then 'CONFIG_TYPE' specifies the name of a key in this Config object to which the given config is appended as a whole object.
		
		:param config: The new config to add; a dict of options. For convenience, a list of dicts can also be given, and they will be added in order. Further, higher-order lists (lists of lists of dicts etc) can also be given and will be fully traversed and added in order. 
		:return: self, for convenience.
		"""
		# If config is a list, we'll recursively call ourself.
		if isinstance(config, list):
			for sub_config in config:
				self.add_config(sub_config)
			# Done.
			return
		
		# Add to our list.
		self.configs.append(config)
		
		# Either merge or append, depending on what sort of config we have.
		if config.get('CONFIG_TYPE', None) is None:
			# Normal config, merge with ourself.
			self.merge_configs(config, self)
		else:
			# 'Multi' config
			try:
				self[config['CONFIG_TYPE']].append(config)
			except KeyError:
				# Key has not yet been set because this is the first time we've seen it.
				self[config['CONFIG_TYPE']] = [config]
			except AttributeError:
				# key has already been set and is not a list.
				raise Config_loader_exception("multi config", "key '{}' has already been set to a non-list type and cannot be extended with multi-config".format(config['CONFIG_TYPE']))
		
		# Give our self back.
		return self
	
	def set_log_level(self, logger_name):
		"""
		Set the logging level of a logger based on the config options in the object.
		
		:param logger_name: The name of a logger to set (this name will be passed to logging.getLogger().
		"""
		# Get our logger.
		logger = logging.getLogger(logger_name)
		
		# Set from log_level first.
		if self['logging']['log_level'] == "OFF":
			logger.setLevel(60)
		else:
			logger.setLevel(self['logging']['log_level'])
		
		# Now adjust with verbosity.
		verbose = self['logging']['verbose']
		if verbose is not None:
			# Set from verbosity.
			new_level = logger.level - verbose * 10
			
			# Don't allow us to reach 0 (because this is actually 'UNSET').
			if new_level <= 0:
				new_level = 10
			
			# And set.
			logger.setLevel(new_level)
	
	@classmethod
	def resolve_config_inheritance(self, config, configs):
		"""
		Resolve the inheritance of a config dict, returning a new dict that is merged with it's parents (if any).
		
		:param config: The config dict to resolve.
		:param configs: List of other config dicts. This will be searched to find dicts to merge with config.
		:return: A resolved config dict.
		"""
		# First get the parent option from our config.
		parent_names = config.get('CONFIG_PARENT', [])
		
		# Loop through each parent.
		for parent_name in parent_names:
			# If we got this far then we have some work to do.
			# First get the references config.
			try:
				parent = self.search_configs(parent_name, configs)
			except Bad_config_exception:
				raise Bad_config_exception("Could not apply config parent '{}' of config '{}'".format(parent_name, config.get('CONFIG_NAME', "(no name)")))
			
			# Clone.
			parent = deepcopy(parent)
			
			# Also resolve this parent config.
			parent = self.resolve_config_inheritance(parent, configs)
			
			# Now combine.
			config = deepcopy(config)
			
			# Remove the parent's name and alias (because we don't want to inherit these).
			parent['CONFIG_NAME'] = None
			parent['CONFIG_ALIAS'] = []
			parent['CONFIG_PARENT'] = []
			config = self.merge_configs(config, parent, none_to_old = False)
		
		# Give it back.
		return config
			
		
	@classmethod
	def search_configs(self, config_name, configs):
		"""
		Get a config dict from a list of configs by its unique 'CONFIG_NAME' or 'CONFIG_ALIAS'.
		
		:raises Bad_config_exception: If the config cannot be found (or more than one was found).
		:param config_name: The name/alias of the config to get (both are searched simultaneously).
		:param configs: List of config dicts.
		"""
		# Get alllll the targets that match (either by name or alias)
		matching_configs = list(filterfalse(lambda conf: conf.get('CONFIG_NAME', None) != config_name and config_name not in conf.get('CONFIG_ALIAS', []), configs) )
		
		# If there's more than one match, the same name/alias has been used more than once.
		if len(matching_configs) > 1:
			raise Bad_config_exception("The name/alias '{}' has been used by more than one config object, name/alias must be unique".format(config_name))
		elif len(matching_configs) == 1:
			return matching_configs[0]
		else:
			# No match.
			raise Bad_config_exception("Unable to locate config with name/alias '{}'".format(config_name))
		
	@classmethod
	def merge_configs(self, new, old, *, none_to_old = True):
		"""
		Recursively merge two dictionaries of options.
		
		You are probably looking for add_config() (unless you are modifying this implementation).
		
		Any options specified in new will overwrite those specified in old.
		
		:param new: A new dictionary to merge into the old.
		:param old: An old dictionary to be overwritten by new.
		:param none_to_old: If True (the default), values of None in the new dictionary will be ignored (and so the value from old will be used). If False, None will be set in the returned dictionary.
		"""
		# This is inspired by https://stackoverflow.com/questions/823196/yaml-merge-in-python.
		# Only attempt to merge if both are dicts (otherwise we just return new straight away).
		if isinstance(new, dict) and isinstance(old, dict):
			for key, value in new.items():
				if key not in old:
					# The option is new, so we can just set without danger.
					old[key] = value
				else:
					# Both old and new have the same key, so we need to merge.
					old[key] = self.merge_configs(value, old[key], none_to_old = none_to_old)
			return old
		# This removed behaviour was intended for the main config file (it allows you to explicitly reset an option to the one in the default.yaml file) but it clashes with the submit config files (where setting None is designed to unset an option set by a CONFIG_PARENT)...
		elif new is None and none_to_old:
			return old
		else:
			return new
		
	@classmethod
	def standard_configs(self):
		"""
		Load the standard config options (from the master and user specific files).
		"""
		configs = []
		# First the 'master' (located in the silico install directory).
		configs.extend(Config_file_loader.list_from_location(Config_file_loader.MASTER_CONFIG_LOCATION()))
		# Next the system (located at /etc/silico)
		configs.extend(Config_file_loader.list_from_location(Config_file_loader.SYSTEM_CONFIG_LOCATION()))
		# Finally the user specific (located at ~/.config/silico)
		configs.extend(Config_file_loader.list_from_location(Config_file_loader.USER_CONFIG_LOCATION()))
		
		return self(configs)
	
	@classmethod
	def get_standard_args(self, argparser):
		"""
		Add standard silico arguments to an argparse ArgumentParser object.
		
		The standard arguments are added as a group called 'general options'.
		
		:param argparser: An argparse ArgumentParser object
		:return: The group is returned for convenience.
		"""
		group = argparser.add_argument_group("general options", "general options that control various aspects of silico")
		group.add_argument("-V", "--verbose", help = "increase verbosity, stack with itself to further increase verbosity (this option overrides log_level)",  action='count', default = None)
		group.add_argument("--log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'OFF'], help = "the level of messages to print", default = None)
		group.add_argument("-K", "--alignment", choices=['MIN', 'FAP', 'AA', 'AAA'], help = "the alignment method to use to align atoms and calculate geometry data. Options are MIN: Minimal, FAP: Furthest Atom Pair, also known as Kebab, AA: Average Angle, also known as Kebab+, AAA: Adjusted Average Angle", default = None)
		group.add_argument("-A", "--angle_units", help = "the units to use to print angles. Options are deg: degrees or rad: radians", default = None)
		group.add_argument("-S", "--setting", help = "set a config option to a value. Options of this type are parsed as if they were a config file (in yaml format) and are then used to set corresponding options, eg -S \"absorption_graph: {fwhm: 100}\"", nargs = "*", default = [], action = "append")
		group.add_argument("--config_files", help = "an additional config file to read from. See the master config file for possible config options. Note that the master config file (at silico/data/config/silico.yaml) and user config file (at ~/.config/silico/silico.yaml) are always read automatically and do not need to be specified here. Multiple files may be given and will be processed in the order specified (the last having highest priority)", nargs = "*", default = [])
		return group
	
	def add_standard_args(self, args):
		"""
		Add the standard command line arguments to this config object.
		
		:param args: An argparse Namespace object.
		"""
		# Build a dictionary of options as they would appear in the config file.
		arg_configs = {
			'logging': {
					'log_level': args.log_level,
					'verbose': args.verbose
				},
			'alignment': args.alignment,
			'angle_units': args.angle_units
			}
		
		# Now merge this new dictionary with ourself.
		self.add_config(arg_configs)
		
		# We will also process any additional config files given on the command line, and add them in order.
		for config_file_name in args.config_files:
			# Load the file.
			config = Config_file_loader.from_file(config_file_name)
			# And add to ourself.
			self.add_config(config)
			
		# Finally, we set any config options.
		# We do this last so that they'll have highest priority.
		# We read each individually because otherwise options that have been set twice won't merge properly (the yaml parser will overwrite the older).
		for config_file in itertools.chain(*args.setting):
			# Now we'll parse it as YAML.
			try:
				#self.add_config(list(yaml.safe_load_all(config_file)))
				self.add_config(Config_file_loader(config_file))
			except Exception:
				raise Config_loader_exception("<command_line>", "failed to parse command-line config options")
			
			
		# All done.

class Config_file_loader(list):
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
	
	def __init__(self, config_file):
		"""
		Constructor for config file loaders.
		
		:param config_file: A loaded config file as a string (or an open file object to read from).
		"""
		super().__init__()
		
		# Load our options if we've been asked to.
		if config_file is not None:
			self.load(config_file)
			
	@classmethod
	def from_file(self, config_file_path):
		"""
		Alternative constructor that first opens a config file path before reading it.
		"""
		with open(config_file_path, "rt") as config_file:
			return self(config_file)
		
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
	
	@classmethod
	def pre_process(self, config):
		"""
		Pre-process a config file immediately after it has been read.
		"""
# 		# If config doesn't have a CONFIG_NAME set but does have a group, use that.
# 		if config.get('CONFIG_NAME', None) is None and config.get('group', None) is not None and config.get('group_sub_name', None) is not None:
# 			config['CONFIG_NAME'] = config['group'] + " " + config['group_sub_name']
		
		# All done.
		return Configurable(config)
				
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
				getLogger(silico.logger_name).debug("Not loading config file '{}'' file not found".format(config_file))
				
		# All done.
		return config_objects
	
class Configurable(dict):
	"""
	Class that represents a configurable dict.
	
	This class is very WIP...
	"""
	
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		
		self._set_CONFIG_NAME()
	
	def __setitem__(self, key, value):
		super().__setitem__(key, value)
		
		self._set_CONFIG_NAME()
	
	def _set_CONFIG_NAME(self):
		"""
		This method is a bit of hack and will be replaced when this class is rewritten.
		"""
		# See if we have a CONFIG_NAME set or not.
		if self.get("CONFIG_NAME", None) is None and self.get('group', None) is not None and self.get('group_sub_name', None) is not None:
			# We have no name set but do have a group name; use that instead.
			self["CONFIG_NAME"] =  self['group'] + " " + self['group_sub_name']
		
	
# 	def get(self, key, default = None):
# 		try:
# 			return self[key]
# 		except KeyError:
# 			return default
# 
# 	def __getitem__(self, key):
# 		if key == "CONFIG_NAME":
# 			# See if we have a CONFIG_NAME set or not.
# 			if super().get("CONFIG_NAME", None) is None and super().get('group', None) is not None and super().get('group_sub_name', None) is not None:
# 				# We have no name set but do have a group name; use that instead.
# 				return self['group'] + " " + self['group_sub_name']
# 		return super().__getitem__(key)


			