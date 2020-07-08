
from silico.config.base import Config
from silico.config.configurable import Configurable_list
from silico.exception.configurable import Configurable_exception

class Silico_options(Config):
	"""
	Class for holding main Silico options from various sources.
	"""
	
	def __init__(self, configs = []):
		"""
		Constructor for Silico_options objects.
		"""
		super().__init__()
			
		# Set Configurable lists.
		self.methods = Configurable_list()
		self.programs = Configurable_list()
		self.calculations = Configurable_list()
		self.basis_sets = Configurable_list()
			
		# Add our list.
		for config in configs:
			self.add_config(config)
			
	def resolve(self):
		"""
		Resolve all Configurables which we contain.
		"""
		self.methods = self.methods.resolve()
		self.programs = self.programs.resolve()
		self.calculations = self.calculations.resolve()
		self.basis_sets = self.basis_sets.resolve()

	def add_config(self, config):
		"""
		Add a set of options to this config object.
		
		The functioning of this method is controlled by the 'TYPE' key in config. If this key is None (or is missing entirely), then config is taken to contain normal config options which are added to this object so that later options will override earlier ones.
		If 'TYPE' is not 'None' or missing, then 'TYPE' specifies the name of a key in this Config object to which the given config is appended as a whole object.
		
		:param config: The new config to add; a dict of options. For convenience, a list of dicts can also be given, and they will be added in order. Further, higher-order lists (lists of lists of dicts etc) can also be given and will be fully traversed and added in order. 
		:return: self, for convenience.
		"""
		# If config is a list, we'll recursively call ourself.
		if isinstance(config, list):
			for sub_config in config:
				self.add_config(sub_config)
			# Done.
			return
		
		# Either merge or append, depending on what sort of config we have.
		if not hasattr(config, 'TYPE'):
			# Normal config, merge with ourself.
			#self.merge_configs(config, self)
			self.merge(config, none_to_old = True)
		else:
			# Configurable, add to the correct list.
			try:
				getattr(self, config.TYPE + "s").append(config)
			except AttributeError:
				# Possibly because TYPE is invalid; see if we have an attribute called TYPE.
				if hasattr(self, config.TYPE + "s"):
					# Something else went wrong, don't handle it here.
					raise
				else:
					# Bad TYPE.
					raise Configurable_exception(config, "configurable TYPE '{}' is not recognised".format(config.TYPE)) from None
		
		# Give our self back.
		return self
	
	def set_log_level(self, logger):
		"""
		Set the logging level of a logger based on the config options in the object.
		
		:param logger: The logger to set (from logging.getLogger()).
		"""		
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
# 	
# 	
# 	@classmethod
# 	def get_standard_args(self, argparser):
# 		"""
# 		Add standard silico arguments to an argparse ArgumentParser object.
# 		
# 		The standard arguments are added as a group called 'general options'.
# 		
# 		:param argparser: An argparse ArgumentParser object
# 		:return: The group is returned for convenience.
# 		"""
# 		group = argparser.add_argument_group("general options", "general options that control various aspects of silico")
# 		group.add_argument("-V", "--verbose", help = "increase verbosity, stack with itself to further increase verbosity (this option overrides log_level)",  action='count', default = None)
# 		group.add_argument("--log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'OFF'], help = "the level of messages to print", default = None)
# 		group.add_argument("-K", "--alignment", choices=['MIN', 'FAP', 'AA', 'AAA'], help = "the alignment method to use to align atoms and calculate geometry data. Options are MIN: Minimal, FAP: Furthest Atom Pair, also known as Kebab, AA: Average Angle, also known as Kebab+, AAA: Adjusted Average Angle", default = None)
# 		group.add_argument("-A", "--angle_units", help = "the units to use to print angles. Options are deg: degrees or rad: radians", default = None)
# 		group.add_argument("-S", "--setting", help = "set a config option to a value. Options of this type are parsed as if they were a config file (in yaml format) and are then used to set corresponding options, eg -S \"absorption_graph: {fwhm: 100}\"", nargs = "*", default = [], action = "append")
# 		group.add_argument("--config_files", help = "an additional config file to read from. See the master config file for possible config options. Note that the master config file (at silico/data/config/silico.yaml) and user config file (at ~/.config/silico/silico.yaml) are always read automatically and do not need to be specified here. Multiple files may be given and will be processed in the order specified (the last having highest priority)", nargs = "*", default = [])
# 		return group
# 	
# 	def add_standard_args(self, args):
# 		"""
# 		Add the standard command line arguments to this config object.
# 		
# 		:param args: An argparse Namespace object.
# 		"""
# 		# Build a dictionary of options as they would appear in the config file.
# 		arg_configs = {
# 			'logging': {
# 					'log_level': args.log_level,
# 					'verbose': args.verbose
# 				},
# 			'alignment': args.alignment,
# 			'angle_units': args.angle_units
# 			}
# 		
# 		# Now merge this new dictionary with ourself.
# 		self.add_config(arg_configs)
# 		
# 		# We will also process any additional config files given on the command line, and add them in order.
# 		for config_file_name in args.config_files:
# 			# Load the file.
# 			config = silico.config.loader.Config_loader.from_file(config_file_name)
# 			# And add to ourself.
# 			self.add_config(config)
# 			
# 		# Finally, we set any config options.
# 		# We do this last so that they'll have highest priority.
# 		# We read each individually because otherwise options that have been set twice won't merge properly (the yaml parser will overwrite the older).
# 		for config_file in itertools.chain(*args.setting):
# 			# Now we'll parse it as YAML.
# 			try:
# 				#self.add_config(list(yaml.safe_load_all(config_file)))
# 				self.add_config(silico.config.loader.Config_loader(config_file))
# 			except Exception:
# 				raise Silico_exception("failed to parse command-line config options")
# 		# All done.