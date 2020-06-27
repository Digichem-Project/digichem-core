from silico.misc import Dynamic_parent
from silico.config import Config
import os
import numpy
from itertools import filterfalse, chain
from silico.exception.base import Configurable_target_exception

class Configurable_target(Dynamic_parent):
	"""
	Top level class for user-configurable submit targets (Calculation types, program types, method types etc.)
	"""
	
	def __init__(self, *, CONFIG_NAME = None, CONFIG_ALIAS = None, CONFIG_TYPE = None, CONFIG_PARENT = None, class_handle, configs):
		"""
		Constructor for submit target objects.
		
		:param CONFIG_NAME: The name of this target. The name should be unique (within the same type_name) as it can be used to refer back to this target, and should succinctly describe this target (eg, SLURM, Gaussian, Orca-DFT, Turbomole-CC, PBE0, etc.) If no name is given then this target is 'hidden'; it cannot be selected by the user but can be referenced (via alias) by other, concrete targets.
		:param CONFIG_ALIAS: A list of alternative names that the target is known by. The target can also be referenced by any of these aliases. A target having no name or alias is allowed but will be unusable.
		:param class: The specific type of this target; this should match (one of) the CLASS_HANDLE(s) of a class in calculation, method or program.
		:param configs: The list of configs from which this Configurable_target is being constructed (used by some constructors).
		"""
		# Save our args.
		self.CONFIG_NAME = CONFIG_NAME
		self.CONFIG_ALIAS = CONFIG_ALIAS if CONFIG_ALIAS is not None else []
		self.CONFIG_TYPE = CONFIG_TYPE
		self.CONFIG_PARENT = CONFIG_PARENT
		self.class_handle = class_handle
	
	@property
	def CONFIG_NAME(self):
		"""
		Get the unique name by which this config target is known (automatically converted to lower case to aid comparison; use _CONFIG_NAME for case insensitive version).
		"""
		return self._CONFIG_NAME.lower() if self._CONFIG_NAME is not None else None
	
	@CONFIG_NAME.setter
	def CONFIG_NAME(self, value):
		"""
		Set the unique name of this config target.
		"""
		self._CONFIG_NAME = value
		
	@property
	def CONFIG_ALIAS(self):
		"""
		Get a list of the unique aliases by which this config target is known (all automatically converted to lower case to aid comparison; use _CONFIG_ALIAS for case insensitive version)
		"""
		return [alias.lower() for alias in self._CONFIG_ALIAS]
	
	@CONFIG_ALIAS.setter
	def CONFIG_ALIAS(self, value):
		"""
		Set the unique aliases of this this config target.
		"""
		self._CONFIG_ALIAS = value
		
	@property
	def names(self):
		"""
		Get a list of all the names that this object can be references by (both CONFIG_NAME and CONFIG_ALIAS combined, all automatically converted to lower case to aid comparison).
		"""
		#name_list = [self.CONFIG_NAME] if self.CONFIG_NAME is not None else []
		#name_list.extend(self.CONFIG_ALIAS)
		#return [name.lower() for name in name_list]
		return list(chain([self.CONFIG_NAME], self.CONFIG_ALIAS))
	
	def get_children(self, target_list):
		"""
		Get all the configurable target objects from a list that contain this object as a parent.
		"""
		return [(index, child) for index, child in enumerate(target_list) if len(set(self.names).intersection(child.submit_parents)) > 0]
		
	@classmethod
	def filter_hidden(self, target_list):
		"""
		Remove hidden targets from a list.
		
		Hidden targets have no CONFIG_NAME set and so are not selectable by the user. They can however be inherited by other, concrete configs creating a basic inheritance chain.
		"""
		return [target for target in target_list if target.CONFIG_NAME != None]
	
	@property
	def all_possible_submit_parents(self):
		"""
		Get all the possible submit parents of this object.
		
		This default implementation raises a Configurable_target_exception. Classes for which this is a valid call should write their own implementation.
		"""
		raise Configurable_target_exception(self, "'all' is not valid for programs/methods for this target")
	
	@property
	def submit_parents(self):
		"""
		Get the list of 'Submit parents' (programs for calculations; methods for programs) that this target supports (all automatically converted to lower case to aid comparison).
		"""
		return [parent.lower() for parent in self._submit_parents]
	
	@submit_parents.setter
	def submit_parents(self, value):
		"""
		Set the list of 'Submit parents' (programs for calculations; methods for programs) that this target supports. 
		"""
		self._submit_parents = value if value != "all" else self.all_possible_submit_parents
	
	def _set_submit_parent(self, parent_type, parent):
		"""
		Convenience method for child classes, sets the 'Submit parent' (program for calculations; method for programs), first checking to make sure the submit_parent is allowed to submit this object.
		
		:raises Configurable_target_exception: If parent is not allowed.
		:param parent_type: String of the attribute name which will be set.
		:param parent: The Configurable_target object to set.
		"""
		# We check to make sure the given parent is actually allowed to be our parent.
		if parent is not None and len(set(parent.names).intersection(self.submit_parents)) == 0:
			raise Configurable_target_exception(self, "{} '{}' is not compatible with {} '{}'".format(parent.CONFIG_TYPE, parent.CONFIG_NAME, self.CONFIG_TYPE, self.CONFIG_NAME))
		# All fine, set.
		#self._program = value
		setattr(self, parent_type, parent)
	
	@classmethod
	def get_available_CPUs(self):
		"""
		Determine how many CPUs there are available to this process.
		
		It is generally not a good idea to rely on this function to determine how many CPUs to use for your calculation, not least because the environment that silico runs in is often very different to that the calculation is performed in. This function mainly exists as a fall-back.
		"""
		try:
			return len(os.sched_getaffinity(0))
			# The python docs mention that sched_getaffinity() is not available on all systems, but it is not clear what exception will be raised in such a case. These two seem most likely.
		except (AttributeError, NotImplementedError):
			return os.cpu_count()
		
	
	@classmethod
	def from_config(self, config, configs, **kwargs):
		"""
		Get a Configurable_target object from a config dict.
		
		This function automatically searches all known children of this class to find one that matches the class_handle of config.
		
		:param config: The config dict to load from.
		:param configs: A list of all configs of this type; used by some Configurable_target constructors.
		:param **kwargs: Additional keyword arguments that will be passed as-is to the new object's constructor.
 		"""
		# Doesn't make sense to call from Configurable_target itself, only its children.
		if self == Configurable_target:
			raise NotImplementedError("from_config() cannot be called from Configurable_target directly; it must be called from a child class")
		
		# First get the class that the config asks for.
		try:
			cls = self.from_class_handle(config['class_handle'])
		except KeyError:
			raise Configurable_target_exception(self, "Config with name '{}' is missing 'class_handle'".format(config.get('CONFIG_NAME')))
		
		# Now init and return.
		try:
			return cls(configs = configs, **config, **kwargs)
		except TypeError:
			raise Configurable_target_exception(self, "Unrecognised or missing key in config '{}'".format(config.get('CONFIG_NAME')))
	
	@classmethod
	def list_from_configs(self, configs, **kwargs):
		"""
		Get a list of Configurable_target objects from a list of config objects.
		"""
		# First, resolve config inheritance.
		#resolved = [Config.resolve_config_inheritance(config, configs) for config in configs]
		resolved = configs
		
		# Now remove those without a name.
		return [self.from_config(config, configs = resolved, **kwargs) for config in resolved if config.get("CONFIG_NAME", None) != None]
	
	@classmethod
	def search_list(self, identifier, target_list):
		"""
		Search a list of Configurable_target objects for one that matches a given name and/or alias.
		"""
		# TODO: This code is duplicated in Config.search_configs...
		if identifier is None:
			raise Configurable_target_exception(self, "Cannot find '{}' identifier is None".format(self.__name__))
		
		# Find which objects match out ID.
		matching_targets = list(filterfalse(lambda target: identifier.lower() not in target.names, target_list))
		
		# If there's more than one match, the same name/alias has been used more than once.
		if len(matching_targets) > 1:
			raise Configurable_target_exception(self, "The name/alias '{}' has been used by more than one config, name/alias must be unique".format(identifier))
		elif len(matching_targets) == 1:
			return matching_targets[0]
		else:
			# No match.
			raise Configurable_target_exception(self, "Unable to locate object with name/alias '{}'".format(identifier))
	
	@classmethod
	def from_name_in_list(self, identifier, target_list):
		"""
		Get a fully configured, ready to go, Configurable_target object from a given name by searching through a list of existing targets.
		"""
		if isinstance(identifier, int) or identifier.isdigit():
			identifier = int(identifier) -1
			try:
				# Don't allow negative indices.
				if identifier < 0:
					raise IndexError("ID out of range (must be 1 or greater)")
				
				return target_list[identifier]
			except IndexError:
				raise Configurable_target_exception(self, "Unknown ID '{}'".format(identifier +1))
		else:
			return self.search_list(identifier, target_list)
	
	@classmethod
	def from_name_in_configs(self, identifier, configs, **kwargs):
		"""
		Get a fully configured, ready to go, Configurable_target object from a given name by creating a list of targets from a config.
		
		:param identifier: The name of the target to get. This can match either CONFIG_NAME or one of CONFIG_ALIAS. Alternatively, if identifier is int-like, identifier can be the index of a target to get (indices start at 1).
		:param configs: A list of config dicts to get from. These configs should all be of the same type (ie, have the same CONFIG_TYPE) else strange things might occur.
		:param **init_args: Optional key-word arguments which will be passed to the constructor of the new object.
		"""
		# Get all known targets.
		#targets = self.list_from_configs(configs, **kwargs)
		
		# Now find the one we want.
		#return self.from_name_in_list(identifier, targets)
		return self.from_config(Config.search_configs(identifier, configs), configs = configs, **kwargs)
		
		
	
	
class Memory():
	"""
	Class for storing memory-amounts.
	"""
	
	# The units that we know about.
	UNITS = {
		'TB': 1000000000000,
		'GB': 1000000000,
		'MB': 1000000,
		'KB': 1000,
		 'B': 1
		}
	
	def __init__(self, value = None, print_decimal = False):
		"""
		"""
		self.value = None
		self.auto = value
		self.print_decimal = print_decimal
	
	@property
	def auto(self):
		"""
		Get the amount of memory, as a string, using an automatically determined suffix.
		"""
		# First, get an ordered list from our known units.
		ordered_units = sorted(self.UNITS.items(), key = lambda item: item[1])
		
		# Now see where our number best fits (thanks numpy).
		suffix_index = numpy.digitize((abs(self.value),), [ordered_unit[1] for ordered_unit in ordered_units])[0]
		
		# Wrap if we're out of bounds and convert to proper index.
		if suffix_index == 0:
			suffix_index += 0
		elif suffix_index == len(ordered_units):
			suffix_index -= 1
		else:
			suffix_index -= 1
		
		try:
			# Now check to see if we have a decimal part which is non zero.
			while not self.print_decimal and (self.value / ordered_units[suffix_index][1]) % 1 != 0:
				# Got a decimal part, use the next smallest unit.
				suffix_index -= 1
		except IndexError:
			# We ran out of units and couldn't remove our fraction.
			raise ValueError("Unable to represent memory value '{}'B as non decimal".format(self.value))
		
		# Now get our value in the correct units.
		value = self.value / ordered_units[suffix_index][1]
		if not self.print_decimal:
			value = int(value)
		
		return "{}{}".format(value, ordered_units[suffix_index][0])
			
		
	@auto.setter
	def auto(self, value):
		"""
		Set the amount of memory, automatically determining the unit suffix (KB, MB, GB etc).
		"""
		if isinstance(value, str):
			# First, determine the suffix (if any).
			suffix = None
			for unit, amount in  sorted(self.UNITS.items(), key = lambda item: item[1], reverse = True):
				if value.lower().endswith(unit.lower()):
					suffix = unit
					# This is a bit of a hack because all suffixes contain 'B' at the end.
					break
					
			# Now remove the suffix (if there is one), convert to float and multiply by the value of suffix (to get bytes).
			if suffix is not None:
				value = float(value[:-len(suffix)]) * self.UNITS[suffix]
		
		# Convert to int (because not sure it makes sense to represent fractions of bytes) and store.
		self.value = int(value)
		
		
	def __float__(self):
		"""
		Floatify this memory amount (returning the number of bytes).
		"""
		return float(self.value)
	
	def __int__(self):
		"""
		Intify this memory amount (returning the number of bytes).
		
		In most cases, the int and float equivalents should be the same (because fractional bytes are not supported).
		"""
		return int(self.value)
	
	def __str__(self):
		"""
		Stringify this memory amount (returning the number of bytes) with an appropriate suffix.
		"""
		return self.auto
		