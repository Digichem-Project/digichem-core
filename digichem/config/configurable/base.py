from copy import deepcopy, copy
import itertools

from silico.config.base import Config
from silico.exception import Configurable_exception, Missing_option_exception
from silico.misc import Dynamic_parent
from silico.exception.base import Silico_exception
from silico.config.configurable.option import Option

class Configurable_list(list):
	"""
	An augmented list the that contains methods for manipulating Configurable objects.
	"""
	
	def filter_hidden(self):
		"""
		Filter out configurables that are hidden, returning a new Configurable_list of concrete configurables only.
		"""
		return type(self)([conf for conf in self if not conf.ABSTRACT])
		
	
	def search(self, *identifiers, max_index = None):
		"""
		Search for Configurable objects that match a given ID.
		
		:param identifier: Either the NAME, ALIAS or index of the configs to get.
		:param max_index: Only return matches that have an index less than this value in configs.
		:return: A Configurable_list of Configurables that match identifier.
		"""
		found = type(self)()
		for identifier in identifiers:
			# Get alllll the targets that match.
			for index, config in [(index, config) for index, config in enumerate(self) if config.match(identifier)]:
				# Only add to found if not already in there, and if index is less than max_index.
				if config not in found and (max_index is None or index < max_index):
					found.append(config)
		
		# Return the match.
		return found
	
	def get_config(self, identifier):
		"""
		Get a Configurable object that matches a given ID.
		
		:raises Configurable_class_exception: If the Configurable cannot be found (or more than one was found).
		:param identifier: Either the NAME, ALIAS or index of the config to get.
		:param configs: List of Configurables to search through.
		"""
		if isinstance(identifier, int) or identifier.isdigit():
			identifier = int(identifier) -1
			try:
				# Don't allow negative indices.
				if identifier < 0:
					raise IndexError("ID out of range (must be 1 or greater)")
				
				return self[identifier]
			except IndexError:
				raise Silico_exception("could not find Configurable with index '{}'".format(identifier +1))
		
		# Get matches
		match = self.search(identifier)		
		
		# Get upset if we got nothing.
		if len(match) == 0:
			raise Silico_exception("could not find Configurable with NAME/ALIAS '{}'".format(identifier))
		
		# Return the match.
		return match[0]
	
	def resolve(self):
		"""
		Resolve this list of Configurables.
		
		Resolving involves merging those Configurables that inherit from each other.
		Note that resolve() partially modifies the list IN PLACE (deepcopy your list first if you want to preserve it), but the final returned list is a new object.

		:return: A modified Configurable_list.
		"""
		# First, resolve duplicates.
		self.resolve_duplicates() 
		
		# Now, resolve parents (do this in place so we're not constantly resolving the same parents).
		for index, config in enumerate(self):
			self[index] = config.resolve_parents(self)
			
		# Resolve splits and return (but only those that are not abstract).
		resolved = type(self)(resolved for resolved in self.resolve_splits() if not resolved.ABSTRACT)
		
		# Finally, sort out names.
		for config in resolved:
			config.resolve_names()
		
		# All done.
		return resolved
	
	def resolve_splits(self):
		"""
		"""
		return type(self)([resolved for config in self for resolved in config.resolve_split()])
		
	
	def configure(self, **kwargs):
		"""
		Call the second init function on all Configurables contained within this list.
		
		:param **kwargs: Key-word arguments that will be passed to configure() of each config in this list.
		"""
		# Set ID's on our children.
		index = 0
		while index < len(self):
			config = self[index]
			# First get the class that the config asks for.
			try:
				typecls = config.from_class_handle(config.TYPE)
				cls = typecls.from_class_handle(config.CLASS)
			except KeyError:
				raise Missing_option_exception(config, 'CLASS')
			except Exception:
				raise Configurable_exception(config, "failed to load class")
			
			# Now init (sets key:values).
			config = cls(config, FILE_NAME = config.FILE_NAME)
			
			# Init properly.
			config.configure(ID = index +1, CONFIGS = self, **kwargs)
# 			try:
# 				#config.post_init(ID = index +1, CONFIGS = self, **kwargs)
# 				config.configure(ID = index +1, CONFIGS = self, **kwargs)
# 			except TypeError:
# 				raise Configurable_exception(config, "unrecognised or missing key")
			
			# Update the list.
			self[index] = config
			index +=1
	
	def resolve_duplicates(self):
		"""
		Resolve duplicate Configurables from this list.
		
		Configurables will overwrite (merge with) other Configurables with the same name/alias that are earlier in the list.
		Resolving will also remove these later duplicates so only one, merged Configurable will remain in the list (at the position one first appeared).
		
		:return: The (possibly shorter) resolved list.
		"""
		# Current position in the list.
		index = 0
		while index < len(self):
			config = self[index]
			# Search for configs with the same name that are earlier in the list.
			duplicates = self.search(*config.NAMES, max_index = index)

			# If there is more than one duplicate, get upset.
			# We can't support multiple duplicates because we replace the duplicate at its current position in the list.
			if len(duplicates) > 1:
				raise Configurable_exception(config, "cannot overwrite multiple duplicates '{}'".format(", ".join([duplicate.description for duplicate in duplicates])))
			elif len(duplicates) == 1:
				# We have a duplicate.
				duplicate = duplicates[0]
				# Merge the lists of PARENT and ALIAS (because lists are not normally merged). This might be unexpected?
				config.PARENT.extend([parent for parent in duplicate.PARENT if parent not in config.PARENT])
				config.ALIAS.extend([alias for alias in duplicate.ALIAS if alias not in config.ALIAS])
				
				# Merge the duplicate with the current config.
				# Because we are not replacing duplicate (it's the same object), this will update the config in the list.
				duplicates[0].merge(config)
				
				# Now we remove config from the list.
				self.pop(index)
				
				# We don't need to increment index (because index already points to the next config and len(self) has decreased).
			else:
				# No duplicates.
				index += 1
		
		# All done.
		return self

class Options_mixin():
	"""
	Mixin class for those that contain configurable options.
	"""
	
	@property
	def OPTIONS(self):
		"""
		Get a list of all Configurable Options of this object.
		"""
		return {getattr(type(self), attr).name: getattr(type(self), attr) for attr in dir(type(self)) if isinstance(getattr(type(self), attr), Option)}

class Configurable(Config, Dynamic_parent, Options_mixin):
	"""
	Class that represents a Configurable.
	
	Configurables are a specific type of Config objects that are used to instantiate a class. Like Config objects, they are enhanced dicts.
	"""
	
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self._configurable_options = {}
		
	# Configurable options.
	_NAME = Option("NAME", help = "The unique name of this Configurable and the name of the folder under which this calculation will run. If left blank and a GROUP_NAME is set, a name will be generated automatically", type = str)
	ALIAS = Option(help = "A list of alternative names for this Configurable", default = [], type = list, no_edit = True)
	GROUP = Option(help = "An ordered list of group names, used for categorisation", default = [], type = list)
	GROUP_NAME = Option(help = "An alternative to NAME; GROUP_NAME is used to set a shorter name that appears as part of a GROUP", default = None, type = str)
	SPLIT = Option(help = "A list of configurable splits; an advanced feature for automatically creating Configurables...", default = [], type = list, no_edit = True)
	PARENT = Option(help = "An ordered list of parent Configurables from which this Configurable will inherit from", default = [], type = list, no_edit = True)
	_ABSTRACT = Option("ABSTRACT", help = "A flag to force a Configurable to be abstract. Abstract Configurables cannot be selected, but can be used as parents for other Configurables", default = False, type = bool, no_edit = True)
	TYPE = Option(help = "The type of this Configurable", choices = ('method', 'program', 'calculation', 'basis_set'), required = True, default = 'calculation', type = str, no_edit = True)
	CLASS = Option(
		help = "The specific sub-type of this Configurable",
		choices = lambda option, configurable: [handle for cls in Configurable.from_class_handle(configurable.TYPE).recursive_subclasses() if hasattr(cls, 'CLASS_HANDLE') for handle in cls.CLASS_HANDLE],
		required = True,
		no_edit = True
	)
		
	def configure(self, *, ID, CONFIGS = None):
		"""
		Configure this configurable.
		
		:param ID: The unique ID of this configurable.
		:param CONFIGS: List of other Configurables (constructed but not necessarily configured).
		"""
		# Save our ID.
		self.ID = ID
		
		# First, get a list of our configurable options.
		options = self.OPTIONS
		
		# Now configure each.
		for name,option in options.items():
			option.configure(self)
		
		# Check exclusions.
		for name,option in options.items():
			for exclusion in option.exclude:
				# TODO: This might need work...
				if not options[exclusion].is_default(self) and not option.is_default(self):
					raise Configurable_exception(self, "options '{}' and '{}' cannot be set at the same time (mutually exclusive)".format(option.name, exclusion))
	
	def resolve_names(self):
		"""
		Resolve names (which can contain formatting characters).
		
		Currently, only NAME and GROUP_NAME are resolved, and resolution only occurs for SPLIT configs.
		This is because SPLIT occurs last (after resolving PARENT, for example), so can't be used for inheritance (and hence resolving ALIAS would be pointless). 
		"""
		if self.get('NAME') is not None:
			self['NAME'] = self['NAME'].format(**self)
		self['GROUP'] = [group.format(**self) for group in self.GROUP]
		if self.get('GROUP_NAME') is not None:
			self['GROUP_NAME'] = self['GROUP_NAME'].format(**self)
		#self['ALIAS'] = [alias.format(**self) for alias in self.ALIAS]
		#self.pop('ALIAS')
			
	
	def match(self, identifier):
		"""
		Determine whether a given identifier matches this Configurable object.
		
		Identifiers match (case insensitively) to either NAME or any of ALIAS.
		"""
		return identifier.lower() in [name.lower() for name in self.NAMES]
						
	
	def resolve_parents(self, configs):
		"""
		Resolve the inheritance of a Configurable, returning a new Configurable that is merged with its parents (if any).
		
		:param configs: List of other Configurables to inherit from.
		:return: A  NEW resolved Configurable.
		"""		
		# Each parent specified will overwrite the options set by the previous parent, and the config itself will overwrite the last parent.
		# Loop through each parent first
		merged_parents = Configurable(FILE_NAME = self.FILE_NAME)
		for parent_name in copy(self.PARENT):
			# If we got this far then we have some work to do.
			# Fetch the Configurable that is referenced.
			try:
				#parent = self.from_ID(parent_name, configs)
				parent = configs.get_config(parent_name)
			except Silico_exception as e:
				# Couldn't find parent.
				raise Configurable_exception(self, "could not find parent") from e		
			
			# Clone.
			parent = deepcopy(parent)
			
			# Also resolve this parent config.
			parent = parent.resolve_parents(configs)
			
			# Remove the parent's name and alias (because we don't want to inherit these).
			parent.pop('NAME', None)
			parent.pop('ALIAS', None)
			parent.pop('PARENT', None)
			parent.pop('ABSTRACT', None)
			merged_parents = merged_parents.merge(parent)
		
		# Now we merge our combined parent with our real config.
		merged_parents.merge(deepcopy(self))
		
		# Finally, delete all PARENTs so we don't resolve this config again.
		merged_parents.pop('PARENT', None)
		
		# All done.
		return merged_parents


	def resolve_split(self):
		"""
		Recursively resolve any SPLITs of this configurable.
		
		SPLIT is used to quickly define children of a Configurable; it is a type of automatic inheritance that can work in conjunction with PARENT.
		SPLIT is a list of dicts; each key in each dict is also a list. The key in each dict specifies fields to split on, the list members are the split values.
		
		:return: A list of Configurables.
		"""
		# If we don't have any splits then we can simply return ourself.
		if len(self.SPLIT) == 0:
			return [self]
		
		
		resolved_splits = Configurable_list()
		
		# Iterate through each SPLIT.
		for split in self.SPLIT:
			# Certain fields we don't split on (because it would either be too complicated or nonsensical.
			split = copy(split)
			# Remove these protected fields.
			protected = {
				'SPLIT': split.pop('SPLIT', [])
			}
			for protected_key in ('NAME', 'GROUP', 'GROUP_NAME'):
				if protected_key in split:
					protected[protected_key] = split.pop(protected_key)
			
			# Next, expand each key (the value of each key is a list, we need a list of one membered dicts instead).
			expanded = ([{key: value} for value in lst] for key, lst in split.items())
			
			# Now we need every combination, thanks itertools!
			# Combinations is list of tuples where each value is a key: value pair that we've been asked to split on.
			combinations = itertools.product(*expanded)
			
			# Now rebuild back into a list of dicts.
			resolved_splits.extend([self.merge_dict(protected, {key: value for pair in combination for key, value in pair.items()}) for combination in combinations])
			
		# Now we have (possibly very big) list of Configurables, but they only contain those values we were asked to split on.
		# Next we need to merge with ourself to fill out other values.
		merged_list = Configurable_list()
		for raw in resolved_splits:
			# First, copy.
			merged = deepcopy(self)
			
			# Delete any alias.
			merged.pop('ALIAS', None)
			
			# We append to GROUP rather than overwriting.
			group = list(self.GROUP)
			group.extend(raw['GROUP'] if 'GROUP' in raw else [])
			raw['GROUP'] = group
						
			# Merge
			merged.merge(raw)
			
			# Finally, resolve NAME and GROUP_NAME.
			#merged.resolve_names()
			
			# Add to list.
			merged_list.append(merged)
						
		# All done, return the list.
		return merged_list.resolve_splits()
			
	@property
	def ABSTRACT(self):
		"""
		Whether or not this Configurable is 'hidden'.
		
		Hidden configurables cannot be selected directly, but they can be referenced by other configurables (vaguely analogous to abstract classes).
		"""
		return self['ABSTRACT'] if 'ABSTRACT' in self else self.NAME is None
		
	@property
	def NAME(self):
		"""
		The unique name by which this config is known.
		
		Note that 'hidden' configurables don't have a NAME by definition, but they will probably have a number of ALIAS set instead.
		
		:return: The NAME (a string) or None if this configurable is hidden.
		"""
		# If we have a real NAME set, we use that.
		if self.get("NAME", None) is not None:
			return self['NAME']
		elif len(self.GROUP) > 0 and self.GROUP_NAME is not None:
		# If we are part of a group, then we might be able to make a name.
			return (" ".join(self.GROUP)) + " " + self.GROUP_NAME
		else:
			# No names we can use.
			return None
			
	@property
	def NAMES(self):
		"""
		A list of all the unique names this Configurable is known by (NAME + ALIAS).
		"""
		names = []
		
		# Add NAME if we have one.
		if self.NAME is not None:
			names.append(self.NAME)
		
		# Add any ALIAS.
		names.extend(self.ALIAS)
		
		# Done
		return names
	
	@property
	def description(self):
		"""
		A string that describes this Configurable object.
		"""
		# Start with the name.
		desc = self.NAME if not self.ABSTRACT else "ABSTRACT"
		
		# Add alias if present.
		if len(self.ALIAS) > 0:
			desc += " ({})".format(" ,".join(self.ALIAS))
			
		return desc
					
	def grouped(self, grouped_dict, *, names = None, number):
		names = self.GROUP if names is None else names
		
		# This is the name we will be adding to in grouped_dict.
		# If it is None, we will be adding to a list (we have no more names to go through).
		curname = names[0] if len(names) > 0 else None
			
		if curname is None:
			# Create a new list if one doesn't exist.
			if curname not in grouped_dict:
				grouped_dict[curname] = []
			
			# Add to list
			grouped_dict[curname].append((number, self))
		else:
			# We have more names to work through.
			# Create an empty dict under the name if not existing.
			if curname not in grouped_dict:
				grouped_dict[curname] = {}
				
			# Recurse.
			self.group(grouped_dict[curname], names[1:])
			
		return grouped_dict
			
		
		