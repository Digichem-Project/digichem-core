from copy import deepcopy, copy

from silico.config.base import Config
from silico.exception import Configurable_exception, Missing_option_exception
from silico.exception.configurable import Configurable_class_exception

class Configurable_list(list):
	"""
	An augmented list the that contains methods for manipulating Configurable objects.
	"""
	
	
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
				if identifier <= 0:
					raise IndexError("ID out of range (must be 1 or greater)")
				
				return self[identifier]
			except IndexError:
				raise Configurable_class_exception(self, "could not find Configurable with index '{}'".format(identifier +1))
		
		# Get matches
		match = self.search(identifier)		
		
		# Get upset if we got nothing.
		if len(match) == 0:
			raise Configurable_class_exception(self, "could not find Configurable with NAME/ALIAS '{}'".format(identifier))
		
		# Return the match.
		return match[0]
	
	def resolve(self):
		"""
		Resolve this list of Configurables.
		
		Resolving involves merging those Configurables that inherit from each other.
		Note that resolve() modifies the list IN PLACE; deepcopy your list first if you do not want this behaviour.

		:return: The modified Configurable_list, for convenience.
		"""
		# First, resolve duplicates.
		self.resolve_duplicates()
		
		# Now, resolve parents (do this in place so we're not constantly resolving the same parents).
		for index, config in enumerate(self):
			self[index] = config.resolve_parents(self)
			
		# Return the list for convenience.
		return self
	
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

class Configurable(Config):
	"""
	Class that represents a Configurable.
	
	Configurables are a specific type of Config objects that are used to instantiate a class. Like Config objects, they are enhanced dicts.
	"""
	
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
		merged_parents = Configurable()
		for parent_name in copy(self.PARENT):
			# If we got this far then we have some work to do.
			# Fetch the Configurable that is referenced.
			try:
				#parent = self.from_ID(parent_name, configs)
				parent = configs.get_config(parent_name)
			except Configurable_class_exception as e:
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
			merged_parents = merged_parents.merge(parent)
		
		# Now we merge our combined parent with our real config.
		merged_parents.merge(deepcopy(self))
		
		# Finally, delete all PARENTs so we don't resolve this config again.
		merged_parents.pop('PARENT', None)
		
		# All done.
		return merged_parents

	
	@property
	def hidden(self):
		"""
		Whether or not this Configurable is 'hidden'.
		
		Hidden configurables cannot be selected directly, but they can be referenced by other configurables (vaguely analogous to abstract classes).
		"""
		return self.NAME is None
		
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
		elif self.get('GROUP', None) is not None and self.get('GROUP_NAME', None) is not None:
			# If we are part of a group, then we might be able to make a name.
			return self['GROUP'] + " " + self['GROUP_NAME']
		else:
			# No names we can use.
			return None
		
	@property
	def ALIAS(self):
		"""
		A list of unique names by which this config is known.
		ALIAS works the same as NAME, except it allows multiple names to be specified and it can be used on both hidden and visible configurables.
		
		:return: A list of aliases (strings) or an empty list if none are set.
		"""
		return self.get('ALIAS', [])
	
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
		desc = self.NAME if not self.hidden else "HIDDEN"
		
		# Add alias if present.
		if len(self.ALIAS) > 0:
			desc += " ({})".format(" ,".join(self.ALIAS))
			
		return desc
	
	@property
	def TYPE(self):
		"""
		The type of this configurable; determines the list of configurables to which this configurable is appended.
		
		:return: The config type.
		"""
		try:
			return self['TYPE']
		except KeyError:
			raise Missing_option_exception(self, "TYPE")
		
	@property
	def PARENT(self):
		"""
		The prents of this Configurable.
		
		:return: PARENT; a list of strings that this Configurable will inherit from.
		"""
		return self.get("PARENT", [])
	
	@property
	def CLASS(self):
		"""
		The class handle of a class that can be instantiated with this configurable.
		
		:return: The class handle (a string).
		"""
		return self.get()
	
	@property
	def GROUP(self):
		"""
		The GROUP this configurable belongs to.
		
		Groups exist mostly for layout purposes.
		
		:return: The group (a string) or None if not group is set.
		"""
		return self.get('GROUP', None)
	
	@property
	def GROUP_NAME(self):
		"""
		An alternative to NAME that can be used in combination with GROUP.
		
		If NAME is not set and if both a GROUP and GROUP_NAME are set, then that configurable's name will be created automatically from a combination of the two.
		:return: The group name (a string) or None if no group name is set.
		"""
		return self.get('GROUP_NAME', None)
					