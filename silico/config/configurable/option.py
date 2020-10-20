from silico.exception.configurable import Configurable_exception

class Option():
	"""
	Class for specifying an option in a configurable.
	
	Options are descriptors that perform type checking and other functionality for Configurables; they expose the options that a certain configurable expects.
	"""
	
	def __init__(self, name = None, *, default = None, help = None, choices = None, validate = None, type = None, rawtype = None, exclude = None, required = False, configure = None, no_edit = False):
		"""
		Constructor for Configurable Option objects.
		
		:param name: The name of this option. If None is given this will be determined automatically from the name of the attribute this option is stored under.
		:param default: Default value for this option. Alternatively, default can be a callable which will be called with 2 arguments: this Option object and the owning Configurable object and should return the default value.
		:param help: Descriptive help string
		:param choices: An iterable of valid choices for this option. Alternatively, choices can be a callable which will be called with 2 arguments: this Option object and the owning Configurable object and should return the list of options.
		:param validate: Function called to check that the given value is valid. The function will be called with 3 arguments: this Option object, the owning Configurable object and the value being set, and should return True or False as appropriate.
		:param type: A callable that is used to set the type of value.
		:param exclude: A list of strings of the names of attributes that this option is mutually exclusive with.
		:param required: Whether this option is required or not.
		:param configure: Callable that will be called with 3 arguments: this Option object, the owning Configurable object and the current value after type conversion and should return the configured value.
		:param no_edit: Flag to indicate that this option shouldn't be edited.
		"""
		self.name = name
		self._default = default
		self.type = type
		self.rawtype = type if rawtype is None else rawtype
		self.help = help
		self._choices = choices
		self.validate = validate if validate is not None else lambda option, configurable, value: True
		self.exclude = exclude if exclude is not None else []
		self.required = required
		self._configure = configure if configure is not None else self.default_configure
		self.no_edit = no_edit
		
	@classmethod
	def default_configure(self, option, configurable, value):
		"""
		The default configure method, does nothing.
		"""
		return value
		
	def choices(self, obj):
		"""
		Get the list of allowed options for this option.
		
		This property will evaluate self._choices if it is a callable.
		
		:return: The list of choices, or None if no choices.
		"""
		if not callable(self._choices):
			return self._choices
		else:
			return self._choices(self, obj)
		
	def __set_name__(self, cls, name):
		"""
		Called automatically during class creation, allows us to know the attribute name we are stored under.
		"""
		self.name = name if self.name is None else self.name
	
	def default(self, obj):
		"""
		Get the default value of this Option.
		"""
		if not callable(self._default):
			return self._default
		else:
			return self._default(self, obj)
		
	def is_default(self, obj):
		"""
		Whether the value of this option is currently the default or not.
		"""
		return not self.name in obj
	
# 	def set_default(self, obj):
# 		"""
# 		Reset this option to default.
# 		"""
# 		obj.pop(self.name)
	
	def getraw(self, obj, dictobj = None):
		"""
		Get the 'raw' value of this Option.
		
		The raw value is that which would appear in the corresponding config file; option.getraw(obj) is roughly equivalent to obj[option.name].
		"""
		dictobj = dictobj if dictobj is not None else obj
		
		if not self.required:
			return dictobj.get(self.name, self.default(obj))
		else:
			try:
				return dictobj[self.name]
			except KeyError:
				raise Configurable_exception(obj, "missing required option '{}'".format(self.name))
			
	def setraw(self, obj, value, dictobj = None):
		"""
		Set the 'raw' value of this Option.
		"""
		dictobj = dictobj if dictobj is not None else obj
		
		dictobj[self.name] = value
		
		
	def __get__(self, obj, cls):
		"""
		Get the 'real' value of this Option.
		"""
		if obj is None:
			return self
		# TODO: There is lots of room for improvement here; the whole Configurables inheriting from dicts can probably be scrapped...
		try:			
			return obj._configurable_options[self.name]
		except KeyError:
			# We haven't been configured yet.
			#return obj.get(self.name, self.default)
			return self.getraw(obj)
		
	def __set__(self, obj, value, dictobj = None):
		"""
		Set the value of this option.
		
		This method will update the underlying raw value (so changes can be preserved) and will also call configure(); validating the value given.
		"""
		old = self.getraw(obj, dictobj = dictobj)
		old_def = self.is_default(obj)
		
		self.setraw(obj, value, dictobj = dictobj)
		try:
			self.configure(obj)
		except Exception:
			if old_def:
				# Pretty sure this is wrong.. (should be dictobj?)
				obj.pop(self.name)
			else:
				self.setraw(obj, old, dictobj = dictobj)
			raise
		
	def __delete__(self, obj):
		"""
		Delete the explicit value of this option, resorting to the default (if one is given).
		
		If this option is required; attempting to delete it will raise a Configurable_exception().
		"""
		del(obj[self.name])
		self.configure(obj)
		
			
	def pre_configure(self, obj, dictobj = None):
		"""
		Set the value of this Option.
		
		:param obj: The Configurable object that is being set.
		:param dictobj: The dictionary object from where the raw value is being retrieved. In most cases this is the same as obj (and this is the default).
		"""
		dictobj = dictobj if dictobj is not None else obj
		value = self.getraw(dictobj)
		
		# Try and cast to our type (if not default).
		#if not self.is_default(dictobj):
		if not self.is_default(dictobj) and value is not None:
			try:
				value = self.type(value) if self.type is not None else value
			except (TypeError, ValueError):
				raise Configurable_exception(obj, "value '{}' of type '{}' is invalid for configurable option '{}'".format(value, type(value).__name__, self.name))
		
		# If we have a list of options, check we chose one.
		choices = self.choices(obj)
		
		if choices is not None:
			# If we a are a list type, we'll check each item in value (rather than value itself).
			try:
				if issubclass(self.type, list) or issubclass(self.type, tuple):
					values = value
				else:
					values = [value]
			except TypeError:
				# issubclass raises this all the time...
				values = [value]
			
			for subvalue in values:
				if subvalue not in choices:
					raise Configurable_exception(obj, "value '{}' is not one of the allowed choices for option '{}'".format(subvalue, self.name))
			
		# Check the value is valid.
		if not self.validate(self, obj, value):
			# Invalid.
			raise Configurable_exception(obj, "value '{}' of type '{}' is invalid for configurable option '{}'".format(value, type(value).__name__, self.name))
		
		return self._configure(self, obj, value)
	
	def configure(self, obj):
		"""
		Configure this Option.
		"""
		# All good.
		obj._configurable_options[self.name] = self.pre_configure(obj)
		