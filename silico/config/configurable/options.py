from silico.config.configurable.option import Option
from silico.config.configurable.base import Options_mixin
from silico.exception.base import Silico_exception
from silico.exception.configurable import Configurable_exception


class Options(Option, dict, Options_mixin):
	"""
	A type of option that expects more options (another dict).
	"""
	# TODO: Options doesn't currently support nested options...
	
	def __init__(self, *args, name = None, help = None, default = None, exclude = None, **kwargs):
		super().__init__(name = name, help = help, rawtype = dict, default = default if default is not None else {}, exclude = exclude)
		dict.__init__(self)
		
		# Go through args and add to kwargs.
		for arg in args:
			if arg.name is None:
				raise Silico_exception("Configurable option given as positional argument to Options must have a name")
			kwargs[arg.name] = arg
		
		self._OPTIONS = {}
		
		# Set names of all kwargs.
		for argname in kwargs:
			kwargs[argname].name = argname
			self._OPTIONS[argname] = kwargs[argname]
	
	@property
	def OPTIONS(self):
		"""
		Get a list of all Configurable Options of this object.
		"""
		#return {getattr(self, attr).name: getattr(self, attr) for attr in dir(self) if attr != 'OPTIONS' and  isinstance(getattr(self, attr), Option)}
		return self._OPTIONS
	
	def pre_configure(self, obj, dictobj = None):
		"""
		Get the 'real' value of this option.
		"""
		dictobj = dictobj if dictobj is not None else obj
		
		value = self.getraw(obj)
		
		for given_option in value:
			# Check we haven't got any extra options.
			if given_option not in self.OPTIONS:
				raise Configurable_exception(obj, "option '{}' is not a valid option in group '{}'".format(given_option, self.name))
		
		return {name:option.pre_configure(obj, self.getraw(obj)) for name,option in self.OPTIONS.items()}
		
	