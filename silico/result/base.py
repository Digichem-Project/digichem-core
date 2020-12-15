from silico.exception.base import Result_unavailable_error
import scipy.constants

class Result_object():
	"""
	Top level class for objects that are designed to hold calculation results.
	"""
	# Some constants.
	#speed_of_light = 299792458 # m/s
	#plancks_constant = 6.62607004e-34
	#electron_volt = 1.602176634e-19 # J
	
	# A dictionary of File_maker objects that can be used to create files (eg. images) from this result object. The keys given here are the references used to access the file maker objects from the object instance.
	#file_classes = {}
	
	def __init__(self):
		"""
		"""
		self._files = {}
	
	def set_file_options(self, **kwargs):
		"""
		Set the options that will be used to create files from this object.
		
		This default implementation is intended to be used by inheriting classes; it takes a list key-word arguments (where the name of the argument is the name of the file and the value is the file object) and sets them as files of this object.
		Inheriting classes are free to specify their own behaviour if they wish.
		"""
		# TODO: Some classes don't use this interface yet (they modify _files directly which makes me sad).
		for file_name, file_obj in kwargs.items():
			self._files[file_name] = file_obj
	
	def cleanup_intermediate_files(self, *args):
		"""
		Remove any intermediate files that may have been created by this object.
		
		This default implementation is intended to be used by inheriting classes; it takes a list of file names (as found in self._files) and attempts to delete them, silently ignoring errors.
		Inheriting classes are free to specify their own behaviour if they wish.
		"""
		for file_name in args:
			try:
				self.get_file(file_name).delete(lazy = True)
			except Result_unavailable_error:
				# We're not bothered if we can't find the given file.
				pass
			
	def del_file_options(self):
		"""
		Remove the options that will be used to create files from this object.
		"""
		self._files = {}
	
	def get_file(self, file_name):
		"""
		Get a File_maker object that can make a file from this result.
		
		:param file_name: The name of the File_maker to get.
		"""
		try:
			return self._files[file_name]
		except KeyError:
			raise Result_unavailable_error("{} file '{}'".format(type(self).__name__, file_name), "Not known or set_file_options() has not yet been called")
	
	def safe_get(self, *attr_names, default = None):
		"""
		Access an attribute of this object, returning default if the attribute is None or is not available (raises Result_unavailable_error).
		"""
		# DO YOU LIKE MY B-TEAM SPAGHETTI CODE!?
		
		# Get the first name, which we'll be getattr'ing.
		first_name = attr_names[0]
		remaining_names = attr_names[1:]
		
		# Get.
		try:
			attr = getattr(self, first_name)
			
			if attr is not None and len(remaining_names) > 0:
				if isinstance(attr, Result_object):
					# We have remaining names to resolve, call attr's safe_get().
					attr = attr.safe_get(*remaining_names, default = default)
				else:
					# We have remaining names to resolve, but our attr is not a Result_set (so we can't use safe_get()).
					for remaining_name in remaining_names:
						attr = getattr(attr, remaining_name)
		except Result_unavailable_error:
			return default
		
		# And done.
		return attr
	
	@classmethod
	def wavelength_to_energy(self, wavelength):
		"""
		Convert a wavelength (in nm) to energy (in eV).
		"""
		# e = (c * h) / λ
		#return ((self.speed_of_light * self.plancks_constant) / (emission_wavelength / 1000000000)) / self.electron_volt
		return ((scipy.constants.speed_of_light * scipy.constants.Planck) / (wavelength / 1000000000)) / scipy.constants.eV
	
	@classmethod
	def energy_to_wavelength(self, energy):
		"""
		Convert an energy (in eV) to wavelength (in nm).
		"""
		# λ = (c * h) / e
		return ((scipy.constants.speed_of_light * scipy.constants.Planck) / (energy * scipy.constants.eV)) * 1000000000
	
	@classmethod
	def wavenumbers_to_energy(self, wavenumbers):
		"""
		Convert wavenumbers (in cm-1) to energy (in eV).
		"""
		return self.wavelength_to_energy((1 / wavenumbers) * 10000000)
	
		
class Result_container(list, Result_object):
	"""
	Top level class for Result_objects that hold a list of results.
	"""
	def __init__(self, *args, **kwargs):
		list.__init__(self, *args, **kwargs)
		Result_object.__init__(self)
		
	
	def __getitem__(self, key):
		try:
			return list.__getitem__(self, key)
		except IndexError:
			raise
			#raise Result_unavailable_error(type(self).__name__, "there is no item at index {}".format(key))
	
		
