from silico.exception.base import Result_unavailable_error
import scipy.constants

class Result_object():
	"""
	Top level class for objects that are designed to hold calculation results.
	"""
	
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
	
		
