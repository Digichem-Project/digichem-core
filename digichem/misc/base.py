import argparse
import math

def date_to_string(datetime_object):
	"""
	Convert a datetime object to a standard string representation.
	
	:param datetime_object: The date and time to convert
	:return: A string representation of datetime_object.
	"""
	return datetime_object.strftime("%d/%m/%Y at %H:%M:%S")

def timedelta_to_string(timedelta_object):
	"""
	Convert a timedelta object to a standard string representation.
	
	:param timedelta_object: The time difference to convert
	:return: A string representation of timedelta_object.
	"""
	hours = math.floor(timedelta_object.seconds / 3600)
	minutes = math.floor((timedelta_object.seconds - hours * 3600) / 60)
	return "{} days, {} hours, {} minutes".format(timedelta_object.days, hours, minutes)

class Dynamic_parent():
	"""
	A mixin class for classes that can recursively get all known children.
	"""
	
	# An iterable of strings that identify this class.
	CLASS_HANDLE = []
	
	@classmethod
	def from_class_handle(self, handle, case_insensitive = True):
		"""
		Get a class that is a child of this class from its human-readable name/handle.
		
		:raises ValueError: If a class with name could not be found.
		:param handle: The handle of the class to get (this is defined by this class itself).
		:paran case_insensitive: If true, the search is performed ignoring the cAsE of handle.
		:return: The class.
		"""
		# Our known classes.
		known_classes = self.recursive_subclasses()
		
		# Convert to lower case if we're doing a case insensitive search.
		if case_insensitive:
			handle = handle.lower()
		
		# Get the class we've been asked for.
		for known_class in known_classes:
			# Get the current classes list of handles:
			class_handles = getattr(known_class, "CLASS_HANDLE", [])
			
			# If the handle is a single string, panic.
			if isinstance(class_handles, str):
				raise TypeError("CLASS_HANDLE of class '{}' is a single string; CLASS_HANDLE should be an iterable of strings".format(known_class.__name__))
			
			# Convert to lower case if we're doing a case insensitive search.
			if case_insensitive:
				class_handles = [cls_handle.lower() for cls_handle in class_handles]
			
			# See if we have a match.	
			if handle in class_handles:
				# Got a match.
				return known_class
			
		# No class.
		raise ValueError("No {} class with name '{}' could be found".format(self.__name__, handle))
	
	@classmethod
	def recursive_subclasses(self):
		"""
		Recursively get all the subclasses of this class.
		
		:return: A set of all the classes that descend from this class.  
		"""
		def get_subclasses_worker(cls):
			return set(cls.__subclasses__()).union(
				[sub_class for top_sub_class in cls.__subclasses__() for sub_class in get_subclasses_worker(top_sub_class)]
			)
			
		return get_subclasses_worker(self)
	

class List_grouper(argparse.Action):
	"""
	Custom action class that groups lists together so we know in what order they were specified.
	"""	
	
	def __init__(self, option_strings, *args, **kwargs):
		"""
		"""
		argparse.Action.__init__(self, option_strings, *args, **kwargs)
		# We'll give ourself a name so we also group no matter which of our option_strings is used.
		self.name = [name for name in option_strings if name[:2] == "--"]
	
	def __call__(self, parser, namespace, values, option_string=None):
		"""
		"""
		grouped_list = {
			'group': self,
			'values': values
		}
		try:
			getattr(namespace, self.dest).append(grouped_list)
		except AttributeError:
			setattr(namespace, self.dest, [grouped_list])