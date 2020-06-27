from silico.misc.base import Dynamic_parent
from silico.exception import Result_unavailable_error, Extractor_error
from silico.misc.file_wrapper import Multi_file_wrapper
import inspect

class _Result_extractor(Dynamic_parent):
	"""
	Top-level class for all extractors.
	"""
	
	def __init__(self, ignore = None):
		self.ignore = ignore if ignore is not None else False
	

class Result_extractor(_Result_extractor):
	"""
	Top-level class for classes that extract data from Result_set objects and present it in various formats.
	"""
	FORCE_CRITERIA = False
	ALLOW_CRITERIA = False
	
	def __init__(self, *criteria, ignore = False):
		"""
		Constructor for Result_extractor objects.
		
		:param *criteria: Parameters that act as filters to select certain subsections of data. The allowed criteria depends on the inheriting class. Some classes do not allow any criteria, and some classes require criteria.
		:param ignore: Whether to tolerate Result_unavailable_error exceptions when extracting results. If False, the exception will be propagated, if True, an appropriate empty value will be returned.
		"""
		super().__init__(ignore = ignore)
		self._criteria = []
		
		if len(criteria) != 0:
			self.criteria = criteria
			
		# Get upset if we weren't given any criteria and we need some.
		if len(self.criteria) == 0 and self.FORCE_CRITERIA:
			raise Extractor_error(self, "criteria is required for this format")
		
	@property
	def criteria(self):
		"""
		The criteria, a list of variables which act a filter for this extractor, modifying the information returned.
		"""
		return self._criteria
	
	@criteria.setter
	def criteria(self, value):
		"""
		Set the criteria, a list of variables which act a filter for this extractor, modifying the information returned.
		
		Note that not all extractors support criteria, attempting to set criteria in such classes will result in Extractor_error.
		"""
		if not self.ALLOW_CRITERIA:
			raise Extractor_error(self, "criteria is not supported for this format")
			
		self._criteria = value
		
	def extract(self, result, **kwargs):
		"""
		Extract data from the given result_set.
		
		:raises Result_unavailable_error: If the data this class operates on is not available and self.ignore is False, if self.ignore is True, an empty dict is returned instead.
		"""
		data = None
		
		try:
			# We call one of our two extract functions depending on whether we were given any criteria.
			if len(self.criteria) > 0:
				# Some criteria.
				try:
					data = self._extract_with_criteria(*self.criteria, result = result, **kwargs)
				except TypeError:
					# This might be because we passed to many args to _extract, we can check by comparing how many args the function expects to how many we passed.
					if len(self.criteria) > len([param for param in inspect.signature(self._extract_with_criteria).parameters.values() if param.kind == param.POSITIONAL_OR_KEYWORD]):
						raise Extractor_error(self, "too many criteria")
						#raise TypeError("too many criteria given to {}".format(self.CLASS_HANDLE[0]))
					else:
						# Something else caused the error.
						raise
			else:
				# No criteria.
				data = self._extract(result, **kwargs)
		except Result_unavailable_error:
			# Data wasn't available.
			# If we've not been told to ignore, raise an exception.
			if not self.ignore:
				raise
		
		return data
	
	def _extract(self, result, **kwargs):
		"""
		This method is called to perform data extraction.
		
		Inheriting classes should write their own implementation.
		"""
		raise NotImplementedError("_extract is not supported by {}".format(type(self).__name__))
	
	def _extract_with_criteria(self, *criteia, result, **kwargs):
		"""
		This alternative extractor method is called to perform extraction when criteria is given.
		
		Inheriting classes that support criteria should write their own implementation.
		"""
		raise NotImplementedError("_extract_with_criteria is not supported by {}".format(type(self).__name__))
		
			
class Result_extractor_group(_Result_extractor):
	"""
 	Top-level class for classes that combine data from several Result_extractor objects.
 	
 	Groups combine data of the same format (ie, different text sections).
 	
 	Result_extractor_groups can also operate on multiple result sets.
 	"""
	
	def __init__(self, *extractor_objects, ignore = False, **kwargs):
		"""
		Constructor for Result_extractor_group objects.
		
		:param *extractor_objects: Result_extractor objects that define what data is extracted by this group. The order given is respected in output. If no objects are given, a default selection will be used.
		:param ignore: Whether to tolerate Result_unavailable_error exceptions when extracting results. If False, the exception will be propagated, if True, an appropriate empty value will be returned.
		"""
		# Get a list of default extractors if none are given.
		if len(extractor_objects) == 0:
			# Also set ignore to True if it was not already set (because we've just added all classes, who knows which ones will actually be able to get data).
			ignore = True if ignore is None else ignore
			extractor_objects = self.get_default_extractors(ignore = ignore, **kwargs)
			
		super().__init__(ignore = ignore)
		
		self.extractor_objects = extractor_objects
	
	@classmethod
	def get_default_extractors(self, **kwargs):
		"""
		Get a list of default extractor objects that can be used to convert a result set to dict format.
		
		:param **kwargs: Keyword args that will be passed as is to each extractor class to construct.
		"""
		return []
	
	def _extract(self, result, **kwargs):
		"""
		Extract data from a Result_set object.
		"""
		data = []
		
		# The order of these two loops changes how the output is grouped...
		for extractor in self.extractor_objects:
			try:
				# Add section.
				data.append(extractor.extract(result, **kwargs))
			except Result_unavailable_error:
				if self.ignore:
					data.append(None)
				else:
					raise
		return self.join(data)
	
	def join(self, extracted_results):
		"""
		Join together results from multiple extractor classes.
		
		:param extracted_results: A list of results extracted by this extractor group (one for each extractor). The format of each extracted result depends on the extractor group class.
		:return: A joined representation of the results. The format will depend on the extractor group class.
		"""
		return "\n".join([result for result in extracted_results if result != None])
	
	def join_results(self, extracted_results):
		"""
		Method called to combine a list of extracted results from multiple result sets.
		
		This default implementation mimics join().
		"""
		return self.join(extracted_results)
	
	def write_single(self, result, output_file_path):
		"""
		Convenience wrapper for write() taking only a single Result_set object and file path as arguments.
		
		:param result: A Result_set object to write.
		:param output_file_path: A path to write to. A single dash '-' is recognised as stdout.
		"""
		return self.write((result,), (output_file_path,))
			
	def write(self, results, output_file_paths):
		"""
		Extract and write data to a list of file_paths.
		
		:param results: An iterable of Result_set objects to write.
		:param output_file_paths: A list of paths to write to. A single dash '-' is recognised as stdout.
		"""	
		# Get our data first.
		data = self.extract(results)
		
		# If it's empty, don't open our file (because this would create a blank file, which is kinda weird.
		if data != None and data != "":
			for output_file_path in output_file_paths:
				with Multi_file_wrapper(output_file_path, "at") as output_file:
					output_file.write(data)
	
	def extract(self, results, **kwargs):
		"""
		Extract data from a list of Result_set objects
		"""
		data = []
		
		for result in results:
			try:
				data.append(self._extract(result, **kwargs))
			except Result_unavailable_error:
				if self.ignore:
					data.append(None)
				else:
					raise
		return self.join_results(data)
		
		
		
		
		
	