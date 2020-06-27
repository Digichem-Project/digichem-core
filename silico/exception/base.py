# Exceptions and other errors.
import inspect
import silico
#import silico.submit.method
#import silico.submit.program

class Silico_exception(Exception):
	"""
	General silico exception.
	"""
	
	def __init__(self, message):
		self.message = message
		
	def __str__(self, *args, **kwargs):
		return self.message
	

class Result_unavailable_error(Silico_exception):
	"""
	Exception for when a requested result is not available (because it could not be found in the calculation results for example).
	"""
	
	def __init__(self, result_name, reason = "result could not be found"):
		"""
		Constructor for Result_unavailable_error objects.
		
		:param result_name: The name of the result that is unavailable.
		:param reason: Optional message explaining why the result is unavailable. If not given, a default message will be used.
		"""
		self.result_name = result_name
		self.reason = reason
		
	def __str__(self, *args, **kwargs):
		"""
		Stringify this error.
		"""
		return "'{}' is not available; {}".format(self.result_name, self.reason)
		
class File_maker_exception(Silico_exception):
	"""
	Exception for when a file cannot be made/rendered for whatever reason.
	"""
	
	def __init__(self, file_maker, reason = ""):
		"""
		Constructor for File_maker_exception objects.
		
		:param file_maker: The file_maker object where the exception occurred.
		:param reason: Optional string describing why the exception occurred.
		"""
		self.file_maker = file_maker
		self.reason = reason
		
	def __str__(self, *args, **kwargs):
		"""
		Stringify this error.
		"""
		return "Error making '{}' file '{}'; {}".format(type(self.file_maker).__name__, self.file_maker.output, self.reason)
	
class Unknown_file_type_exception(Silico_exception):
	"""
	Exception for when a file is given but its type cannot be determined.
	"""
	
	def __init__(self, file_path, expected = None):
		"""
		Constructor for Unknown_file_type_exception objects.
		
		:param file_path: String-like path of the file that is unrecognised.
		:param expected: An optional string-like representing the type of file that was expected. 
		"""
		self.file_path = file_path
		self.expected = expected
		
	def __str__(self, *args, **kwargs):
		"""
		Stringify this error.
		"""
		err_str = "Unknown file type '{}'".format(self.file_path)
		if self.expected is not None:
			err_str = "{}; expected file of type '{}'".format(err_str, self.expected)
		return err_str
	
class Config_loader_exception(Silico_exception):
	"""
	Exceptions that occur during reading of a config file.
	"""
	
	def __init__(self, config, reason):
		self.config = config
		self.reason = reason
		
	def __str__(self, *args, **kwargs):
		"""
		Stringify this error.
		"""
		return "Unable to load config file '{}'; {}".format(self.config, self.reason)
		
class Bad_config_exception(Silico_exception):
	"""
	Exceptions for when a config file can be read/parsed, but there is some higher level error that prevents it from being used.
	"""
	
	def __init__(self, reason):
		self.reason = reason
		
	def __str__(self, *args, **kwargs):
		"""
		Stringify this error.
		"""
		return "Bad config; {}".format(self.reason)
	
class Configurable_target_exception(Silico_exception):
	"""
	Exceptions for missing/wrong/conflicting arguments being supplied to Configurable_target objects.
	"""
		
	def __init__(self, object_or_class, reason):
		self.config_type = object_or_class.__name__ if inspect.isclass(object_or_class) else "{}/{}".format(object_or_class.CONFIG_TYPE, object_or_class.CONFIG_NAME)
		self.reason = reason
		
	def __str__(self, *args, **kwargs):
		"""
		Stringify this error.
		"""
		return "Error in  '{}'; {}".format(self.config_type, self.reason)
	
class Submission_error(Silico_exception):
	"""
	Exceptions for when an error occurs during calculation submission.
	"""
	
	def __init__(self, calculation, reason, file_name = None):
		"""
		Constructor for Submission_error exception objects.
		
		:param calculation: The calculation that was in process of being submitted when the error occured. This can be any one of a Method_target, Program_target or Calculation_target.
		:param reason: String describing why the error occured.
		"""		
		# Do some quick type checking.
		if calculation.CONFIG_TYPE == "submit_methods":
			# 'Calculation' is actually a Method_target.
			calculation = calculation.program.calculation
		elif calculation.CONFIG_TYPE == "submit_programs":
			# 'Calculation' is actually a Program_target.
			calculation = calculation.calculation
		
		# Decide on file name.
		if file_name is None:
			file_name = calculation.name
		self.file_name = file_name
		
		self.calculation = calculation
		self.reason = reason
		
	def __str__(self, *args, **kwargs):
		"""
		Stringify this error.
		"""
		return "Error submitting file '{}' to '{}/{}/{}'; {}".format(self.file_name, self.calculation.program.method._CONFIG_NAME, self.calculation.program._CONFIG_NAME, self.calculation._CONFIG_NAME, self.reason)
	
class Extractor_error(Silico_exception):
	"""
	Exceptions for when an error occurs during result extraction.
	"""
	
	def __init__(self, extractor, reason):
		"""
		Constructor for Extractor_error exception objects.
		"""
		self.extractor = extractor
		self.reason = reason
		
	def __str__(self):
		"""
		Stringify this error.
		"""
		return "{} ({}); {}".format(type(self.extractor).__name__, ", ".join(getattr(self.extractor, 'CLASS_HANDLE', [])), self.reason)
	