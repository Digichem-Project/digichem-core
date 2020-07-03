import logging
import textwrap
from subprocess import CalledProcessError


class Variable_formatter(logging.Formatter):
	"""
	The logging formatter used by silico, the format changes depending on a logger's log_level.
	"""
	
	default_format = '%(name)s: %(levelname)s: %(message)s'
	
	def __init__(self, logger, fmt=None, datefmt=None, style='%'):
		if fmt is None:
			fmt = self.default_format
		super().__init__(fmt, datefmt, style)
		# Save our logger.
		self.logger = logger
				
	def formatException(self, exc_info):
		"""
		Format an exception.
		
		In Variable_formatter, how we do this depends on our log level.
		
		:return: The formatted exception.
		"""
		if self.logger.getEffectiveLevel() > logging.DEBUG:
			return self.exception_to_str(exc_info[1])
		else:
			return textwrap.indent(super().formatException(exc_info), "  ")
		
	@classmethod
	def process_output_to_str(self, process_exception):
		"""
		Get additional descriptive text from a CalledProcessError exception.
		
		
		"""
		output = ""
		if process_exception.stdout is not None and process_exception.stdout != "\n":
			output += "\n" + process_exception.stdout.strip()
		if process_exception.stderr is not None and process_exception.stderr != "\n":
			output += "\n" + process_exception.stderr.strip()
		return output
		
	@classmethod
	def exception_to_str(self, exception):
		"""
		Recursively get a string representation of an exception.
		
		:return: A string representation of exception and any prior exceptions.
		"""
		excstr = textwrap.indent("{}: {}".format(type(exception).__name__, exception), "\t")
		
		# If the exception is a CalledProcessError, we'll also append the stdout/stderr.
		if isinstance(exception, CalledProcessError):
			excstr += textwrap.indent(self.process_output_to_str(exception), "\t\t")
		
		if exception.__cause__ is not None:
			excstr += "\n{}".format(self.exception_to_str(exception.__cause__))
		elif exception.__context__ is not None and not exception.__suppress_context__:
			excstr += "\n{}".format(self.exception_to_str(exception.__context__))
		
		return excstr