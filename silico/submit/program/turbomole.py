from pathlib import Path

from silico.submit.program.base import Program_target


class Turbomole(Program_target):
	"""
	Top level class for submitting calculations to Gaussian.
	"""
	
	CLASS_HANDLE = ("turbomole",) 

	@property
	def coord_file_path(self):
		"""
		Path to the (ready-to-go) input coord file.
		"""
		return Path(self.method.calc_dir.input_directory, "coord")