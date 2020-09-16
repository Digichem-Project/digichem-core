import getpass
from pathlib import Path
from logging import getLogger
import copy
import silico

from silico.submit.structure.flag import Flag
from silico.exception.base import Submission_error
from silico.file.babel import Openbabel_converter
from silico.submit import Configurable_target, Memory
from silico.config.configurable.option import Option
from silico.config.configurable.options import Options

def _merge_silico_options(option, configurable, silico_options):
	"""
	Helper function, called to merge specific silico options with the global set.
	"""
	# Deep copy silico_options (because we're going to merge it).
	combined_silico_options = copy.deepcopy(configurable.global_silico_options)
	
	# Merge silico_options with the global options.
	combined_silico_options = configurable.merge_dict(silico_options, combined_silico_options)
				
	return combined_silico_options

class Calculation_target(Configurable_target):
	"""
	Abstract top-level class for calculation targets.
	"""
	# Top level Configurable for calculations.
	CLASS_HANDLE = ("calculation",)
	
	# Configurable options.
	programs = Option(help = "A list of programs that this calculation is compatible with", required = True, type = list)
	
	def configure(self, available_basis_sets, silico_options, **kwargs):
		"""
		Configure this calculation.
		
		:param available_basis_sets: List (possibly empty) of known external basis sets.
		:param silico_options: Global Silico options (note that this is not the per-calculation Option of the same name, but the actual global Silico options with which the former will be merged).
		"""
		self.available_basis_sets = available_basis_sets
		self.global_silico_options = silico_options
		super().configure(**kwargs)
	
	@property
	def program(self):
		"""
		Get the Program_target object that is going to run this calculation.
		"""
		return self._program
	
	@program.setter
	def program(self, value):
		"""
		The Program_target object that is going to run this calculation.
		
		:raises Configurable_target_exception: If the given program is not compatible with this calculation.
		"""
		self._set_submit_parent("_program", value)
		
	@property
	def submit_parents(self):
		"""
		Convenience property to get the 'submit parents' (programs for calculations; methods for programs) that this target supports.
		"""
		return self.programs
	
	@classmethod
	def safe_name(self, file_name):
		"""
		Get a filename safe for Gaussian (and other programs).
		
		What constitutes a safe name from Gaussian's point of view is not entirely clear, to play it safe we'll only allow alpha-numeric characters, dots and underscores.
		
		:param file_name: The file name to make safe, note that a path (containing path separators) will not maintain its structure after a call to safe_name().
		:return: The safe path name.
		"""
		# Adapted from https://stackoverflow.com/questions/7406102/create-sane-safe-filename-from-any-unsafe-string
		safe_chars = "._"
		return "".join([char if char.isalnum() or char in safe_chars else "_" for char in file_name])

	@classmethod	
	def prepare_list(self, calculation_list):
		"""
		Prepare a number of Calculation_target objects for submission.
		
		:param calculation_list: A list of 3-membered tuples (method, program, calculation) that are to be prepared. 
		"""
		first = None
		previous = None
				
		for method, program, calculation in calculation_list:
			# Expand calculation.
			for expanded_calculation in calculation.prepare():
				# Keep track of the first.
				if first is None:
					first = expanded_calculation
					
				# Setup the parent/child relationships.
				expanded_calculation.program = program
				program.method = method
				
				# Add to the chain.
				if previous is not None:
					previous.next = expanded_calculation
					
				previous = expanded_calculation
		
		# Return the first calculation in the chain.
		return first

class Concrete_calculation(Calculation_target):
	"""
	Top-level class for real calculations.
	"""
	
	CLASS_HANDLE = ()
	
	# A list of strings describing the expected input file types (file extensions) for calculations of this class. The first item of this list will be passed to obabel via the -o flag. 
	INPUT_FILE_TYPES = []
	
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		# Next is a linked list of Calculation_targets. We will call submit() on next once this object has been submitted.
		self.next = None
		self.name = None
	
	# Configurable options.
	memory = Option(help = "The amount of memory to use for the calculation", required = True, type = Memory, rawtype = str)
	num_CPUs = Option(help = "An integer specifying the number of CPUs to use for this calculation", default = 1, type = int)
	scratch_options = Options(
		help = "Options that control the use of the scratch directory",
		use_scratch = Option(help = "Whether to use a scratch directory. False will disable the scratch directory, and is not recommended", default = True, type = bool),
		path = Option(help = "Path to the top of the scratch directory. For each calculation, a sub-directory will be created under this path", default = "/scratch", type = str),
		use_username = Option(help = "Whether to create a subdirectory for each user", default = True, type = bool),
		keep = Option(help = "Whether to copy any leftover files from the scratch directory once the calculation has completed successfully", default = False, type = bool),
		rescue = Option(help = "Whether to copy files from the scratch directory if the calculation fails or stops unexpectedly", default = True, type = bool),
		force_delete = Option(help = "Whether to always delete the scratch directory at the end of the calculation, even if output files could not be safely copied", default = False, type = bool),
		all_output = Option(help = "Whether to output all files to the scratch directory. If False, only intermediate files will be written to scratch (.log will be written directly to output directory for example)", default = False, type = bool)
	)
	write_summary = Option(help = "Whether to write Silico summary text files to the 'Results' folder at the end of the calculation", default = True, type = bool)
	write_report = Option(help = "Whether to write a Silico PDF report to the 'Report' folder at the end of the calculation", default = True, type = bool)
	custom_silico_options = Option(
		"silico_options",
		help = "Silico options to overwrite for this calculation",
		default = {},
		type = dict
	)

	@property
	def silico_options(self):
		"""
		Get the specific silico options to this calculation.
		
		This property is a speed-hack; it combines global and specific silico_options the first time it is called and caches the results.
		"""
		try:
			return self._silico_options
		except AttributeError:
			# First time.
			# Deep copy silico_options (because we're going to merge it).
			self._silico_options = copy.deepcopy(self.global_silico_options)
			
			# Merge silico_options with the global options.
			self._silico_options = self.merge_dict(self.custom_silico_options, self._silico_options)
						
			return self._silico_options
	
	def prepare(self):
		"""
		Prepare this calculation target for submission.
		
		Objects should return a number (possibly 0) of calculation targets that will actually be submitted to; these targets need not be the object itself.
		
		This default implementation simply returns the same object.
		
		:return: A list of ready-to-go calculation targets.
		"""
		return (self,)
		
	@property
	def scratch_directory(self):
		"""
		Path to the scratch directory to use for this calculation. This will return None if no scratch is to be used.
		"""
		# Return None if we're not using scratch.
		if not self.scratch_options['use_scratch']:
			return None
		
		# Start with the root path.
		directory = self.scratch_options['path']
		
		# Add username if we've been asked to.
		if self.scratch_options['use_username']:
			directory += "/" + getpass.getuser()
		
		# Make a path and return.
		return Path(directory, self.program.method.unique_name)
		
	@property
	def num_CPUs(self):
		"""
		The number of CPUs to use for the calculation, this (roughly ?) translates to the number of worker processes that will be used.
		
		This property will translate 'auto' into a number of CPUs, use _num_CPUs if this is not desirable.
		"""
		if self._num_CPUs == "auto":
			return self.get_available_CPUs()
		else:
			return self._num_CPUs
	
	@num_CPUs.setter
	def num_CPUs(self, value):
		"""
		Set the number of CPUs to use for the calculation. In addition to an exact integer amount, the string "auto" can also be supplied, in which case all available CPUs will be used.
		"""
		# Set.
		self._num_CPUs = value
					
	@property
	def descriptive_name(self):
		"""
		Get a name that describes the calculation and file together.
		"""
		return "{} {}".format(self.name, self.NAME)
	
	def submit_from_file(self, output, input_file_path, *, convert = "auto", gen3D = None):
		"""
		Submit a calculation by first reading in an input file.
		
		:param output: Path to perform the calculation in.
		:param input_file_path: Path (string or pathlib.Path) to the file to submit.
		:param convert: If True, the given file will be automatically converted to an appropriate input file type. If "auto", conversion is only done if the file is not already of the correct type.
		:param gen3D: If True (and convert is True or "auto") and the loaded molecule does not have 3D coordinates, obabel will be used to generate them.
		"""
		input_format = input_file_path.suffix[1:]
		input_str = None
		
		# If we're auto converting, decide whether we should or not.
		if convert == "auto":
			if input_format == "":
				getLogger(silico.logger_name).warning("Cannot convert input file '{}' because file has no suffix (cannot determine format); the file will be submitted without conversion".format(input_file_path))
				convert = False
			else:
				convert = input_format.lower() not in self.INPUT_FILE_TYPES
				
		# Try and convert the format if we've been asked to.
		if convert:
			if input_format == "":
				raise Submission_error(self, "cannot convert input file because file has no suffix (cannot determine format)", file_name = input_file_path)
				
			try:
				input_str = Openbabel_converter.from_file(input_file_path, output_file_type = self.INPUT_FILE_TYPES[0], input_file_type = input_format, gen3D = gen3D).convert()
			except Exception:
				raise Submission_error(self, "failed to convert input file (format may not be supported)", file_name = input_file_path)
		
		# Read our file in by hand if we haven't done so far.
		if input_str is None:
			with open(input_file_path, "rt") as input_file:
				input_str = input_file.read()
		
		# Submit normally.
		return self.submit(output, input_str, input_file_path.with_suffix("").name)
		
	def submit(self, output, input_str, name):
		"""
		Submit a calculation input file to this Calculation_target.
		
		The program attribute of this Calculation_target should already be set (as should the method of that Program_target).
		"""
		self.submit_init(output, input_str, name)
		self.submit_pre()
		self.submit_proper()
		self.submit_post()
		
	def submit_init(self, output, input_str, name):
		"""
		Step 1/4 of the submission process, this method is called to set-up submission.
		Importantly, submit_init() will return before any resumable submission methods have resumed, meaning the environment during submit_init() may not reflect the final submission environment. You should typically avoid writing files that you will need later, because they may not be available in later submission methods.
		
		Inheriting classes should avoid overriding this method directly. Instead, override _submit_init().
		
		It is important to note that the normal submission order is reversed for submit_init(); the order is calculation -> program -> method.
		
		
		:param output: Path to perform the calculation in.
		:param input_str: String containing a calculation file to submit, this should be in a format appropriate for this Calculation_target.
		:param name: Name of the submitted file (should exclude extension if any). This is used eg as the base name for files created during the calculation.
		"""
		self._submit_init(output, input_str, name)
		self.program.submit_init(self)
		
	def _submit_init(self, output, input_str, name):
		"""
		Step 1/4 of the submission process, this method is called to set-up submission.
		
		Inheriting classes should override this method to perform init.
		
		This default implementation saves output and name to attributes of the same name; call super()._submit_init() in your implementation if you want to preserve this behaviour.
		If you do not call super()._submit_init(), know that several other classes will expect the output and name attributes to exist. 
		
		:param output: Path to perform the calculation in.
		:param input_str: String containing a calculation file to submit.
		:param name: Name of the submitted file (should exclude extension if any). This is used eg as the base name for files created during the calculation.
		"""
		self.output = output
		self.name = name
		self.input_file = input_str
		
	def submit_pre(self):
		"""
		Step 2/4 of the submission process, this method is called before submission begins.
		
		Inheriting classes should avoid overriding this method directly. Instead, override _submit_pre().
		
		Note the order of submission; which is method -> program -> calculation.
		"""
		try:
			self.program.submit_pre()
			self._submit_pre()
		except Exception:
			self.program.method.calc_dir.set_flag(Flag.ERROR)
			self.program.method.calc_dir.set_flag(Flag.DONE)
			raise
		
	def _submit_pre(self):
		"""
		Step 2/4 of the submission process, this method is called before submission begins.
		
		This default implementation does nothing.
		"""
		pass
		
	def submit_proper(self):
		"""
		Step 3/4 of the submission process, this method is called to perform submission.
		
		Inheriting classes should avoid overriding this method directly. Instead, override _submit_proper().
		
		Note the order of submission; which is method -> program -> calculation.
		
		The calculation will occur during this method and will have completed before this method returns (for blocking Method_targets).
		"""
		try:
			self.program.submit_proper()
			self._submit_proper()
		except Exception:
			self.program.method.calc_dir.set_flag(Flag.ERROR)
			self.program.method.calc_dir.set_flag(Flag.DONE)
			raise
		
	def _submit_proper(self):
		"""
		Step 3/4 of the submission process, this method is called to perform submission.
		
		For Calculation_targets, the calculation will already have completed before this method is called (this may change in the future).
		
		This default implementation does nothing.
		"""
		pass
		
	def submit_post(self):
		"""
		Step 4/4 of the submission process, this method is called after submission.
		
		Inheriting classes should avoid overriding this method directly. Instead, override _submit_post().
		
		Note the order of submission; which is method -> program -> calculation.
		"""
		try:
			self.program.submit_post()
		except Exception:
			self.program.method.calc_dir.set_flag(Flag.ERROR)
			raise
		finally:
			self.program.method.calc_dir.set_flag(Flag.DONE)
			
		# We don't wrap _submit_post() in flags because this method submits the next calculation; any errors that occur relate to the new submission, not this one which has just finished.
		self._submit_post()
		
	def _submit_post(self):
		"""
		Step 4/4 of the submission process, this method is called after submission.
		
		If there are any more calculations in the chain; we will call submit() on the next here.
		"""
		if self.next is not None:
			# We have another calculation to do.
			# Submit our output file.
			try:
				self.next.submit_from_file(self.output, self.program.calc_output_file_path, gen3D = False)
			except Exception:
				raise Submission_error(self, "failed to submit to next calculation in chain",  file_name = self.program.calc_output_file_path)
	
	
	
				