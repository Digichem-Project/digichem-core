from silico.submit import Configurable_target, Memory
import getpass
from pathlib import Path
from silico.exception.base import Submission_error
from silico.file.babel import Openbabel_converter
from logging import getLogger
import silico
import copy
from silico.submit.structure.flag import Flag
from silico.config.configurable import Configurable

class Calculation_target(Configurable_target):
	"""
	Top-level class for classes that implement a particular calculation
	"""
	
	# A list of strings describing the expected input file types (file extensions) for calculations of this class. The first item of this list will be passed to obabel via the -o flag. 
	INPUT_FILE_TYPES = []
	
	def __init__(
			self,
			program = None,
			*args,
			programs,
			num_CPUs = None,
			memory = None,
			scratch_options = None,
			write_summary = None,
			write_report = None,
			silico_options,
			available_basis_sets,
			**kwargs
		):
		super().__init__(*args, **kwargs)
		self.programs = programs
		self.program = program
		self.memory = Memory(memory) if memory is not None else None
		self.num_CPUs = num_CPUs
		default_scratch_options = {
			'path': '/scratch',
			'include_username': True,
			'all_output': True,
			'keep': False,
			'rescue': True,
			'force_delete': False
		}
		self.scratch_options = Configurable.merge_dict(scratch_options, default_scratch_options) if scratch_options is not None else None
		self.write_summary = write_summary if write_summary is not None else True
		self.write_report = write_report if write_report is not None else True
		self.silico_options = silico_options
		self.available_basis_sets = available_basis_sets
		# Next is a linked list of Calculation_targets. We will call submit() on next once this object has been submitted.
		self.next = None
		self.name = None
		
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
		# Start with the root path.
		directory = self.scratch_options['path']
		
		# Add username if we've been asked to.
		if self.scratch_options['include_username']:
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
				
		
	@classmethod
	def from_config(self, config, *, silico_options, available_basis_sets, **kwargs):
		"""
		Get a Configurable_target object from a config dict.
		
		Unlike other configurable targets, the from_config() method of Calculation_targets takes an addition, keyword argument.
		This argument is a config dictionary of all silico options. It will be merged with the 'silico_options' key in config, so each calculation can overwrite certain options if desired.
		This is useful, for example, for specifying an alternative alignment method for analysis for some calculations but not others.
		
		:param config: The config dict to load from.
		:param silico_options: A config dict containing all config options for configuring silico.
		:param available_basis_sets: A list of configured Basis_set targets.
 		"""
		# Deep copy both config and silico_options (because we're going to merge them).
		config = copy.deepcopy(config)
		silico_options = copy.deepcopy(silico_options)
		
		# Merge silico_options if given in the specific calc config.
		if config.get('silico_options') is not None:
			silico_options = Configurable.merge_dict( config.get('silico_options'), silico_options)
			
		# Clone our config so we don't permanently overwrite something.
		config['silico_options'] = silico_options
		
		# Continue as normal.
		return super().from_config(config, available_basis_sets = available_basis_sets, **kwargs) 
		
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
	def programs(self):
		"""
		Get a list of Program_target objects that are supported by this calculation.
		"""
		return self.submit_parents
	
	@programs.setter
	def programs(self, value):
		"""
		Set the list of Program_target objects that are supported by this calculation.
		"""
		self.submit_parents = value
		
	@property
	def descriptive_name(self):
		"""
		Get a name that describes the calculation and file together.
		"""
		return "{} {}".format(self.name, self._NAME)
	
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
			if input_format is "":
				getLogger(silico.logger_name).warning("Cannot convert input file '{}' because file has no suffix (cannot determine format); the file will be submitted without conversion".format(input_file_path))
				convert = False
			else:
				convert = input_format.lower() not in self.INPUT_FILE_TYPES
				
		# Try and convert the format if we've been asked to.
		if convert:
			if input_format is "":
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
		
		The input file to be submitted can be specified in one of two ways.
			- input_file_path: a string/path to a file that will be read and submitted.
			- input_str: a string containing a calculation file, useful if the file has already been loaded into memory.
		
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
		
	def submit_pre(self):
		"""
		Step 2/4 of the submission process, this method is called before submission begins.
		
		Inheriting classes should avoid overriding this method directly. Instead, override _submit_pre().
		
		Note the order of submission; which is method -> program -> calculation.
		"""
		self.program.submit_pre()
		self._submit_pre()
		
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
			self.program.method.calc_dir.set_flag(Flag.DONE)
		
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
		finally:
			self.program.method.calc_dir.set_flag(Flag.DONE)
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
	
	
	
				