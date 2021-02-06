import getpass
from pathlib import Path
import copy

from silico.submit.structure.flag import Flag
from silico.exception.base import Submission_error
from silico.submit import Configurable_target, Memory
from silico.config.configurable.option import Option
from silico.config.configurable.options import Options
from silico.file.convert import Silico_input

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
	parents = Option(
		"programs",
		help = "A list of programs that this calculation is compatible with",
		required = True,
		type = list)
	
	def configure(self, available_basis_sets, silico_options, **kwargs):
		"""
		Configure this calculation.
		
		:param available_basis_sets: List (possibly empty) of known external basis sets.
		:param silico_options: Global Silico options (note that this is not the per-calculation Option of the same name, but the actual global Silico options with which the former will be merged).
		"""
		self.available_basis_sets = available_basis_sets
		self.global_silico_options = silico_options
		super().configure(**kwargs)
	
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
	
	def expand(self):
		"""
		Expand this calculation target if it represents multiple real calcs.
		
		Objects should return a number (possibly 0) of calculation targets that will actually be submitted to; these targets need not be the object itself.
				
		:return: A list of ready-to-go calculation targets.
		"""
		raise NotImplementedError()

	@classmethod
	def link(self, calculation_list):
		"""
		Prepare a number of Calculation_target objects for submission by creating an ordered, linked list.
		
		:param calculation_list: A list of 3-membered tuples (method, program, calculation) that are to be prepared. 
		:return: A 3-membered tuple of (method, program, calculation) that is to be submitted first.
		"""
		first = None
		previous = None
		
		# These objects are class templates.
		for method_t, program_t, calculation_t in calculation_list:
			# Finalize the method and program.
			# We don't force here.
			# TODO: We shouldn't be finalizing here...
			method_t.finalize(False)
			program_t.finalize(False)
			
			# Expand calculation (because the 'calculation' could actually be a meta-calc representing multiple real calcs.
			for expanded_calculation_t in calculation_t.expand():
				# Finalize the calc.
				expanded_calculation_t.finalize(False)
				
				# Init the method, prog and calc.
				# This also links the three together.
				method = method_t()
				prog = program_t(method)
				calc = expanded_calculation_t(prog)
				
				# Keep track of the first.
				if first is None:
					first = calc
									
				# Add to the chain.
				if previous is not None:
					previous.next = calc
					
				previous = calc
		
		# Return the first calculation in the chain.
		return first

class Concrete_calculation(Calculation_target):
	"""
	Top-level class for real calculations.
	"""
	
	CLASS_HANDLE = ()
	DIRECTORY_CALCULATION = False
	
	# A list of strings describing the expected input file types (file extensions) for calculations of this class. The first item of this list will be passed to obabel via the -o flag. 
	INPUT_FILE_TYPES = []
	
	# Configurable options.
	memory = Option(help = "The amount of memory to use for the calculation", required = True, type = Memory, rawtype = str)
	_num_CPUs = Option("num_CPUs", help = "An integer specifying the number of CPUs to use for this calculation", default = 1, type = int)
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
	structure = Options(
		help = "Options that control the calculation folder structure",
		program_sub_folder = Option(help = "Whether to use a separate subfolder for each calculation program", default = False, type = bool),
		prepend_program_name = Option(help = "Whether to prepend the calculation program name to the calculation folder", default = True, type = bool),
		append_program_name = Option(help = "Whether to append the calculation program name to the calculation folder", default = False, type = bool),
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
	
	def expand(self):
		"""
		Expand this calculation target if it represents multiple real calcs.
		
		Objects should return a number (possibly 0) of calculation targets that will actually be submitted to; these targets need not be the object itself.
		
		This default implementation simply returns the same object.
		
		:return: A list of ready-to-go calculation targets.
		"""
		return (self,)
	
	
	############################
	# Class creation mechanism #
	############################
	
	class _actual(Calculation_target._actual):
		"""
		Inner class for calculations.
		"""
		
		def __init__(self, program):
			"""
			Constructor for calculation objects.
			
			:param program: A Program_target_actual object to submit to.
			"""	
			# Next is a linked list of Calculation_targets. We will call submit() on next once this object has been submitted.
			self.next = None
			self.output = None
			self.input_coords = None
			self.program = program
			self.validate_parent(program)
			self.program.calculation = self
		
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
		def molecule_name(self):
			"""
			The name of the molecule under study.
			
			This name is 'safe' for Gaussian and other sensitive programs.
			"""
			return self.safe_name(self.input_coords.name)
		
		@property
		def descriptive_name(self):
			"""
			Get a name that describes the calculation and file together.
			"""
			return "{} {}".format(self.molecule_name, self.NAME)
		
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
			
		def prepare(self, output, input_coords):
			"""
			Prepare this calculation for submission.
			
			:param output: Path to a directory to perform the calculation in.
			:param input_coords: A Silico_input object containing the coordinates on which the calculation will be performed.
			"""
			# Because we normally run the program in a different environment to where we are currently, we need to load all input files we need into memory so they can be pickled.
			self.output = output
			self.input_coords = input_coords
			
		def prepare_from_calculation(self, output, calculation):
			"""
			Alternative submission constructor that gets the input coordinates from a previously finished calc.
			
			:param output: Path to a directory to perform the calculation in.
			:param calculation: A previously submitted calculation.
			"""
			self.prepare_from_file(
				output,
				calculation.program.next_coords,
				input_format = calculation.OUTPUT_COORD_TYPE,
				molecule_name = calculation.molecule_name,
				molecule_charge = calculation.input_coords.charge,
				molecule_multiplicity = calculation.input_coords.multiplicity,
				# Don't gen3D or add H (we want to use of output coordinates exactly as-is).
				gen3D = False
			)
			
		def prepare_from_file(self,
			output,
			input_file_path,
			*,
			input_format = None,
			gen3D = None,
			molecule_name = None,
			molecule_charge = None,
			molecule_multiplicity = None):
			"""
			Alternative submission constructor that first reads in an input file.
			
			:param output: Path to a directory to perform the calculation in.
			:param input_file_path: Path (string or pathlib.Path) to the file to submit. This file should contain input coordinates for the system under study.
			:param gen3D: If True (and convert is True or "auto") and the loaded molecule does not have 3D coordinates, obabel will be used to generate them.
			:param molecule_name: Optional molecule name to use for the new calculation. If None, a name will be chosen based on the given file.
			:param molecule_charge: optional molecule charge to use for the new calculation. If None, a charge will be taken from the given file or otherwise a default will be used.
			:param molecule_multiplicity: optional molecule multiplicity to use for the new calculation. If None, a multiplicity will be taken from the given file or otherwise a default will be used.
			"""	
			# First, try and convert our given input file to the universal silico input format.
			try:
				# Load file.
				input_coords = Silico_input.from_file(input_file_path, input_format, name = molecule_name, charge = molecule_charge, multiplicity = molecule_multiplicity, gen3D = gen3D)
			except Exception:
				raise Submission_error(self, "failed to prepare input file (input format may not be supported)", file_name = input_file_path)
			
			# Prep normally.
			self.prepare(output, input_coords)
			
			
		def submit(self):
			"""
			Submit a calculation to this Calculation_target.
			"""			
			# Watch for errors so we can set flags appropriately.
			try:				
				# Do submit.
				self.program.submit()
				
			except Exception:
				# Something went wrong, set the error flag.
				self.program.method.calc_dir.set_flag(Flag.ERROR)
				self.program.method.calc_dir.set_flag(Flag.DONE)
				raise
			
			# Can't put this in finally because it will catch Submission_paused...
			self.program.method.calc_dir.set_flag(Flag.DONE)
			
			# We don't wrap our own submit_post() in flags because this method submits the next calculation; any errors that occur relate to the new submission, not this one which has just finished.		
			self.post()
			
		def post(self):
			"""
			Step 4/4 of the submission process, this method is called after submission.
			
			If there are any more calculations in the chain; we will call submit() on the next here.
			"""
			if self.next is not None:
				# We have another calculation to do.
				# Submit our output file.
				try:					
					# Prepare.
					#self.next.prepare_from_file(self.output, self.program.calc_output_file_path)
					self.next.prepare_from_calculation(self.output, self)
					
					# Go.
					self.next.submit()
					
				except Exception:
					raise Submission_error(self, "failed to submit to next calculation in chain")
	
	
class Directory_calculation_mixin():
	"""
	A class for calculations that take an existing calculation directory as input, instead of a molecule/geometry.
	"""
	
	DIRECTORY_CALCULATION = True
	
	############################
	# Class creation mechanism #
	############################
	
	class _actual():
		"""
		Inner class for calculations.
		"""
		
		def __init__(self, *args, **kwargs):
			"""
			Constructor for calculation objects.
			
			:param program: A Program_target_actual object to submit to.
			"""	
			super().__init__(*args, **kwargs)
			
			# We don't have input coords.
			del(self.input_coords)
			self.input = None
		
		@property
		def molecule_name(self):
			"""
			The name of the molecule under study.
			
			This name is 'safe' for Gaussian and other sensitive programs.
			"""
			return self.safe_name(self._molecule_name)
		
		@molecule_name.setter
		def molecule_name(self, value):
			"""
			Set the name of the molecule under study.
			"""
			self._molecule_name = value
		
		def prepare(self, output, input, *, molecule_name):
			"""
			Prepare this calculation for submission.
			
			:param output: Path to a directory to perform the calculation in.
			:param input: A directory containing existing calculation files to run.
			:param molecule_name: A name that refers to the system under study (eg, Benzene etc).
			"""
			# Because we normally run the program in a different environment to where we are currently, we need to load all input files we need into memory so they can be pickled.
			self.output = output
			self.input = input
			self.molecule_name = molecule_name 
		
		def prepare_from_file(self,
			output,
			input,
			*,
			molecule_name = None):
			"""
			Alternative submission constructor that first reads in an input file.
			
			:param output: Path to a directory to perform the calculation in.
			:param input: A directory containing existing calculation files to run.
			:param molecule_name: A name that refers to the system under study (eg, Benzene etc).s
			"""	
			# Prep normally.
			self.prepare(output, input, molecule_name = molecule_name)
			
		def prepare_from_calculation(self, output, calculation):
			"""
			Alternative submission constructor that gets the input coordinates from a previously finished calc.
			
			:param output: Path to a directory to perform the calculation in.
			:param calculation: A previously submitted calculation.
			"""
			self.prepare_from_file(
				output,
				calculation.program.method.calc_dir.output_directory,
				molecule_name = calculation.molecule_name
			)
	