from silico.result.result import Result_set
#from silico.file.cube import Fchk_to_cube
from silico.file.fchk import Fchk_maker
from pathlib import Path
from silico.exception import Unknown_file_type_exception, Silico_exception
import silico.file.types as file_types
from logging import getLogger
from itertools import chain
from silico.result.molecular_orbitals import Molecular_orbital_list
from silico.result.base import Result_object
import silico
from silico.parser import get_parser
from silico.result.alignment.base import Alignment

class Report():
	"""
	An enhanced report set object that contains graphics and other objects for building graphical reports.
	"""
	
	def __init__(self,
			results,
			*,
			gaussian_log_file,
			chk_file_path = None,
			fchk_file_path = None,
			options
		):
		"""
		Constructor for Report objects.
		
		:param chk_file_path: Optional path to a chk file, which can be converted to a fchk file,
		:param fchk_file_path: Optional path to an fchk file that will be used to generate cube files that are missing.
		:param options: A silico Config dictionary which contains various options that control the appearance of this report.
		"""	
		# Save our result set object.
		self.results = results
		
		# Save our image maker options.
		self.options = options
		
		# Save our aux inputs.
		self.gaussian_log_file = Path(gaussian_log_file)
		self.chk_file_path = Path(chk_file_path) if chk_file_path is not None else None
		self.fchk_file_path = Path(fchk_file_path) if fchk_file_path is not None else None
		
		# Decide which extra orbitals we want.
		self._init_get_orbitals_to_render(
			orbital_distances = options['report']['orbital_image']['orbital_distances'],
			beta_distances = options['report']['orbital_image']['beta_distances'],
			orbital_levels = options['report']['orbital_image']['orbital_levels'],
			beta_levels = options['report']['orbital_image']['beta_levels'],
			excited_state_transition_threshold = options['report']['orbital_image']['et_transition_threshold']
		)
		
		
	def _init_get_orbitals_to_render(self,
			orbital_distances = None,
			beta_distances = None,
			orbital_levels = None,
			beta_levels = None,
			excited_state_transition_threshold = None
		):
		"""
		Init helper function which decides which orbitals to render later.
		"""
		orbital_distances = orbital_distances if orbital_distances is not None else []
		beta_distances = beta_distances if beta_distances is not None else []
		orbital_levels = orbital_levels if orbital_levels is not None else []
		beta_levels = beta_levels if beta_levels is not None else []
		excited_state_transition_threshold = excited_state_transition_threshold if excited_state_transition_threshold is not None else 2
		
		try:
			# Init our list.
			self.orbitals_to_render = Molecular_orbital_list()
			
			# Now add those orbitals we've been asked to find.
			self.orbitals_to_render.extend(
				chain(
					# First alpha.
					chain.from_iterable(self.molecular_orbitals.search(HOMO_difference = HOMO_difference) for HOMO_difference in orbital_distances),
					chain.from_iterable(self.molecular_orbitals.search(level = level) for level in orbital_levels),
					
					# Now beta.
					chain.from_iterable(self.beta_orbitals.search(HOMO_difference = HOMO_difference) for HOMO_difference in beta_distances),
					chain.from_iterable(self.beta_orbitals.search(level = level) for level in beta_levels),
					
					# Also add orbitals that are involved in excited state transitions.
					chain.from_iterable(
						(transition.starting_mo, transition.ending_mo) for excited_state in self.excited_states for transition in excited_state.transitions if transition.probability >= excited_state_transition_threshold
					)
				)
			)
			
			# Now remove duplicates and reorder.
			self.orbitals_to_render = self.orbitals_to_render.ordered()
			
		except Exception:
			self.orbitals_to_render = []
			raise
	
	@classmethod
	def from_calculation_files(
			self,
			*input_files,
			gaussian_log_file = None,
			chk_file_path = None,
			fchk_file_path = None,
			discover_additional_inputs = True,
			alignment_class_name = None,
			options,
			**kwargs):
		"""
		A more intelligent constructor that can automatically determine file type(s).
		
		Keyword arguments can be used to specify files with a given type to avoid guesswork.
		
		Note that all files given should be from the same calculation, else bizarre behaviour may occur.
		
		This method is designed to be called at the start of programs and so is quite verbose.
		
		:param *input_files:  Calculation result files to be analysed, file type will be determined automatically (from file extension in most cases).
		:param gaussian_log_file: A Gaussian .log file to analyse.
		:param chk_file_path: A Gaussian .chk (checkpoint) file to analyse. This file will be safely ignored if a .fchk file is also specified.
		:param fchk_file_path: A Gaussian .fchk (formatted checkpoint) file to analyse.
		:param discover_additional_inputs: If True, the input directory will be searched for missing input files.
		:param alignment_class_name: Optional string matching the class handle of an alignment class to use for analysis. If None, this will be interpreted from options.
		:param options: A silico Config dictionary which contains various options that control the appearance of this report.
		"""
		# A dictionary of known files types.
		files = {
			file_types.gaussian_log_file: gaussian_log_file,
			file_types.gaussian_chk_file: chk_file_path,
			file_types.gaussian_fchk_file: fchk_file_path
		}
		
		# Go through our generic input_files and determine their type.
		for input_file in input_files:
			# Check type.
			key_name = None
			if file_types.gaussian_log_file.check(input_file):
				key_name = file_types.gaussian_log_file
			elif file_types.gaussian_chk_file.check(input_file):
				key_name = file_types.gaussian_chk_file
			elif file_types.gaussian_fchk_file.check(input_file):
				key_name = file_types.gaussian_fchk_file
			else:
				# Get upset if we don't recognise the given file.
				raise Unknown_file_type_exception(input_file, expected = "calculation result file")
			
			if files[key_name] is not None:
				# Give a warning if multiple of the same file type have been given.
				getLogger(silico.logger_name).warning("Ignoring input file '{}'; a {} file has already been specified: {}".format(input_file, key_name, files[key_name]))
			else:
				# Save the new file.
				files[key_name] = input_file
		
		# If we're allowed, try and get any missing inputs.
		if discover_additional_inputs:
			# Get a list of files that we actually know about.		
			concrete_files = [input_file for input_file in files.values() if input_file is not None]
			
			# See which inputs are missing.
			for file_type in files:
				if files[file_type] is None:
					# Missing, see if we can find one.
					files[file_type] = self.find_additional_inputs(concrete_files, file_type)
					
					# Print a message if we found something.
					if files[file_type] is not None:
						getLogger(silico.logger_name).info("Found '{}' in input directory; using as {} file".format(files[file_type], file_type))
		
		# Check we have enough files to actually continue.
		if files[file_types.gaussian_log_file] is None:
			raise Silico_exception("Missing required file type '{}'".format(file_types.gaussian_log_file))
		
		# Get a result set.
		# First decide on which alignment class we're using.
		alignment_class = Alignment.from_class_handle(options['alignment'] if alignment_class_name is None else alignment_class_name)
		results = get_parser(files[file_types.gaussian_log_file]).process(alignment_class)
		
		# Load emission results if we're given file names instead of results.
		for emission in ['vertical_emission_ground_result', 'adiabatic_emission_ground_result', 'emission_excited_result']:
			if kwargs.get(emission, None) is not None and not isinstance(kwargs.get(emission, None), Result_set):
				# This emission 'result' is not a result (assume it is a path); try and load it.
				try:
					#kwargs[emission] = Result_set.from_calculation_file(kwargs[emission], alignment_class_name = alignment_class_name)
					kwargs[emission] = get_parser(kwargs[emission]).process(alignment_class)
				except Exception:
					raise Silico_exception("Error loading emission result file '{}'".format(kwargs[emission]))
		
		# Add emission energies to result.
		results.add_emission(
			**kwargs
		)
				
		# Use our proper constructor.
		report = self(results, gaussian_log_file = files[file_types.gaussian_log_file], chk_file_path = files[file_types.gaussian_chk_file], fchk_file_path = files[file_types.gaussian_fchk_file], options = options)

		# All done.
		return report
		
	def __getattr__(self, attrname):
		"""
		Provide access to the result set.
		"""
		return getattr(self.results, attrname)
					
	@classmethod
	def find_additional_inputs(self, input_files, file_type):
		"""
		Search the directories of a list of input_files for additional input files of a given type.
		
		:param input_files: An iterable of input_files. The directory in which each of these files resides will be searched. None values in the list are not tolerated.
		:param file_type: A File_type object describing what file type to search for.
		"""
		# Loop through each input file we've been given.
		for input_file in input_files:
			# And loop through each possible suffix of the given file type (most only have 1...)
			for suffix in file_type.extensions:
				aux_path = Path(input_file).with_suffix(suffix)
				# See if this file exists.
				if aux_path.exists():
					return aux_path
		
		# Nothing good found.
		return None
	
	def set_file_options(self, output_dir, output_name, output_base, **kwargs):
		"""
		Set the options that will be used to create images from this object.
		
		:param output_dir: A pathlib Path object to the directory within which our files should be created.
		:param output_name: A string that will be used as the start of the file name of the files we create.
		"""
		# First, get our fchk file.
		Result_object.set_file_options(self,
			fchk_file = Fchk_maker(
				Path(output_dir, output_name + ".fchk"),
				chk_file = self.chk_file_path,
				fchk_file = self.fchk_file_path,
				output_base = output_base
			)
		)
		
		# Get a dictionary of key-word arguments we'll pass to our children.
		image_maker_options = {
			'fchk_file': self.fchk_file,
			'translations': self.alignment.translations,
 			'rotations': self.alignment.rotations,
 			'output_base': output_base,
 			'options': self.options
 		}
		
		# Call our result (sets spin density picture).
		self.results.set_file_options(output_dir, output_name, **image_maker_options)
		
		# Set our MO pictures.
		self.molecular_orbitals.set_file_options(output_dir, output_name, **image_maker_options)
		self.beta_orbitals.set_file_options(output_dir, output_name, **image_maker_options)
		
		# First set our atoms because we can then reuse the same cube file.
		self.atoms.set_file_options(output_dir, output_name, **image_maker_options)
		
		# Get said cube file.
		Result_object.set_file_options(self, structure_cube = self.atoms.get_file('structure_cube'))
		
		# Also set our alignment object.
		self.alignment.set_file_options(output_dir, output_name, cube_file = self.structure_cube, **image_maker_options)
		
		# Then dipole.
		if self.dipole_moment is not None:
			self.dipole_moment.set_file_options(output_dir, output_name, cube_file = self.structure_cube, **image_maker_options)
			
		# Now our excited states (this also sets all our TDMs, hence why we give it the cube file).
		#ground_state = Ground_state.from_results(self)
		self.excited_states.set_file_options(output_dir, output_name, cube_file = self.structure_cube, ground_state = self.ground_state, **image_maker_options)
		
		# And emission energy.
		if self.vertical_emission is not None:
			self.vertical_emission.set_file_options(output_dir, output_name, ground_state = self.ground_state, **image_maker_options)
		if self.adiabatic_emission is not None:
			self.adiabatic_emission.set_file_options(output_dir, output_name, ground_state = self.ground_state, **image_maker_options)
		
		# Convergence graphs.
		self.CC_energies.set_file_options(output_dir, output_name, **image_maker_options)
		self.MP_energies.set_file_options(output_dir, output_name, **image_maker_options)
		self.SCF_energies.set_file_options(output_dir, output_name, **image_maker_options)
		
		# And vibrations.
		self.vibrations.set_file_options(output_dir, output_name, **image_maker_options)
		
	@property
	def fchk_file(self):
		return self.get_file('fchk_file')
	
	@property
	def structure_cube(self):
		return self.get_file('structure_cube')
	
	def cleanup_intermediate_files(self):
		"""
		Remove any intermediate files that may have been created by this object.
		"""
		# Get rid of our fchk and structure cube.
		Result_object.cleanup_intermediate_files(self, 'fchk_file', 'structure_cube')
		
		# Cleanup super.
		self.results.cleanup_intermediate_files()
		
		# Now cleanup our children too.
		for attr_name, sub_object in vars(self).items():
			if isinstance(sub_object, Result_object):
				sub_object.cleanup_intermediate_files()
				
	def _write(self, output, **kwargs):
		"""
		Write the various elements of this report to file.
		
		:param output: A filename/path (can be a pathlib.Path object) that describes the base filename to use for writing files. Each file will use this name with an appropriate suffix (eg, '.HOMO.png', '.html' etc). If output has a suffix (.pdf, .html), it will be replaced in all created files.
		
		"""
		# Base directory for our images.
		image_dir = Path(output.parent, "image")
		# The base name for our images.
		#image_base_name = output.with_suffix("").name
		image_base_name = Path(self.metadata.name).with_suffix("").name
		
		# Make our output directory.
		try:
			output.parent.mkdir(parents = True)
		except FileExistsError:
			# This will happen when the dir already exists, which is fine.
			pass
		
		# And our image dir.
		try:
			image_dir.mkdir(parents = True)
		except FileExistsError:
			# This will happen when the dir already exists, which is fine.
			pass
		
		
		# Now set up all our images.
		self.set_file_options(image_dir, image_base_name, output.parent)
		
		
		
	def write(self, output, **kwargs):
		"""
		Write this HTML_report to file.
		"""
		self._write(output, **kwargs)
		
		# If we've been asked to remove intermediate files, do so.
		if self.options['report']['cleanup']:
			self.cleanup_intermediate_files()
