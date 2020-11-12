import subprocess
from subprocess import CalledProcessError
from silico.exception.base import Silico_exception
from logging import getLogger
import silico
import re
import os
import copy
from silico.file.convert.gaussian import Gaussian_input_parser

# Try and load openbabel bindings.
HAVE_BINDINGS = False
try:
	from openbabel import pybel
	HAVE_BINDINGS = True
except ModuleNotFoundError:
	# No bindings, carry on.
	getLogger(silico.logger_name).debug("Could not load python pybel bindings; falling back to obabel executable", exc_info = True)
except Exception:
	# Some other error occurred; print an error but continue.
	getLogger(silico.logger_name).error("Found but could not load python pybel bindings; falling back to obabel executable", exc_info = True)
	

class Openbabel_converter():
	"""
	Top level class for openbabel wrappers.
	
	We support both the python interface (pybel) and running obabel directly.
	"""
	
	def __init__(self, *, input_file = None, input_file_path = None, input_file_type = None, gen3D = None):
		"""
		Constructor for the OpenBabel converter.
		
		:param input_file: A string (unicode or byte) in the format given by input_file_type that should be converted.
		:param input_file_path: Alternatively, a Path to a file that should be converted. (input_file and input_file_path are mutually exclusive).
		:param gen3D: If True and the loaded molecule does not have 3D coordinates, these will be generated (this will scramble atom coordinates).
		"""
		# Do some arg checking.
		if (input_file is None and input_file_path is None) or (input_file is not None and input_file_path is not None):
			raise TypeError(type(self).__name__ + "; exactly one of input_file or input_file_path must be given as argument")
		
		# Check we have an input type.
		if input_file_type is None and input_file_path is not None:
			input_file_type = input_file_path.suffix[1:]

		
		self.input_file = input_file
		self.input_file_path = input_file_path
		self.input_file_type = input_file_type
		self.gen3D = gen3D if gen3D is not None else False
		
	def from_com(self):
		"""
		"""
		
	@property
	def input_name(self):
		"""
		A descriptive name of the file we are converting. Works even if converting from memory.
		"""
		if self.input_file_path is not None:
			return self.input_file_path
		else:
			return "(file loaded from memory)"
		
	@classmethod
	def from_file(self, input_file_path, input_file_type, gen3D = None, **kwargs):
		"""
		A more powerful constructor that automatically decides which concrete class to use.
		"""
		# First, decide which class
		cls = self.get_cls(input_file_type)
		# Normally we use input_file_path, not input_file.
		input_file = None
		
		# Obabel can't natively read Gaussian input files (.com, .gau, .gjc, .gjf etc), but these are just fancy .xyz files so it's trivial to implement them.
		if input_file_type.lower() in ["com", "gau", "gjc", "gjf"]:
			with open(input_file_path, "rt") as gaussian_file:
				# Get and parse our gaussian file.
				gaussian_parser = Gaussian_input_parser(gaussian_file.read())
				
				# Because we've now loaded into memory, we'll unset input_file_path and use input_file instead.
				input_file_path = None
				input_file = gaussian_parser.xyz
				
				# And our file type has changed.
				input_file_type = "xyz"
		
		# And return, if gen3D is not set, we turn it on for cdx files (which have to use the naive Obabel_converter and are always 2D).
		return cls(input_file_path = input_file_path, input_file = input_file, input_file_type = input_file_type, gen3D = True if gen3D is None and input_file_type.lower() == "cdx" else gen3D, **kwargs)
		
		
	def convert(self, output_file_type):
		"""
		Convert the input file wrapped by this class to the designated output_file_type.
		
		Inheriting classes should write their own implementation.
		
		:param output_file_type: The file type to convert to.
		"""
		raise NotImplementedError("Abstract class Babel_converter does not have a convert() method defined (inheriting classes should write their own)")
		
	@classmethod
	def get_cls(self, input_file_type):
		"""
		Automatically get a concrete Babel_converter class that can be used to convert a file.
		
		If the pybel bindings are available and loaded successfully; then a Pybel_converter object will be returned,
		otherwise, the Obabel_converter will be returned (this requires openbabel to be installed and obabel to be in the path).
		
		The only exception is for the cdx format for which Obabel_converter is always returned (because of bug https://github.com/openbabel/openbabel/issues/1690 which still seems to be plaguing us in mid-2020).
		"""
		if not HAVE_BINDINGS or input_file_type.lower() == "cdx":
			return Obabel_converter
		else:
			return Pybel_converter
		
			

# Babel_convert doesn't inherit from File_converter because we are only interested in reading to/from stdin/out (which File_converter doesn't support)
#TODO: Add stdin/stdout support to File_converter
if HAVE_BINDINGS:
	class Pybel_converter(Openbabel_converter):
		"""
		Wrapper class for pybel
		
		"""
			
		def __init__(self, *args, gen3D = None, **kwargs):
			"""
			Constructor for the Pybel converter.
			
			:param gen3D: If True (default) and the loaded molecule does not have 3D coordinates, these will be generated (this will scramble atom coordinates).
			"""
			super().__init__(
				*args,
				# For Pybel, gen3D defaults to True, because we'll only use gen3d if not already in 3D.
				gen3D = gen3D if gen3D is not None else True,
				**kwargs
				)
			
		
		def convert(self, output_file_type):
			"""
			Convert the input file wrapped by this class to the designated output_file_type.
			
			:param output_file_type: The file type to convert to.
	 		:return: The converted file.
			"""
			# Get upset if input_file_type is empty (because openbabel acts weird when it is).
			if self.input_file_type is None or self.input_file_type == "":
				raise TypeError("Cannot convert file; input_file_type '{}' is None or empty".format(self.input_file_type))
			
			# Read in the molecule(s) in the given file.
			try:
				# This is a generator.
				# Use a different func depending on whether we're reading from file or memory.
				if self.input_file_path is not None: 
					# Readfile gives us an iterator of molecules...
					molecules = pybel.readfile(self.input_file_type, str(self.input_file_path))
					
					# ...but we're only ever interested in one.
					# Try and get the first molecule.
					try:
						molecule = next(molecules)
					except StopIteration:
						raise Silico_exception("Cannot read file '{}'; file does not contain any molecules".format(self.input_name))
					
				else:
					molecule = pybel.readstring(self.input_file_type, str(self.input_file))
			except Exception as e:
				raise Silico_exception("Failed to parse file '{}'".format(self.input_name)) from e
						
			# If we got a 2D (or 1D) format, convert to 3D (but warn that we are doing so.
			if molecule.dim != 3 and self.gen3D:
				# We're missing 3D coords.
				getLogger(silico.logger_name).warning("Generating 3D coordinates from {}D file '{}'; this will scramble atom coordinates".format(molecule.dim, self.input_name))
				molecule.localopt()
			
			# Now convert and return
			return molecule.write(output_file_type)
			
class Obabel_converter(Openbabel_converter):
	"""
	Wrapper class for openbabel.
	
	Unlike Babel_converter (which wraps the openbabel python interface), Openbabel_wrapper uses subprocess.run() 
	"""
	
	# The regex we'll use to check obabel converted successfully.
	#obabel_success = re.compile(r"\b(?!0\b)\d*\b molecules? converted")
	obabel_fail = re.compile(r"\b0 molecules converted")
	
	# 'Path' to the obabel executable.
	obabel_execuable = "obabel"
		
	def __init__(self, *args, gen3D = None, **kwargs):
		"""
		Constructor for the OpenBabel converter.
		
		:param gen3D: If True (default) and the loaded molecule does not have 3D coordinates, these will be generated (this will scramble atom coordinates).
		"""
		super().__init__(
			*args,
			# For Obabel, gen3D defaults to False, because we can't determine ahead of time whether we're in 3D or not.
			gen3D = gen3D if gen3D is not None else False,
			**kwargs
			)
			
		
	
	def convert(self, output_file_type):
		"""
		Convert the input file wrapped by this class to the designated output_file_type.
		 		
 		:param output_file_type: The file type to convert to.
 		:return: The converted file.
		"""
		# Run
		return self.run_obabel(output_file_type)
		
		
	def run_obabel(self, output_file_type):
		"""
 		Run obabel, converting the input file wrapped by this class to the designated output_file_type.
		 		
 		:param output_file_type: The file type to convert to. 		 		
 		:return: The converted file.
 		"""
 		# The signature we'll use to run obabel.
		sig = [self.obabel_execuable]
		
		# Add the input path if we're reading from file.
		if self.input_file_path is not None:
			sig.append(str(self.input_file_path))
			
		# Now add the input and output switches.
		sig.extend([
 			"-o", output_file_type,
 			"-i", self.input_file_type
		])
		
		if self.gen3D:
			getLogger(silico.logger_name).warning("Generating 3D coordinates from file '{}'; this will scramble atom coordinates".format(self.input_name))
			sig.append("--gen3D")
		
		# There are several openbabel bugs re. the chem draw format; one of them occurs when we are frozen and have set the BABEL_LIBDIR env variable.
		# The workaround is to temp unset BABEL_LIBDIR.
		# Get our current environment.
		env = copy.copy(os.environ)
		# Now delete BABEL_LIBDIR if we are frozen.
		if silico.frozen:
			try:
				del env['BABEL_LIBDIR']
			except KeyError:
				# The BABEL_LIBDIR isn't set.
				pass
		
		# Give our input_file as stdin if we're not reading from file.
		inputs = self.input_file
		
		# GO.
		done_process = subprocess.run(
 			sig,
 			input = inputs,
 			stdout = subprocess.PIPE,
 			stderr = subprocess.PIPE,
 			# TODO: Using universal newlines is probably not safe here; some formats are binary (.cdx etc...)
 			universal_newlines = True,
 			check = True,
 			env = env
 		)
		
		# Sadly, openbabel doesn't appear to make use of return codes all the time.
		# We'll do basic error checking on whether our output contains a certain string.
		#if not self.obabel_success.search(done_process.stderr):
		if self.obabel_fail.search(done_process.stderr):
			raise Silico_exception("obabel command did not appear to complete successfully") from CalledProcessError(done_process.returncode, " ".join(done_process.args), done_process.stdout, done_process.stderr)
		
		# Return our output.
		return done_process.stdout
