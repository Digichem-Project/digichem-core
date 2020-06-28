import subprocess
from subprocess import CalledProcessError
from silico.exception.base import Silico_exception
from logging import getLogger
import silico
import re
from mako.lookup import TemplateLookup
import os
import copy

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
	
	def __init__(self, input_file_path, output_file_type, *, input_file_type = None, gen3D = None):
		"""
		Constructor for the OpenBabel converter.
		
		:param gen3D: If True and the loaded molecule does not have 3D coordinates, these will be generated (this will scramble atom coordinates).
		"""
		self.input_file_path = input_file_path
		self.input_file_type = input_file_type if input_file_type is not None else input_file_path.suffix[1:]
		self.output_file_type = output_file_type
		self.gen3D = gen3D if gen3D is not None else False
		self.output_str = None
		
	@classmethod
	def from_file(self, *args, input_file_type, gen3D = None, **kwargs):
		"""
		A more powerful constructor that automatically decides which concrete class to use.
		"""
		# First, decide which class
		cls = self.get_cls(input_file_type)
		
		# And return, if gen3D is not set, we turn it on for cdx files (which have to use the naive Obabel_converter and are always 2D).
		return cls(*args, input_file_type = input_file_type, gen3D = True if gen3D is None and input_file_type.lower() == "cdx" else gen3D, **kwargs)
		
		
	def convert(self):
		"""
		Convert the input file wrapped by this class to the designated output_file_type.
		
		The converted file is saved under the output_str attribute.
 		
 		:return: The converted file (for convenience).
		"""
		self.output_str = self._convert()
		return self.output_str
		
	def _convert(self):
		"""
		Method that should actually perform file conversion.
		
		Inheriting classes should write their own implementation.
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
			
		def __init__(self, input_file_path, output_file_type, *, input_file_type, gen3D = None):
			"""
			Constructor for the OpenBabel converter.
			
			:param gen3D: If True (default) and the loaded molecule does not have 3D coordinates, these will be generated (this will scramble atom coordinates).
			"""
			super().__init__(
				input_file_path,
				output_file_type,
				input_file_type = input_file_type,
				gen3D = gen3D if gen3D is not None else True
				)
			
		
		def _convert(self):
			"""
			Convert the input file wrapped by this class to the designated output_file_type.
			
	 		:return: The converted file.
			"""
			# Get upset if input_file_type is empty (because openbabel acts weird when it is).
			if self.input_file_type is None or self.input_file_type == "":
				raise TypeError("Cannot convert file; input_file_type '{}' is None or empty".format(self.input_file_type))
			
			# Read in the molecule(s) in the given file.
			try:
				# This is a generator.
				molecules = pybel.readfile(self.input_file_type, str(self.input_file_path))
			except Exception:
				raise Silico_exception("Failed to read file '{}'".format(self.input_file_path))
			
			# Try and get the first molecule.
			try:
				molecule = next(molecules)
			except StopIteration:
				raise Silico_exception("Cannot read file '{}'; file does not contain any molecules".format(self.input_file_path))
			
			# This has been removed for now because it results in openbabel writing an error directly to output for some (all?) gaussian log files
			# Also try and get one more molecule, so we can emit a warning if we're ignoring something.
	# 		try:
	# 			next_mol = next(molecules)
	# 		except StopIteration:
	# 			# This is fine.
	# 			pass
	# 		else:
	# 			# Got another molecule, emit warning.
	# 			getLogger(silico.logger_name).warning("Molecule file '{}' contains multiple molecules; ignoring all but the first".format(self.input_file_path))
			
			# If we got a 2D (or 1D) format, convert to 3D (but warn that we are doing so.
			if molecule.dim != 3 and self.gen3D:
				# We're missing 3D coords.
				getLogger(silico.logger_name).warning("Generating 3D coordinates from {}D file '{}'; this will scramble atom coordinates".format(molecule.dim, self.input_file_path))
				molecule.localopt()
			
			# Now convert and return
			return molecule.write(self.output_file_type)
			
class Obabel_converter(Openbabel_converter):
	"""
	Wrapper class for openbabel.
	
	Unlike Babel_converter (which wraps the openbabel python interface), Openbabel_wrapper uses subprocess.run() 
	"""
	
	# The regex we'll use to check obabel converted successfully.
	obabel_success = re.compile(r"\b(?!0\b)\d*\b molecules? converted")
	
	# 'Path' to the obabel executable.
	obabel_execuable = "obabel"
		
	def __init__(self, input_file_path, output_file_type, *, gen3D = None, **kwargs):
		"""
		Constructor for the OpenBabel converter.
		
		:param gen3D: If True 3D coordinates will be generated by obabel (this will scramble atom coordinates). Because we cannot automatically determine whether the loaded file is 3D or not using obabel (we can with pybel), this option defaults to off.
		"""
		super().__init__(
			input_file_path,
			output_file_type,
			gen3D = gen3D if gen3D is not None else False,
			**kwargs
			)
			
		
	
	def _convert(self):
		"""
		Convert the input file wrapped by this class to the designated output_file_type.
		
		The converted file is saved under the output_str attribute.
 		
 		:return: The converted file (for convenience).
		"""
		# Run
		return self.run_obabel()
		
		
	def run_obabel(self):
		"""
 		Run obabel, converting the input file wrapped by this class to the designated output_file_type.
 		 		
 		:return: The converted file.
 		"""
# 		obabel_wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/obabel_wrapper.mako").render_unicode(executable = self.obabel_execuable, input_file = str(self.input_file_path), output_format = self.output_file_type, gen3D = self.gen3D)
# 		
# 		if self.gen3D:
# 			getLogger(silico.logger_name).warning("Generating 3D coordinates from file '{}'; this will scramble atom coordinates".format(self.input_file_path))
# 		
# 		done_process = subprocess.run(
# 			['bash'],
# 			input = obabel_wrapper_body,
# 			universal_newlines = True,
# 			check = True,
# 			# Capture output.
# 			stdout = subprocess.PIPE,
# 			stderr = subprocess.PIPE
# 			)
		
		sig = [
 			self.obabel_execuable,
 			str(self.input_file_path),
 			"-o", self.output_file_type
 		]
		
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
		
		if self.gen3D:
			getLogger(silico.logger_name).warning("Generating 3D coordinates from file '{}'; this will scramble atom coordinates".format(self.input_file_path))
			sig.append("--gen3D")
		
		done_process = subprocess.run(
 			sig,
 			stdout = subprocess.PIPE,
 			stderr = subprocess.PIPE,
 			universal_newlines = True,
 			check = True,
 			env = env
 		)
		
		# Sadly, openbabel doesn't appear to make use of return codes all the time.
		# We'll do basic error checking on whether our output contains a certain string.
		if not self.obabel_success.search(done_process.stderr):
			raise Silico_exception("Obabel command did not appear to complete successfully") from CalledProcessError(done_process.returncode, " ".join(done_process.args), done_process.stdout, done_process.stderr)
		
		# Return our output.
		return done_process.stdout
