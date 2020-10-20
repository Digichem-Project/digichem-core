from pathlib import Path

from silico.submit.program.base import Program_target
from silico.config.configurable.option import Option
import subprocess
from silico.exception.base import Submission_error
from subprocess import TimeoutExpired, CalledProcessError
import os
from mako.lookup import TemplateLookup
import silico
from silico.misc.directory import copytree
import shutil
import re


class Turbomole(Program_target):
	"""
	Top level class for submitting calculations to Turbomole.
	"""
	
	CLASS_HANDLE = ("Turbomole",) 

	define_executable = Option(help = "Name/path of the define executable", default = "define", type = str)
	root = Option(help = "Path to the directory where turbomole is installed", required = True, type = str)
	init_file = Option(help = "Path to the turbomole Config_turbo_env script which is run to set-up Turbomole", required = True, type = str)
	
	# Regex for matching the UFF section.
	UFF_SECTION = r"(\$uff\n(.|\n)+?)\n\$.+"
	
	############################
	# Class creation mechanism #
	############################
	
	class _actual(Program_target._actual):
		"""
		Inner class for programs.
		"""
		
		@property
		def environ(self):
			"""
			Environmental variables for turbomole.
			"""
			env = os.environ.copy()
			# Add the turbomole directory.
			env['TURBODIR'] = self.turbomole_root
			# Also extend path.
			
		
		@property
		def coord_file_path(self):
			"""
			Path to the input coord file.
			"""
			return Path(self.method.calc_dir.input_directory, "coord")
		
		@property
		def define_output_path(self):
			"""
			Path to the file in which define output is written.
			"""
			return Path(self.method.calc_dir.output_directory, "define.out")
		
		@property
		def turbomole_output_path(self):
			"""
			Path to the file in which turbomole output is written.
			"""
			return Path(self.method.calc_dir.output_directory, self.calculation.molecule_name).with_suffix(".log")
				
		@property
		def next_coords(self):
			"""
			Path to the output coordinate file that should be used for any subsequent calculations.
			"""
			return Path(self.method.calc_dir.output_directory, "coord")
		
		@property
		def scratch_base(self):
			"""
			The basename of the directories to use as scratch.
			When running in multi-process mode (MPI), the scratch directory given to turbomole is only used as a suggestion; each node will create its own folder by APPENDING "-001" etc to the given file name.
			
			This is annoying, so when in MPI, we make sure each of these folders is created inside the specified scratch dir.
			"""
			if self.calculation.scratch_directory is None:
				return None
			elif self.calculation.parallel_mode == "MPI":
				return Path(self.calculation.scratch_directory, "node")
			else:
				return self.calculation.scratch_directory
			
		@property
		def scratch_output(self):
			"""
			Path to the scratch folder in which the calculation will be run.
			
			When scratch is on for turbomole, two directories are considered:
			 - TURBOTMPDIR: The 'real' scratch location as understood by turbomole, possibly ignored when not in a parallel mode (SMP/MPI)?
			 - Output dir: The directory where 'control' is located, turbomole will be run with this directory as the CWD.
			Most scratch files are written according to TURBOTMPDIR, but many other files are written to the the output dir, including some that could be arguably described as scratch files.
			As such, we will set the output dir to be inside the scratch directory if all_output is True.
			This is the recommended option, but may brake MPI...
			"""
			if self.calculation.scratch_directory is None or not self.calculation.scratch_options['all_output']:
				return None
			else:
				return Path(self.calculation.scratch_directory, "Output")
		
		@property
		def define_input_path(self):
			"""
			Path to the input file used to power the define input generator.
			"""
			return Path(self.method.calc_dir.input_directory, "define.input")
		
		def define(self):
			"""
			Run setup for turbomole.
			
			Normally this involves running define, but some turbomole calcs have a different setup (eg, UFF).
			"""
			# Write define input file.
			with open(self.define_input_path, "wt") as define_input:
				define_input.write(self.calculation.define_input)
			
			# Get our wrapper script.
			wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/define_wrapper.mako").render_unicode(program = self)
						
			# Run control to generate input.
			try:
				subprocess.run(
					("bash",),
					input = wrapper_body,
					stdout = subprocess.PIPE,
					stderr = subprocess.STDOUT,
					universal_newlines = True,
					timeout = self.calculation.define_timeout,
					cwd = self.method.calc_dir.input_directory,
					check = True
				)
				
			except TimeoutExpired as e:
				# Ran out of time, probably got stuck.
				raise Submission_error(self, "Program 'define' failed to finish executing in {} s, check output file '{}' for errors".format(self.calculation.define_timeout, self.define_output_path)) from e	
			
			except CalledProcessError as e:
				# Something went wrong.
				e.__context__ = None
				raise Submission_error(self, "define did not exit successfully:\n{}".format(e.stdout)) from e
			
			except Exception as e:
				# Something else.
				raise e from None
			
		def UFF_define(self):
			"""
			Run setup for turbomole UFF.
			
			Unlike real turbomole calcs, UFF doesn't use define for setup (but we still take some options from the control file weirdly...)
			"""
			# Get our wrapper script.
			wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/turbomole_wrapper.mako").render_unicode(program = self)
						
			# Run control to generate input.
			try:
				subprocess.run(
					("bash",),
					input = wrapper_body,
					stdout = subprocess.PIPE,
					stderr = subprocess.STDOUT,
					universal_newlines = True,
					cwd = self.method.calc_dir.input_directory,
					check = True
				)
				
			except CalledProcessError as e:
				# Something went wrong.
				e.__context__ = None
				raise Submission_error(self, "UFF (setup) did not exit successfully:\n{}".format(e.stdout)) from e
			
			except Exception as e:
				# Something else.
				raise e from None
			
			# Get our custom uff section.
			uff_input = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/uff.mako").render_unicode(program = self)
			
			# We now need to manually alter some options in control (sigh...).
			# Open control.
			with open(Path(self.method.calc_dir.input_directory, "control"), "r+") as control_file:
				# Read in existing control.
				control = control_file.read()
				
				# Replace existing uff section with our custom one.
				control = re.sub(self.UFF_SECTION, uff_input, control)
				
				# Seek back to start of file.
				control_file.seek(0)
				
				# Write modified control.
				control_file.write(control)
				
				# Remove any leftover data.
				control_file.truncate()
			
		
		def pre(self):
			"""
			Pre-calculation setup for Turbomole.
			"""
			# Call parent for setup first.
			super().pre()
					
			# Write our input file to our calculation Input directory.
			with open(self.coord_file_path, "wt") as coord_file:
				coord_file.write(self.calculation.input_coords)
			
			# Run define.
			if "Turbomole-UFF" in self.calculation.CLASS_HANDLE:
				# Use alternative UFF setup.
				self.UFF_define()
			else:
				# Normal setup.
				self.define()
			
			# If we're using a scratch output dir, create it now.
			if self.scratch_output is not None:
				try:		
					# Copy our input dir to the scratch version.
					copytree(self.method.calc_dir.input_directory, self.scratch_output)
				except Exception as e:
					raise Submission_error(self, "Failed to make scratch subdirectory") from e
				
				# We're using scratch to write all our output.
				# So that we can keep track of the calc as it runs, we'll temporarily turn our real output path into a symlink to the scratch output dir.
				#os.symlink(Path(self.method.calc_dir.output_directory, "job.last").resolve(), Path(self.scratch_output, "job.last"), target_is_directory = False)
			else:
				# Copy input to normal output folder.
				copytree(self.method.calc_dir.input_directory, self.method.calc_dir.output_directory)
			
		def calculate(self):
			"""
			Main submission method; the calculation will be run here (for which this method will block, possibly for hours+).
			"""										
			# Get our wrapper script.
			wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/turbomole_wrapper.mako").render_unicode(program = self)
			
			# Decide on where we are running.
			if self.scratch_output is not None:
				cwd = self.scratch_output
			else:
				cwd = self.method.calc_dir.output_directory
			
			
			# Run Turbomole!
			subprocess.run(
				("bash",),
				input = wrapper_body,
				stdout = subprocess.PIPE,
				stderr = subprocess.STDOUT,
				universal_newlines = True,
				cwd = cwd,
				check = True
			)
	
		def cleanup(self, success):
			"""
			Cleanup files and directory once the calculation has finished (successfully or otherwise).
			
			:param success: True if the calculation finished normally, false otherwise.
			"""
			# If we were using scratch output, copy back now.
			if self.scratch_output is not None:
				# Delete the job.last symlink.
# 				try:
# 					os.unlink(Path(self.scratch_output, "job.last"))
# 				except FileNotFoundError:
# 					# ok.
# 					pass
				
				# Copy.
				copytree(self.scratch_output, self.method.calc_dir.output_directory)
				
				# Delete the scratch output.
				shutil.rmtree(self.scratch_output)
				
		def post(self):
			# Not implemented for turbomole yet.
			pass
		
			# If we got an error during the calc, re-raise it now.
			if self.error is not None:
				raise self.error
			
			
# class Turbomole_UFF(Turbomole):
# 	"""
# 	Top level class for submitting calculations to Turbomole.
# 	"""
# 	
# 	CLASS_HANDLE = ("Turbomole-UFF",) 
# 
# 	define_executable = Option(help = "Name/path of the define executable", default = "define", type = str)
# 	root = Option(help = "Path to the directory where turbomole is installed", required = True, type = str)
# 	init_file = Option(help = "Path to the turbomole Config_turbo_env script which is run to set-up Turbomole", required = True, type = str)
# 	
# 	# The regex for setting the UFF max iterations,
# 	MAXITER_MATCH = re.compile(r"\$uff\n +[0-9]+ +[0-9]+ +[0-9]+")
# 	
# 	# The regex for setting the convergence criteria.
# 	CONV_MATCH = re.compile(r"\$uff\n +[0-9]+ +[0-9]+ +[0-9]+")
# 	
# 	# Regex for matching the UFF section.
# 	UFF_SECTION = r"(\$uff\n(.|\n)+?)\n\$.+"
# 	
# 	############################
# 	# Class creation mechanism #
# 	############################
# 	
# 	class _actual(Program_target._actual):
# 		"""
# 		Inner class for programs.
# 		"""
# 		
# 		def define(self):
# 			"""
# 			Run setup for turbomole UFF.
# 			
# 			Unlike real turbomole calcs, UFF doesn't use define for setup (but we still take some options from the control file weirdly...)
# 			"""
# 			# Get our wrapper script.
# 			wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/turbomole_wrapper.mako").render_unicode(program = self)
# 						
# 			# Run control to generate input.
# 			try:
# 				subprocess.run(
# 					("bash",),
# 					input = wrapper_body,
# 					stdout = subprocess.PIPE,
# 					stderr = subprocess.STDOUT,
# 					universal_newlines = True,
# 					cwd = self.method.calc_dir.input_directory,
# 					check = True
# 				)
# 				
# 			except CalledProcessError as e:
# 				# Something went wrong.
# 				e.__context__ = None
# 				raise Submission_error(self, "UFF (setup) did not exit successfully:\n{}".format(e.stdout)) from e
# 			
# 			except Exception as e:
# 				# Something else.
# 				raise e from None
# 			
# 			# Get our custom uff section.
# 			uff_input = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/uff.mako").render_unicode(program = self)
# 			
# 			# We now need to manually alter some options in control (sigh...).
# 			# Open control.
# 			with open(Path(self.method.calc_dir.input_directory, "control"), "r+") as control_file:
# 				# Read in existing control.
# 				control = control_file.read()
# 				
# 				# Replace existing uff section with our custom one.
# 				control = re.sub(self.UFF_SECTION, uff_input, control)
# 				
# 				# Seek back to start of file.
# 				control_file.seek(0)
# 				
# 				# Write modified control.
# 				control_file.write(control)
# 				
# 				# Remove any leftover data.
# 				control_file.truncate()
# 	
			