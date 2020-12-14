from silico.submit.program import Program_target
from pathlib import Path
import subprocess
import silico
from mako.lookup import TemplateLookup
from logging import getLogger
from silico.file.fchk import Fchk_maker
from silico.config.configurable.option import Option

class Gaussian(Program_target):
	"""
	Top level class for submitting calculations to Gaussian.
	"""
	
	CLASS_HANDLE = ("gaussian",)
	
	# Configurable options.
	executable = Option(help = "Name/path of the main Gaussian executable", required = True, type = str)
	root_environ_name = Option(help = "The name of the environmental variable Gaussian looks for to find 'gaussian_root'", required = True, type = str)
	gaussian_root = Option(help = "Path to the directory one above where gaussian is installed", required = True, type = str)
	gaussian_init_file = Option(help = "Path to the gaussian .profile script which is run to set-up gaussian", required = True, type = str)
		
		
	############################
	# Class creation mechanism #
	############################
	
	class _actual(Program_target._actual):
		"""
		Inner class for programs.
		"""
		
		def __init__(self, *args, **kwargs):
			"""
			Constructor for Gaussian programs.
			"""
			super().__init__(*args, **kwargs)
			
		@property
		def com_file_path(self):
			"""
			Path to the (ready-to-go) input file. Note that although this is known as the com file, it may infact have any extension (.com and .gjf are most common).
			"""
			return Path(self.method.calc_dir.input_directory, self.calculation.com_file_name)
	
		@property
		def log_file_path(self):
			"""
			Default path to the .log output file written to by Gaussian, see log_file_path for where the log file is currently.
			"""
			return Path(self.method.calc_dir.output_directory, self.calculation.com_file_name).with_suffix(".log")
	
		@property
		def calc_output_file_path(self):
			"""
			Path to the main calculation output file.
			
			For Gaussian calcs this is the .log file.
			"""
			return self.log_file_path
				
		@property
		def next_coords(self):
			"""
			Path to the output coordinate file that should be used for any subsequent calculations.
			"""
			return self.log_file_path
	
	
		@property
		def chk_file_path(self):
			"""
			Path to the Gaussian checkpoint .chk file written to by Gaussian.
			"""
			return Path(self.method.calc_dir.output_directory, self.calculation.chk_file_name)
		
		@property
		def rwf_file_path(self):
			"""
			Path to the Gaussian read-write .rwf file written to by Gaussian.
			"""
			return Path(self.method.calc_dir.output_directory, self.calculation.rwf_file_name)
			
		@property
		def fchk_file_path(self):
			"""
			Path to the formatted checkpoint .fchk file.
			"""
			return Path(self.method.calc_dir.output_directory, self.calculation.chk_file_name).with_suffix(".fchk")
	
		def pre(self):
			"""
			Pre-calculation setup for Gaussian.
			"""
			# Call parent for setup first.
			super().pre()
			
# 			# Set locations.
# 			if self.calculation.scratch_directory is not None:			
# 				# Set our chk location.
# 				self.chk_file_path = Path(self.calculation.scratch_directory, self.calculation.chk_file_name)
# 				
# 				# If we're writing everything to scratch, set our log file there too.
# 				if self.calculation.scratch_options['all_output']:
# 					self.log_file_path = Path(self.calculation.scratch_directory, self.calculation.com_file_name).with_suffix(".log")
			
			# Write our input file to our calculation Input directory.
			with open(self.com_file_path, "wt") as com_file:
				com_file.write(self.calculation.com_file_body)
	
		def calculate(self):
			"""
			Main submission method; the calculation will be run here (for which this method will block, possibly for hours+).
			"""							
			# First, get our wrapper script.
			gaussian_wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/gaussian/wrapper.mako").render_unicode(program = self)
			
			# Run Gaussian!
			subprocess.run(
				['bash'],
				input = gaussian_wrapper_body,
				universal_newlines = True,
				cwd = self.working_directory,
				check = True,
				# Capture output.
				stdout = subprocess.PIPE,
				stderr = subprocess.STDOUT
				)
		
		def cleanup(self, success):
			"""
			Cleanup files and directory once the calculation has finished (successfully or otherwise).
			
			:param success: True if the calculation finished normally, false otherwise.
			"""
# 			# Check to see if our main output files are in their proper (default) locations or not. Move them if necessary.
# 			for out_file, default_location in [
# 					('log_file_path', 'default_log_file_path'),
# 					('chk_file_path', 'default_chk_file_path')
# 				]:
# 				# Check to see if the file is in scratch or not.
# 				if getattr(self, out_file).resolve() != getattr(self, default_location).resolve():
# 					# Try and move.
# 					try:
# 						# Smart move.
# 						shutil.move(getattr(self, out_file), getattr(self, default_location))
# 					except FileNotFoundError:
# 						# This is safe to ignore, the file simply wasn't written.
# 						pass
# 						# We halt on other errors so we don't do something dangerous. Perhaps the most common reason why this copying might fail is because dst is out of file space, in such an instance it would a shame to contine and delete the (possibly comopleted) calc files.
# 					
# 					# And update (reset) our attribute.
# 					setattr(self, out_file, None)
		
		def post(self):
			"""
			Post submission method.
			"""
			# Chk/fchk management. Do this before making report (in super()) to avoid making fchk twice.
			try:
				# Create an fchk file if asked.
				if self.calculation.convert_chk:
					fchk_file = Fchk_maker(self.fchk_file_path, chk_file = self.chk_file_path)
					fchk_file.get_file()
			except Exception:
				getLogger(silico.logger_name).error("Failed to create fchk file", exc_info = True)
			else:
				try:
					# Now delete the chk file if we were asked to.
					if not self.calculation.keep_chk:
						self.chk_file_path.unlink()
				except FileNotFoundError:
					# We can ignore not finding the file; we were trying to get rid of it anyway.
					pass
				except Exception:
					getLogger(silico.logger_name).error("Failed to delete chk file", exc_info = True)
			
			# Use our parent to create result and report files.
			super().post()
			
			# Remove rwf if we've been asked to.
			try:
				if not self.calculation.keep_rwf:
					self.rwf_file_path.unlink()
			except FileNotFoundError:
				# We can ignore not finding the file; we were trying to get rid of it anyway.
				pass
			except Exception:
				getLogger(silico.logger_name).error("Failed to delete rwf file", exc_info = True)
		