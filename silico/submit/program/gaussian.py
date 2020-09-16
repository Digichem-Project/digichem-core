from silico.submit.program import Program_target
from pathlib import Path
from silico.exception import Submission_error
import subprocess
import silico
from mako.lookup import TemplateLookup
from logging import getLogger
import shutil
from silico.file.fchk import Fchk_maker
from silico.exception.uncatchable import Signal_caught
from silico.submit.structure.flag import Flag
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

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.log_file_path = None
		self.chk_file_path = None
		
	@property
	def com_file_path(self):
		"""
		Path to the (ready-to-go) input file. Note that although this is known as the com file, it may infact have any extension (.com and .gjf are most common).
		"""
		return Path(self.method.calc_dir.input_directory, self.calculation.com_file_name)
	
	@property
	def default_log_file_path(self):
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
	def log_file_path(self):
		"""
		Current path to the Gaussian .log output file. Note that this location can move throughout the submission process.
		"""
		if self._log_file_path is not None:
			return self._log_file_path
		else:
			return self.default_log_file_path
		
	@log_file_path.setter
	def log_file_path(self, value):
		"""
		Change the Gaussian .log output file location, set None to reset to default.
		"""
		self._log_file_path = value
	
	@property
	def default_chk_file_path(self):
		"""
		Default path to the Gaussian checkpoint .chk file written to by Gaussian, see chk_file_path for where the chk file is currently.
		"""
		return Path(self.method.calc_dir.output_directory, self.calculation.chk_file_name)
	
	@property
	def chk_file_path(self):
		"""
		Current path to the Gaussian checkpoint .chk file. Note that this location can move throughout the submission process.
		"""
		if self._chk_file_path is not None:
			return self._chk_file_path
		else:
			return self.default_chk_file_path
		
	@chk_file_path.setter
	def chk_file_path(self, value):
		"""
		Change the Gaussian checkpoint .chk file location, set None to reset to default.
		"""
		self._chk_file_path = value
		
	@property
	def fchk_file_path(self):
		"""
		Path to the formatted checkpoint .fchk file.
		"""
		return Path(self.method.calc_dir.output_directory, self.calculation.chk_file_name).with_suffix(".fchk")
	
	def _submit_proper(self):
		"""
		Step 3/4 of the submission process, this method is called to perform submission.
		
		Main submission method; the calculation will be run here (for which this method will block, possibly for hours+).
		"""
		super()._submit_proper()
		
		# Write our input file to our calculation Input directory.
		with open(self.com_file_path, "wt") as com_file:
			com_file.write(self.calculation.com_file_body)
		
		# Make our scratch directory if we're using scratch.
		if self.calculation.scratch_directory is not None:
			try:
				# Make the folder.
				self.calculation.scratch_directory.mkdir(parents = True)
			except FileExistsError:
				# The scratch folder already existing is actually pretty serious; it's supposed to be unique to us, so if it already exists something's probably gone wrong.
				# However, some SLURM implementations automatically create our scratch dir for us; because this scenario is common (?), we currently print a warning and continue.
				# This may change in the future.
				getLogger(silico.logger_name).warning("Could not create scratch directory '{}' because it already exists; continuing anyway".format(self.calculation.scratch_directory)) 
			except Exception:
				raise Submission_error(self, "unable to create scratch directory")
			
			# Set our chk location.
			self.chk_file_path = Path(self.calculation.scratch_directory, self.calculation.chk_file_name)
			
			# If we're writing everything to scratch, set our log file there too.
			if self.calculation.scratch_options['all_output']:
				self.log_file_path = Path(self.calculation.scratch_directory, self.calculation.com_file_name).with_suffix(".log")
				
			# Note that chk_file_path is only for our reference, this is always relative to gaussian's working directory.
			# This is because the Gaussian %Chk option we use to set the chk file path is fragile and will crash on 'unusual' characters (possibly even just whitespace.)
			# So to be safe, we always use a sanitized chk path with no path separators at all, which has already been set by the Calculation_target object.
		
		# Now get our wrapper script.
		gaussian_wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/gaussian_wrapper.mako").render_unicode(program = self)
				
		# Run Gaussian (script).
		self.calc_start()
		try:
			subprocess.run(
				['bash'],
				input = gaussian_wrapper_body,
				universal_newlines = True,
				cwd = self.method.calc_dir.output_directory if self.calculation.scratch_directory is None else self.calculation.scratch_directory,
				check = True,
				# Capture output.
				stdout = subprocess.PIPE,
				stderr = subprocess.STDOUT
				)
		except Signal_caught as e:
			# We we've been told to stop (probably by SLURM because we went over time etc).
			self.calc_end(False)
			
			# Do cleanup.
			self.calculation_cleanup(False)
			
			# Continue stopping.
			raise
			
		except Exception as e:
			self.calc_end(False)
			# Something went wrong.
			# We always try and rescue any output files that might have been written.
			self.calculation_cleanup(False)
			
			# Raise.
			#raise Submission_error(self, "Error executing Gaussian '{}'".format(self.executable)) from e
			# Store for later so we can try and generate PDFs and results.
			self.error = Submission_error(self, "Error executing Gaussian '{}'".format(self.executable))
		else:
			self.calc_end()
			# All good, do cleanup.
			self.calculation_cleanup(True)
		
	def calculation_cleanup(self, success):
		"""
		Cleanup files and directory once the calculation has finished (successfully or otherwise).
		
		:param success: True if the calculation finished normally, false otherwise.
		"""
		# Set our cleanup flag.
		self.method.calc_dir.set_flag(Flag.CLEANUP)
		
		# Check to see if our main output files are in their proper (default) locations or not. Move them if necessary.
		try:
			for out_file, default_location in [
					('log_file_path', 'default_log_file_path'),
					('chk_file_path', 'default_chk_file_path')
				]:
				# Check to see if the file is in scratch or not.
				if getattr(self, out_file).resolve() != getattr(self, default_location).resolve():
					# Try and move.
					try:
						# Smart move.
						shutil.move(getattr(self, out_file), getattr(self, default_location))
					except FileNotFoundError:
						# This is safe to ignore, the file simply wasn't written.
						pass
						# We halt on other errors so we don't do something dangerous. Perhaps the most common reason why this copying might fail is because dst is out of file space, in such an instance it would a shame to contine and delete the (possibly comopleted) calc files.
					
					# And update (reset) our attribute.
					setattr(self, out_file, None)
			# If we've been asked to, we'll also try and save what remains of the scratch directory (might contain something useful).
			if self.calculation.scratch_options['use_scratch'] and ((self.calculation.scratch_options['rescue'] and not success) or self.calculation.scratch_options['keep']):
				# Check to see if there's anything in scratch.
				scratch_content = -1
				try:
					scratch_content = len(list(self.calculation.scratch_directory.iterdir()))
				finally:
					# Move our scratch dir (or at least try to) if it's not empty (or we couldn't determine if it's empty).
					if scratch_content != 0:
						shutil.move(str(self.calculation.scratch_directory), str(self.method.calc_dir.scratch_directory))
				
		finally:
			# Remove the scratch directory, forcibly if need be.
			if self.calculation.scratch_options['use_scratch'] and (success or self.calculation.scratch_options['force_delete']):
			#if success or (self.calculation.scratch_options['use_scratch'] and self.calculation.scratch_options['force_delete']):
				try:
					shutil.rmtree(self.calculation.scratch_directory)
				except FileNotFoundError:
					# This is ok, means the scratch has already been deleted (or moved).
					pass
				except Exception:
					# Some other error occurred that prevented us from deleting the scratch.
					# This is odd, annoying and possibly evidence of a bug, but shouldn't stop us from continuing.
					getLogger(silico.logger_name).warning("Failed to delete scratch directory '{}'".format(self.calculation.scratch_directory), exc_info = True)
				
			# Done cleanup.
			self.method.calc_dir.del_flag(Flag.CLEANUP)
		
	def _submit_post(self):
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
				# We can ignore not finding the file; we we're trying to get rid of it anyway.
				pass
			except Exception:
				getLogger(silico.logger_name).error("Failed to delete chk file", exc_info = True)
		
		# Use our parent to create result and report files.
		super()._submit_post()
		