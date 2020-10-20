from silico.submit import Configurable_target
import silico
from pathlib import Path
from silico.report.pdf import PDF_report
from silico.extract.text import Text_summary_group_extractor
from silico.extract.csv import Long_CSV_group_extractor
from silico.extract.long import Atoms_long_extractor, Orbitals_long_extractor,\
	Beta_long_extractor, SCF_long_extractor, MP_long_extractor, CC_long_extractor,\
	Excited_state_long_extractor, Excited_state_transitions_long_extractor,\
	TDM_long_extractor, Vibrations_long_extractor,\
	Absorption_spectrum_long_extractor, Absorption_energy_spectrum_long_extractor,\
	IR_spectrum_long_extractor
from logging import getLogger
from timeit import default_timer as timer
import datetime
from silico import misc
from silico.submit.structure.flag import Flag
from silico.config.configurable.option import Option
from silico.file.babel import Openbabel_converter
from silico.exception.base import Submission_error
from silico.exception.uncatchable import Signal_caught
import shutil

class Program_target(Configurable_target):
	"""
	Top-level class for classes that implement submission to a calculation program.
	"""
	
	CLASS_HANDLE = ("program",)
	
	# Configurable options.
	parents = Option("methods", help = "A list of methods that this program is compatible with", required = True, type = list)
	
	############################
	# Class creation mechanism #
	############################
	
	class _actual(Configurable_target._actual):
		"""
		Inner class for programs.
		"""
		
		def __init__(self, method):
			"""
			Constructor for program objects.
			
			:param method: A Method_target_actual object that is going to be submitted to.
			"""
			# Set our method.
			self.method = method
			self.validate_parent(method)
			# Let our method know who we are.
			self.method.program = self
			# We don't have a calculation yet.
			self.calculation = None
		
			# A Result_set object that will be populated once this calculation has completed.
			self.result = None
			self.start_time = None
			self.end_time = None
			self.duration = None
			
			# An exception caught during the calculation. This will be re-raised once cleanup and analysis has been finished (attempted).
			self.error = None
			
		def submit(self):
			"""
			Submit this program; running the specified calculation.
			
			Inheriting classes can normally avoid overwriting this method by defining a calculate() method instead.
			
			:param calculation: A Calculation_target_actual object that is going to be submitted.
			"""
			# Call method submit first to setup the environment.
			self.method.submit()
			
			
			# Run Program (script).
			self.start()			
				
			try:
				# Pre-calc (write input files etc).
				self.pre()
				
				# Go.
				self.calculate()
				
			except (Signal_caught, KeyboardInterrupt):
				# We've been told to stop (probably by SLURM because we went over time etc).
				self.end(False)
				
				# Continue stopping.
				raise
				
			except Exception as e:
				# Something went wrong.
				self.end(False)
				
				# Raise.
				# Store for later so we can try and generate PDFs and results.
				self.error = Submission_error(self, "Error executing calculation program")
				self.error.__cause__ = e
				
			else:
				# Finished normally.
				self.end(True)
				
			# Post calc (write result files etc).
			self.post()
	
		@property
		def success(self):
			"""
			Get whether the calculation has completed successfully.
			
			This property has one of 3 values:
				- True, if the calculation finished successfully.
				- False, if the calculation finished normally but the quantum chemistry program signalled an error (eg, failed convergence, unknown keyword etc).
				- None, if the calculation has not yet finished or if the calculation failed catastrophically (and we could not load enough results to determine success or not). 
			"""
			if self.result is None:
				return None
			else:
				return self.result.safe_get('metadata', 'calc_success')
						
			
		def start(self):
			"""
			Signals the start of a calculation, this method is called automatically by submit().
			"""
			# Set the start and running flags.
			self.method.calc_dir.set_flag(Flag.STARTED)
			self.method.calc_dir.set_flag(Flag.RUNNING)
			
			# Save our start time.
			self.start_timer = timer()
			self.start_date = datetime.datetime.now()
			
			# Log.
			getLogger(silico.logger_name).info("Calculation start on {} ".format(misc.date_to_string(self.start_date)))
		
		def end(self, success = True):
			"""
			Signals the end of a calculation, this method is called automatically by submit().
			
			:param success: Should be False if the program finished with error.
			"""
			try:		
				# Set our error flag if something went wrong
				if not success:
					self.method.calc_dir.set_flag(Flag.ERROR)
				
				self.end_timer = timer()
				self.end_date = datetime.datetime.now()
				
				# Assemble our string.
				message = "Calculation end on {}" if success else "Abnormal calculation end on {}"
				message = message.format(misc.date_to_string(self.end_date)) 
				
				if success:
					getLogger(silico.logger_name).info(message)
				else:
					getLogger(silico.logger_name).error(message)
					
				# Work out how much time has passed.
				self.duration = datetime.timedelta(seconds = self.end_timer - self.start_timer)
				getLogger(silico.logger_name).info("Calculation duration: {} ({} total seconds)".format(misc.timedelta_to_string(self.duration), self.duration.total_seconds()))
				
				###########
				# Cleanup #
				###########
				
				# Unset our running flag.
				self.method.calc_dir.del_flag(Flag.RUNNING)
				# Set our cleanup flag.
				self.method.calc_dir.set_flag(Flag.CLEANUP)
			
				# Call user specified cleanup.
				self.cleanup(success)
				
				# If we've been asked to, we'll try and save what remains of the scratch directory (might contain something useful).
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
		
		@property
		def calc_output_file_path(self):
			"""
			Path to the main calculation output file.
			"""
			raise NotImplementedError()
		
		def pre(self):
			"""
			Perform pre-setup for a calculation.
			"""
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
	
		def calculate(self):
			"""
			This method should be implemented in child classes to perform the designated calculation.
			"""
			raise NotImplementedError()
		
		def cleanup(self, success):
			"""
			This method should be implemented in child classes to perform cleanup after the calculation.
			
			:param success: True if the calculation finished normally.
			"""
			pass
		
		def post(self):
			"""
			Perform post analysis and cleanup, this method is called after a calculation has finished.
			"""
			# Set Flag.
			self.method.calc_dir.set_flag(Flag.POST)
			
			# First, (try) and load our results.
			# We'll actually load a report object because it is the same as a Result_set but with some extra methods which we might need to write a report. Saves us loading results twice.
			# We need to know whether the calculation was successful or not, so we make no effort to catch exceptions here.
			try:
				self.result = PDF_report.from_calculation_files(
					date = self.end_date,
					duration = self.duration,
					gaussian_log_file = self.calc_output_file_path,
					prog_version = silico.version,
					options = self.calculation.silico_options
					)
			except Exception:
				# No good.
				self.method.calc_dir.set_flag(Flag.ERROR)
				raise
			else:
				# See if our calculation was successful or not.
				if self.result.metadata.calc_success and self.error is None:
					self.method.calc_dir.set_flag(Flag.SUCCESS)
				else:
					# No good.
					self.method.calc_dir.set_flag(Flag.ERROR)
					
				# Also check optimisation convergence.
				if self.result.metadata.optimisation_converged is not None and "Optimisation" in self.result.metadata.calculations:
					if self.result.metadata.optimisation_converged:
						# Converged.
						self.method.calc_dir.set_flag(Flag.CONVERGED)
					else:
						self.method.calc_dir.set_flag(Flag.NOT_CONVERGED)
			
			# If we've been asked to write result files, do so.
			try:
				if self.calculation.write_summary:
					self.write_summary_files()
			except Exception:
				getLogger(silico.logger_name).warning("Failed to write calculation result summary files", exc_info = True)
				
			# Write XYZ file.
			try:
				self.write_XYZ()
			except Exception:
				getLogger(silico.logger_name).warning("Failed to write XYZ result file", exc_info = True)
				
			# Similarly, if we've been asked to write a report, do that.
			# TEMP: Don't write reports if we fail, see #29.
			if self.error is None:
				try:
					if self.calculation.write_report:
						self.write_report_files()
				except Exception:
					getLogger(silico.logger_name).warning("Failed to write calculation report", exc_info = True)
			else:
				getLogger(silico.logger_name).info("Skipping report generation because calculation did not finish successfully")
				
			# Delete Flag.
			self.method.calc_dir.del_flag(Flag.POST)
				
			# If we got an error during the calc, re-raise it now.
			if self.error is not None:
				raise self.error
		
		def write_summary_files(self):
			"""
			Write text result files (like with cresult) from this calculation.
			"""
			# First, make our result directory.
			try:
				self.method.calc_dir.result_directory.mkdir()
			except FileExistsError:
				pass
			
			# First, write a text summary.
			Text_summary_group_extractor(ignore = True, config = self.calculation.silico_options).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name +".summary"))
			
			# Atoms.
			Long_CSV_group_extractor(Atoms_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".atoms.csv"))
			
			# Alpha. We'll use a different file name depending on whether we are restricted or unrestricted.
			Long_CSV_group_extractor(Orbitals_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".{}.csv".format("orbitals" if self.result.metadata.orbital_spin_type == "restricted" else "alpha")))
			# And beta (which can only be for unrestricted.
			Long_CSV_group_extractor(Beta_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".beta.csv"))
			
			# Energies.
			Long_CSV_group_extractor(SCF_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".SCF.csv"))
			Long_CSV_group_extractor(MP_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".MP.csv"))
			Long_CSV_group_extractor(CC_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".CC.csv"))
			
			# Excited states.
			Long_CSV_group_extractor(Excited_state_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".ES.csv"))
			Long_CSV_group_extractor(Excited_state_transitions_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".transitions.csv"))
			Long_CSV_group_extractor(TDM_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".TDM.csv"))
			Long_CSV_group_extractor(Absorption_spectrum_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".UV-Vis.csv"))
			Long_CSV_group_extractor(Absorption_energy_spectrum_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".absorptions.csv"))
			
			# And vibrations.
			Long_CSV_group_extractor(Vibrations_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".vibrations.csv"))
			Long_CSV_group_extractor(IR_spectrum_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".IR.csv"))
			
		def write_report_files(self):
			"""
			Write report files (like with creport) from this calculation.
			"""
			# The full report.
			self.result.write(Path(self.method.calc_dir.report_directory, self.calculation.molecule_name + ".pdf"))
			# And atoms.
			self.result.write(Path(self.method.calc_dir.report_directory, self.calculation.molecule_name + ".atoms.pdf"), report_type = "atoms")
			
		def write_XYZ(self):
			"""
			Write an XYZ file from the finished calculation results.
			"""
			# Get our converter.
			conv = Openbabel_converter.from_file(self.calc_output_file_path, output_file_type = "xyz", input_file_type = self.calc_output_file_path.suffix[1:])
			
			# Open output file.
			with open(Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".xyz"), "wt") as xyz_file:
				# Write.
				xyz_file.write(conv.convert())
		
		