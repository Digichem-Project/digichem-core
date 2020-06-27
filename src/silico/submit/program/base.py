from silico.submit import Configurable_target
import silico
from pathlib import Path
from silico.report.pdf import PDF_report
from silico.extract.text import Text_summary_group_extractor
from silico.extract.csv import Long_CSV_group_extractor
from silico.extract.long import Atoms_long_extractor, Orbitals_long_extractor,\
	Beta_long_extractor, SCF_long_extractor, MP_long_extractor, CC_long_extractor,\
	Excited_state_long_extractor, Excited_state_transitions_long_extractor,\
	TDM_long_extractor, Vibrations_long_extractor
from logging import getLogger
from timeit import default_timer as timer
import datetime
from silico import misc

class Program_target(Configurable_target):
	"""
	Top-level class for classes that implement submission to a calculation program.
	"""
	
	def __init__(self, method = None , *args, methods, **kwargs):
		super().__init__(*args, **kwargs)
		self.methods = methods
		self.method = method
		
	@property
	def method(self):
		"""
		Get the Method_target object that is going to run this program.
		"""
		return self._method
	
	@method.setter
	def method(self, value):
		"""
		The Method_target object that is going to run this program.
		
		:raises Configurable_target_exception: If the given method is not compatible with this program.
		"""
		self._set_submit_parent("_method", value)
		
	@property
	def methods(self):
		"""
		Get a list of Method_target objects that are supported by this program.
		"""
		return self.submit_parents
	
	@methods.setter
	def methods(self, value):
		"""
		Set the list of Method_target objects that are supported by this program.
		"""
		self.submit_parents = value
		
	@property
	def silico_options(self):
		"""
		Get the global silico options dictionary (this is actually found under our Calculation_target).
		"""
		return self.calculation.silico_options
		
	@property
	def calc_output_file_path(self):
		"""
		Path to the main calculation output file.
		"""
		raise NotImplementedError()
	
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
		
	def submit_init(self, calculation):
		"""
		Step 1/4 of the submission process, this method is called to set-up submission.
		Importantly, submit_init() will return before any resumable submission methods have resumed, meaning the environment during submit_init() may not reflect the final submission environment. You should typically avoid writing files that you will need later, because they may not be available in later submission methods.
		
		Inheriting classes should avoid overriding this method directly. Instead, override _submit_init().
		
		It is important to note that the normal submission order is reversed for submit_init(); the order is calculation -> program -> method.
		
		:param calculation: A Calculation_target object that is going to be submitted. This Calculation_target object will have completed submit_init() before this method is called.
		"""
		self._submit_init(calculation)
		self.method.submit_init(self)
	
	def _submit_init(self, calculation):
		"""
		Step 1/4 of the submission process, this method is called to set-up submission.
		
		Inheriting classes should override this method to perform init.
		
		This default implementation saves calculation to an attribute of the same name. Call super()._submit_init() in your implementation if you want this behaviour.
		If you do not call super()._submit_init(), know that several other classes will expect the calculation attribute to exist. 
		
		:param calculation: A Calculation_target object that is going to be submitted. This Calculation_target object will have completed submit_init() before this method is called.
		"""
		self.calculation = calculation
		
		# A Result_set object that will be populated once this calculation has completed.
		self.result = None
		self.start_time = None
		self.end_time = None
		self.duration = None
		
	def submit_pre(self):
		"""
		Step 2/4 of the submission process, this method is called before submission begins.
		
		Inheriting classes should avoid overriding this method directly. Instead, override _submit_pre().
		
		Note the order of submission; which is method -> program -> calculation.
		"""
		self.method.submit_pre()
		self._submit_pre()
		
	def _submit_pre(self):
		"""
		Step 2/4 of the submission process, this method is called before submission begins.
		
		This default implementation does nothing.
		"""
		pass
	
	def calc_start(self):
		"""
		Implementing classes should call this function immediately before running the computational chemistry software.
		"""
		self.start_timer = timer()
		self.start_date = datetime.datetime.now()
		
		getLogger(silico.logger_name).info("Calculation start on {} ".format(misc.date_to_string(self.start_date)))
		
	def calc_end(self, success = True):
		"""
		Implementing classes should call this function immediately after running the computational chemistry software.
		
		:param success: Should be False if the program finished with error.
		"""
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
# 		hours = math.floor(self.duration.seconds / 3600)
# 		minutes = math.floor((self.duration.seconds - hours * 3600) / 60)
# 		getLogger(silico.logger_name).info("Calculation duration: {} days, {} hours, {} minutes  ({} total seconds)".format(self.duration.days, hours, minutes, self.duration.total_seconds()))
		getLogger(silico.logger_name).info("Calculation duration: {} ({} total seconds)".format(misc.timedelta_to_string(self.duration), self.duration.total_seconds()))
	
	def submit_proper(self):
		"""
		Step 3/4 of the submission process, this method is called to perform submission.
		
		Inheriting classes should avoid overriding this method directly. Instead, override _submit_proper().
		
		Note the order of submission; which is method -> program -> calculation.
		
		The calculation will occur during this method and will have completed before this method returns (for blocking Method_targets).
		"""
		self.method.submit_proper()
		self._submit_proper()
		
	def _submit_proper(self):
		"""
		Step 3/4 of the submission process, this method is called to perform submission.
		
		For Program_targets, the calculation should be performed here.
		
		This default implementation does nothing.
		"""
		pass
		
	def submit_post(self):
		"""
		Step 4/4 of the submission process, this method is called after submission.
		
		Inheriting classes should avoid overriding this method directly. Instead, override _submit_post().
		
		Note the order of submission; which is method -> program -> calculation.
		"""
		self.method.submit_post()
		self._submit_post()
		
	def _submit_post(self):
		"""
		Step 4/4 of the submission process, this method is called after submission.
		"""
		# First, (try) and load our results.
		# We'll actually load a report object because it is the same as a Result_set but with some extra methods which we might need to write a report. Saves us loading results twice.
		# We need to know whether the calculation was successful or not, so we make no effort to catch exceptions here.
		self.result = PDF_report.from_calculation_files(
			date = self.end_date,
			duration = self.duration,
			gaussian_log_file = self.calc_output_file_path,
			prog_version = silico.version,
			options = self.silico_options
			)
		
		# If we've been asked to write result files, do so.
		try:
			if self.calculation.write_summary:
				self.write_summary_files()
		except Exception:
			getLogger(silico.logger_name).warning("Failed to write calculation result summary files", exc_info = True)
			
		# Similarly, if we've been asked to write a report, do that.
		try:
			if self.calculation.write_report:
				self.write_report_files()
		except Exception:
			getLogger(silico.logger_name).warning("Failed to write calculation report", exc_info = True)
		
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
		Text_summary_group_extractor(ignore = True).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.name +".summary"))
		
		# Atoms.
		Long_CSV_group_extractor(Atoms_long_extractor(ignore = True)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.name + ".atoms.csv"))
		
		# Alpha. We'll use a different file name depending on whether we are restricted or unrestricted.
		Long_CSV_group_extractor(Orbitals_long_extractor(ignore = True)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.name + ".{}.csv".format("orbitals" if self.result.metadata.orbital_spin_type == "restricted" else "alpha")))
		# And beta (which can only be for unrestricted.
		Long_CSV_group_extractor(Beta_long_extractor(ignore = True)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.name + ".beta.csv"))
		
		# Energies.
		Long_CSV_group_extractor(SCF_long_extractor(ignore = True)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.name + ".SCF.csv"))
		Long_CSV_group_extractor(MP_long_extractor(ignore = True)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.name + ".MP.csv"))
		Long_CSV_group_extractor(CC_long_extractor(ignore = True)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.name + ".CC.csv"))
		
		# Excited states.
		Long_CSV_group_extractor(Excited_state_long_extractor(ignore = True)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.name + ".ES.csv"))
		Long_CSV_group_extractor(Excited_state_transitions_long_extractor(ignore = True)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.name + ".transitions.csv"))
		Long_CSV_group_extractor(TDM_long_extractor(ignore = True)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.name + ".TDM.csv"))
		
		# And vibrations.
		Long_CSV_group_extractor(Vibrations_long_extractor(ignore = True)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.name + ".vibrations.csv"))
		
	def write_report_files(self):
		"""
		Write report files (like with creport) from this calculation.
		"""
		# The full report.
		self.result.write(Path(self.method.calc_dir.report_directory, self.calculation.name + ".pdf"))
		# And atoms.
		self.result.write(Path(self.method.calc_dir.report_directory, self.calculation.name + ".atoms.pdf"), report_type = "atoms")
		
		
		
		