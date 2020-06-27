from silico.exception.uncatchable import Submission_paused
import pickle
from pathlib import Path

class Resumable_method():
	"""
	Mixin class for submit_methods that are resumable.
	
	Resumable here means that program execution stops during the submission process. Submission is then 'resumed' in a new process immediately before the calculation proper begins.
	This mechanism is required by several methods, for example SLURM (calculation has to occur from the SLURM node) and SSH (calculation occurs on a different machine entirely).
	
	The 'resume' is achieved via pickle, so you class should be picklable.
	"""
	
	def __init__(self):
		"""
		Constructor for Resumable_method objects.
		
		:param output: Path to the file that we will write to resume from.
		"""
		# A flag keeping track of which side of the pickle reload we are.
		self._resumed = False
		
	def _submit_pre(self):
		"""
		
		"""
		# If we've already resumed, then there's nothing for us to do.
		if not self._resumed:
			# Pause here.
			self.pause()
			
			# Use a special 'exception' to prevent normal submission.
			raise Submission_paused()
		
	@property
	def resume_file_path(self):
		"""
		Path to our pickled resume file. 
		"""
		return Path(self.calc_dir.input_directory, "silico.resume.pickle")
		
	def pause(self):
		"""
		The first part of the 'resume' mechanism, pause() should set-up the class for resuming later. Typically this involves writing to a pickle file.
		
		#This default implementation is a helper function for inheriting classes, calling it will pickle this object to the path given. Normal implementations will take no arguments.
		
		#:param output: The path to pickle to. Note that this argument is required, it is set as optional so a more descriptive error message can be given.
		"""
		#if output is None:
		#	raise Silico_exception("default Resumable_method.pause() implementation called without output argument; if you are writing your own Resumable_method class, you should write your own implementation of pause()")
		
		with open(self.resume_file_path, "bw") as pickle_file:
			pickle.dump(self, pickle_file)
		
	def resume(self):
		"""
		The second part of the 'resume' mechanism, resume() is called after the pick file has been re-loaded. Execution should continue from here.
		"""
		self._resumed = True
		
		# Continue submitting.
		self.program.calculation.submit_pre()
		self.program.calculation.submit_proper()
		
		# We now rest our flag because we might re-use the same SLURM method for the next calc, which hasn't resumed yet.
		self._resumed = False
		self.program.calculation.submit_post()
		
		
		
		
		