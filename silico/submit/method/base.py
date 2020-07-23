from silico.submit import Configurable_target
from silico.submit.structure.directory import Calculation_directory
from uuid import uuid4
from silico.config.configurable.option import Option

class Method_target(Configurable_target):
	"""
	Top-level class for classes that implement a method of submitting calculations.
	
	'Methods' define how and where a calculation is run, but don't manage the calculation software itself (see silico.submit.program for those definitions).
	Methods, for example, handle submission to a scheduling software (SLURM, TORQUE etc), running as a daemon, submission to a networked server etc.
	"""
	
	CLASS_HANDLE = ("method",)
	
	# Configurable Options.
	warning = Option(help = "A warning message to display when this method is chosen", default = None, type = str)
	
	@property
	def unique_name(self):
		"""
		Get a name that is unique for this calculation instance.
		
		Some methods may provide their own get_unique_name() methods; this default implementation returns a random string with a very low collision chance.
		"""
		if getattr(self, "_unique_name", None) is None:
			self._unique_name = uuid4().hex
			
		return self._unique_name
	
	@property
	def silico_options(self):
		"""
		Get the global silico options dictionary (this is actually found under our Calculation_target).
		"""
		return self.program.calculation.silico_options
	
	@property
	def status(self):
		"""
		This method is called to get 'status' about this method.
		
		Status is always a string, but what it contains depends entirely on the Method_target.
		It is used, for example, to report the number of free nodes for SLURM.
		
		This default implementation raises NotImplementedError.
		"""
		raise NotImplementedError
	
	def submit_init(self, program):
		"""
		Step 1/4 of the submission process, this method is called to set-up submission.
		Importantly, submit_init() will return before any resumable submission methods have resumed, meaning the environment during submit_init() may not reflect the final submission environment. You should typically avoid writing files that you will need later, because they may not be available in later submission methods.
		
		Inheriting classes should avoid overriding this method directly. Instead, override _submit_init().
		
		It is important to note that the normal submission order is reversed for submit_init(); the order is calculation -> program -> method.
		
		:param program: A Program_target object that is going to be submitted. This Program_target object will have completed submit_init() before this method is called.
		"""
		self._submit_init(program)
		
	def _submit_init(self, program):
		"""
		Step 1/4 of the submission process, this method is called to set-up submission.
		
		Inheriting classes should override this method to perform init.
		
		This default implementation saves program to an attribute of the same name. Call super()._submit_init() in your implementation if you want this behaviour.
		If you do not call super()._submit_init(), know that several other classes will expect the program attribute to exist. 
		
		:param program: A Program_target object that is going to be submitted. This Program_target object will have completed submit_init() before this method is called.
		"""
		self.program = program
		
		# We'll set our output directory here, but we won't make it yet (so resumable methods can decide when to create it themselves.
		self.calc_dir = Calculation_directory.from_calculation(self.program.calculation)

	def submit_pre(self):
		"""
		Step 2/4 of the submission process, this method is called before submission begins.
		
		Inheriting classes should avoid overriding this method directly. Instead, override _submit_pre().
		
		Note the order of submission; which is method -> program -> calculation.
		"""
		self._submit_pre()
	
	def _submit_pre(self):
		"""
		Step 2/4 of the submission process, this method is called before submission begins.
		
		This default implementation does nothing.
		"""
		pass
	
	def submit_proper(self):
		"""
		Step 3/4 of the submission process, this method is called to perform submission.
		
		Inheriting classes should avoid overriding this method directly. Instead, override _submit_proper().
		
		Note the order of submission; which is method -> program -> calculation.
		"""
		self._submit_proper()
	
	def _submit_proper(self):
		"""
		Step 3/4 of the submission process, this method is called to perform submission.
		
		This default implementation does nothing.
		"""
		pass
	
	def submit_post(self):
		"""
		Step 4/4 of the submission process, this method is called after submission.
		
		Inheriting classes should avoid overriding this method directly. Instead, override _submit_post().
		
		Note the order of submission; which is method -> program -> calculation.
		"""
		self._submit_post()
	
	def _submit_post(self):
		"""
		Step 4/4 of the submission process, this method is called after submission.
		"""
		pass
	
	