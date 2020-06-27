from silico.submit.method import Method_target
from silico.exception.base import Configurable_target_exception,\
	Submission_error, Silico_exception
from mako.lookup import TemplateLookup
from pathlib import Path
import silico
import os
import stat
from silico.submit import Memory
from silico.submit.method.resume import Resumable_method
from silico.exception.uncatchable import Submission_paused
import subprocess
from subprocess import CalledProcessError

class SLURM(Method_target, Resumable_method):
	"""
	Implementation to allow submission to SLURM, a popular scheduling system.
	"""
	
	CLASS_HANDLE = "SLURM"
	
	# The name of the script which we pass to sbatch.
	SBATCH_SCRIPT_NAME = "sbatch.submit"
	
	def __init__(
			self,
			*args,
			partition = None,
			time = None,
			num_tasks = None,
			CPUs_per_task = None,
			mem_per_CPU = None,
			common_directory = None,
			sbatch_command = None,
			sinfo_command = None,
			options = None,
			**kwargs
		):
		super().__init__(*args, **kwargs)
		Resumable_method.__init__(self)
		
		self.partition = partition
		self.time = time
		self.num_tasks = num_tasks if num_tasks is not None else 1
		self.CPUs_per_task = CPUs_per_task if CPUs_per_task is not None else "auto"
		self.mem_per_CPU = mem_per_CPU if mem_per_CPU is not None else "auto"
		self.options = options if options is not None else {}
		if common_directory is False:
			raise Configurable_target_exception(self, "common_directory = False is currently not supported")
		self.common_directory = common_directory if common_directory is not None else True
		self.sbatch_command = sbatch_command if sbatch_command is not None else "sbatch"
		self.sinfo_command = sinfo_command if sinfo_command is not None else "sinfo"
		
	def get_num_free_nodes(self):
		"""
		Get the current number of free nodes for this partition.
		
		:raises Silico_exception: If the number of nodes could not be determined.
		:return: The number of nodes.
		"""
		# The signature we'll use to call sinfo.
		sig = [
			self.sinfo_command,
			"-p", self.partition,
			"-t", "IDLE",
			"-o", "%D",
			"-h"
		]

		try:
			# Call sinfo.
			done = subprocess.run(
				sig,
				stdout = subprocess.PIPE,
				stderr = subprocess.PIPE,
				universal_newlines = True,
				check = True,
			)
			
			# The number of nodes will be our output (unless there were none, in which case nothing is outputed).
			num_nodes = int(done.stdout) if done.stdout != "" else 0
		except Exception:
			raise Silico_exception("Failed to retrieve sinfo for partition {}".format(self.partition))
		
		return num_nodes
	
	def get_CPU_info(self):
		"""
		Get information on the current allocation of CPUs in this partition.
		
		:return: A tuple of integers specifying the current CPU allocation in the format (allocated, idle, other, total).
		"""
		# The signature we'll use to call sinfo.
		sig = [
			self.sinfo_command,
			"-p", self.partition,
			"-o", "%C",
			"-h"
		]

		try:
			# Call sinfo.
			done = subprocess.run(
				sig,
				stdout = subprocess.PIPE,
				stderr = subprocess.PIPE,
				universal_newlines = True,
				check = True,
			)
			
			# Our output is a string of the form alloc/idle/other/total (unless the partition isn't valid, in which case we get nothing).
			cpu_info = done.stdout.split("/")
			
			# This will trigger an error if the format is weird, which is what we want.
			cpu_info = (int(cpu_info[0]),int(cpu_info[1]),int(cpu_info[2]),int(cpu_info[3]))
		except Exception:
			raise Silico_exception("Failed to retrieve sinfo for partition {}".format(self.partition))
		
		return cpu_info
	
	@property
	def status(self):
		"""
		This method is called to get 'status' about this method.
		
		For SLURM, we use sinfo to get the number of free nodes for this partition.
		
		:raises Silico_exception: If the number of nodes could not be determined.
		:return: Info string.
		"""
		cpu_info = self.get_CPU_info()
		return "{} idle nodes, {} ({:0.0f}%) idle CPUs".format(self.get_num_free_nodes(), cpu_info[1], (cpu_info[1]/cpu_info[3])*100 if cpu_info[3] != 0 else 0)
		#return "{} free nodes, {} ({:0.0f}%) free CPUs".format(5, 224, (224/2848)*100)
		
	
	@property
	def CPUs_per_task(self):
		"""
		Get the number of CPUs to assign for this calculation.
		
		This property will resolve 'auto' to an actual number of processors, use _CPUs_per_task if you do not want this behaviour.
		"""
		if self._CPUs_per_task == "auto":
			return self.program.calculation.num_CPUs if self.program.calculation.num_CPUs is not None else 1
		else:
			return self._CPUs_per_task
		
	@CPUs_per_task.setter
	def CPUs_per_task(self, value):
		"""
		Set the number of CPUs to assign for this calculation.
		"""
		self._CPUs_per_task = value
		
	@property
	def mem_per_CPU(self):
		"""
		Get the amount of memory to assign (per CPU).
		
		This property will resolve 'auto' to an actual amount of memory, use _mem_per_CPU if you do not want this behaviour.
		"""
		if self._mem_per_CPU == "auto":
			return SLURM_memory(round(float(self.program.calculation.memory)) / self.CPUs_per_task)
		else:
			return SLURM_memory(float(self._mem_per_CPU))
	
	@mem_per_CPU.setter
	def mem_per_CPU(self, value):
		"""
		Set the amount of memory to assign (per CPU).
		"""
		self._mem_per_CPU = value
	
	@property
	def unique_name(self):
		"""
		Get a name that is unique for this calculation instance.
		
		SLURM tries to get a unique name based on our allocated SLURM ID, but if we have not yet been submitted to SLURM this method will fallback to a different method.
		"""
		if getattr(self, "_unique_name", None) is None:
			try:
				self._unique_name = os.environ['SLURM_JOB_ID']
			except KeyError:
				# SLURM_JOB_ID isn't set, means we haven't been submitted to SLURM yet.
				# Fallback.
				self._unique_name = super().unique_name
		
		return self._unique_name			
		
	def _submit_init(self, program):
		"""
		Step 1/4 of the submission process, this method is called to set-up submission.
				
		:param program: A Program_target object that is going to be submitted. This Program_target object will have completed submit_init() before this method is called.
		"""
		# Call our parent; this creates our directory structure for us.
		super()._submit_init(program)

	def write_sbatch_script(self):
		"""
		Write the control script which is passed to sbatch to file.
		
		:param path: The path to where the file should be written (this should point to a directory; the filename is appended automatically).
		"""
		# Get and load our template.
		template_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/slurm.mako").render_unicode(SLURM_target = self)
		
		with open(self.sbatch_script_path, "wt") as sbatch_file:
			sbatch_file.write(template_body)
			
		# Make it executable.
		os.chmod(self.sbatch_script_path, os.stat(self.sbatch_script_path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
	
	@property
	def sbatch_script_path(self):
		"""
		Path to the bash script file that we will write and pass to sbatch.
		"""
		return Path(self.calc_dir.input_directory, self.SBATCH_SCRIPT_NAME)
	
	def _submit_pre(self):
		"""
		Method called at the start of submission.
		
		SLURM is a resumable method; this method will get called twice (automatically) during the submission process. This method can raise Submission_paused exceptions as part of this proces, you should not go out of you way to catch these exceptions unless you know what you are doing (and you should almost certainly be re-raising once you are done).
		"""
		# First, call our parent (currently does nothing).
		super()._submit_pre()
		
		# If we have not yet resumed, create our directory structure.
		# We need to do this before the resume because this is where we'll write out SLURM batch file.
		if not self._resumed:
			try:
				self.calc_dir.create_structure(True)
			except Exception:
				raise Submission_error(self.program.calculation, "could not create directory structure; try setting a different output directory ('-o')")
					
			# Write the control file.
			self.write_sbatch_script()
		
		# Now call our Resumable_method
		try:
			Resumable_method._submit_pre(self)
		except Submission_paused as paused:
			# We are before the resume, call sbatch.
			try:
				subprocess.run(
					[self.sbatch_command, self.sbatch_script_path],
					stdout = subprocess.PIPE,
					stderr = subprocess.STDOUT,
					check = True,
					universal_newlines = True,
				)
			except CalledProcessError as e:
				# Something went wrong.
				e.__context__ = None
				raise Submission_error(self, "{} did not exit successfully:\n{}".format(self.sbatch_command, e.stdout)) from e
			except FileNotFoundError as e:
				# Couldn't find sbatch.
				e.__context__ = None
				raise Submission_error(self, "unable to locate sbatch executable '{}'".format(self.sbatch_command)) from e
			except Exception as e:
				raise e from None
			
			# Continue exiting.
			raise paused
		
		# If we get this far, then we have resumed and can continue as normal.
		
		# This is not actually true; slurm does not change working dir.
		# When we resume, our working directory will have changed (to inside the input directory in fact). This means our paths are no longer correct.
		# Update our directory object.
		
		# DEBUGGING ONLY
		#print("DEBUGGING BREAKPOINT HANDLE")
		#os.chdir("/home/oliver/ownCloud/Chemistry/St. Andrews PhD/Test Molecules/Benzene/Opt Freq PBE0_6-31G(d,p) UltraFine/Input")

		
		#self.calc_dir.molecule_directory.path = "../../"
		
		
		
	
class SLURM_memory(Memory):
	
	# SLURM has its own set of units...
	UNITS = {
		'T': 1000000000000,
		'G': 1000000000,
		'M': 1000000,
		'K': 1000,
		}
			
			