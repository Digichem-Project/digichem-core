from silico.submit.method import Method_target
from silico.submit.method.resume import Resumable_method
import subprocess
from subprocess import CalledProcessError
from silico.exception.base import Submission_error

class Series(Method_target):
	"""
	Implementation to allow local submission.
	
	'Local' here indicates that the calculation is to be performed right-here, right-now, in the current-process.
	Local is useful if an external scheduling manager is being used etc.
	
	Note this this method is blocking (not resumable), and will carry out each calculation in series, one after another.
	"""
	
	CLASS_HANDLE = ("series",)
	
	
class Parallel(Resumable_method):
	"""
	Implementation to allow local submission.
	
	'Local' here indicates that the calculation is to be performed right-here, right-now.
	Local is useful if an external scheduling manager is being used etc.
	
	Unlike Series, Parallel will create a separate subprocess for each calculation.
	"""
	
	CLASS_HANDLE = ("parallel",)
	
	
	############################
	# Class creation mechanism #
	############################
	
	class _actual(Resumable_method._actual):
		"""
		Inner class for Parallel.
		"""
		
		def post(self):
			"""
			Submission method called after pausing.
			"""
			# Create a new silico instance.
			try:
				# Open output log.
				# We don't wrap this in a with statement; Popen will close the file object on its own?
				stdout = open(self.calc_dir.log_file, "wt")
				
				subprocess.Popen(
					["silico", "resume", self.resume_file_path],
					stdout = stdout,
					stderr = subprocess.STDOUT,
					start_new_session = True,
					universal_newlines = True
				)
			except CalledProcessError as e:
				# Something went wrong.
				e.__context__ = None
				
				raise Submission_error(self, "failed to execute silico child:\n{}".format(e.stdout)) from e
			except FileNotFoundError as e:
				# Couldn't find sbatch.
				e.__context__ = None
				raise Submission_error(self, "unable to locate silico executable '{}'") from e
			except Exception as e:
				raise e from None
			
			
		