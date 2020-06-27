from pathlib import Path
import os

# Based on https://stackoverflow.com/questions/431684/how-do-i-change-the-working-directory-in-python/24176022#24176022
class cd:
	"""
	Context manager for temporarily changing the working directory.
	"""
	
	def __init__(self, directory):
		"""
		Constructor for cd.
		
		:param directory: Path to the directory to change to.
		"""
		self.new_directory = Path(directory)
		self.old_directory = None
		
	def __enter__(self):
		# Save our current working directory.
		self.old_directory = os.getcwd()
		
		# And change to the new directory.
		os.chdir(str(self.new_directory))
		
	def __exit__(self, exc_type, exc_value, traceback):
		# Restore our old directory.
		os.chdir(self.old_directory)