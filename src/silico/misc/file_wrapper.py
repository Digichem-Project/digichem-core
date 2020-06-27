import sys

class Multi_file_wrapper():
	"""
	A class that can be used to transparently handle both 'normal' files and stdin/stdout.
	"""
	
	def __init__(self, file, mode = "r", *args, **kwargs):
		# If the filename is the special symbol '-' (a dash), then 'open' stdout as appropriate).
		if file == '-':
			if 'w' in mode or 'a' in mode:
				# We are writing.
				if 'b' in mode:
					# We are a binary format, use a different stdout that accepts binary output.
					object.__setattr__(self, 'file', sys.stdout.buffer)
					#self.file = sys.stdout.buffer
				else:
					object.__setattr__(self, 'file', sys.stdout)
					#self.file = sys.stdout
			elif '+' in mode:
				# We are updating (not allowed).
				raise ValueError("Invalid mode specified '{}', updating is not valid for stdin/stdout".format(mode))
			else:
				# We are reading.
				if 'b' in mode:
					object.__setattr__(self, 'file', sys.stdin.buffer)
					#self.file = sys.stdin.buffer
				else:
					object.__setattr__(self, 'file', sys.stdin)
					#self.file = sys.stdin
			# We don't close stdin/stdout.
			object.__setattr__(self, 'should_close_file', False)
			#self.should_close_file = False
		else:
			# Normal file.
			
			object.__setattr__(self, 'file', open(file, mode, *args, **kwargs))
			object.__setattr__(self, 'should_close_file', True)
			#self.file = open(file, mode, *args, **kwargs)
			#self.should_close_file = True
			
	def __getattr__(self, name):
		"""
		Magic method so we can pretend to be a real file.
		"""
		return getattr(self.file, name)
	
	def __setattr__(self, name, value):
		"""
		Magic method so we can pretend to be a real file.
		"""
		return setattr(self.file, name, value)
		
	def __delattr__(self, name):
		"""
		Magic method so we can pretend to be a real file.
		"""
		if name in self.__dict__:
			object.__delattr__(self, name)
		else:
			return delattr(self.file, name)
			
			
	def close(self):
		if self.should_close_file:
			self.file.close()
		else:
			self.file.flush()
			
	def __enter__(self):
		"""
		Magic function for the 'with' keyword.
		"""
		return self
		
	def __exit__(self, etype, value, traceback):
		"""
		Magic function for the 'with' keyword, called at the end of the block.
		"""
		# Close our file.
		self.close()