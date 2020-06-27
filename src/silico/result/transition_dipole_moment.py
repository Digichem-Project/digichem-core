from silico.result.dipole_moment import Dipole_moment
from silico.exception.base import Result_unavailable_error
from silico.file.cube import Cube_maker
from pathlib import Path
from silico.image.vmd import Dipole_image_maker

class Transition_dipole_moment(Dipole_moment):
	"""
	Class that represents a transition dipole moment.
	
	Note that this class is almost identical to the standard (permanent) dipole moment, but the way in which we fetch the data is different.
	"""
	
	# The text we look for that tells us where the transition dipole data is printed.
	t_dipole_section_header = " Ground to excited state transition electric dipole moments (Au):\n"
	
	def __init__(self, state_level, *args, **kwargs):
		"""
		Transition_dipole_moment constructor.
		
		:param state: The excited state level this dipole transitions to.
		"""
		super().__init__(*args, **kwargs)
		
		# An int describing the excited state we belong to. This is always available.
		self.state_level = state_level
		
		# The excited state object that we belong to. This may be None. We can't set this here because then both the Excited_state and Transition_dipole_moment constructors would depend on each other.
		self.excited_state = None
		
		# Save a name describing which dipole we are (permanent vs transition etc).
		self.dipole_type = "transition"
	
	@property
	def name(self):
		"""
		Name that describes this dipole moment.
		"""
		return "{} TDM".format(self.excited_state.state_symbol)
	
	def set_excited_state(self, excited_state):
		"""
		Set the excited state object that we belong to.
		"""
		# Not sure what we should do if the given excited_state has a different level to what we expect.
		if excited_state.level != self.state_level:
			pass
		self.excited_state = excited_state
		
	def set_file_options(self, output_dir, output_name, cube_file = None, **kwargs):
		"""
		Set the options that will be used to create images from this object.
		
		:param output_dir: A pathlib Path object to the directory within which our files should be created.
		:param output_name: A string that will be used as the start of the file name of the files we create.
		"""
		# Work out what we'll name our files.
		file_name = "{}_dipole".format(self.excited_state.state_symbol)
		sub_dir_name = "{} Transition Dipole Moment".format(self.excited_state.state_symbol)
		
		# Get ourselves a cube file maker if we need one.
		if cube_file is None:
			# We'll just use the HOMO to get our cube, as it almost certainly should exist.
			cube_file = Cube_maker.from_image_options(Path(output_dir, sub_dir_name, output_name + ".{}.cube".format(file_name)), cubegen_type = "MO", orbital = "HOMO", **kwargs)
			
		# Get our image.
		self._files['dipole_image'] = Dipole_image_maker.from_image_options(Path(output_dir, sub_dir_name, output_name + ".{}.jpg".format(file_name)), cube_file = cube_file, dipole_moment = self, **kwargs)
	
	def cleanup_intermediate_files(self):
		"""
		Remove any intermediate files that may have been created by this object.
		"""
		# Remove our cube file.
		super().cleanup_intermediate_files('cube_file')
	
	@classmethod
	def list_from_log(self, log_file_path, atoms = None):	
		"""
		Construct a list of transition Dipole_moment objects from a Gaussian logfile.
		
		Currently cclib does not support extracting transition dipole moment data from Gaussian files, so we'll have to fetch it ourselves.
		:param log_file_path: Path to a logfile to read data from.
		:param atoms: An Atom_list object that is passed to the Dipole_moment's constructor. Used for getting alignment info about the dipole etc.
		:return: The list of transition dipole moment objects. This list will be empty (len == 0) if no data is available.
		"""
		# First, try and get our list of raw dipole moment data.
		try:
			# Get our raw list.
			raw_dipole_list = self.get_t_dipole_section(log_file_path)
			
			# Sort our list (we probably don't need to do this, but it doesn't hurt.)
			raw_dipole_list.sort(key = lambda dipole: dipole['state_level'])
			
			# Build our object list.
			return [self(state_level = raw_dipole['state_level'], origin_coords = raw_dipole['origin_coords'], vector_coords = raw_dipole['vector_coords'], atoms = atoms) for raw_dipole in raw_dipole_list]
		
		except Result_unavailable_error:
			return []
		
	@classmethod
	def get_t_dipole_section(self, log_file_path):
		"""
		Get and return the whole transition dipole moment section without further formating.
		
		:param logfile: Path to a logfile to read data from.
		:return: The t dipole section as an unformatted string
		"""
		transition_dipole_groups = []
		
		# Open our file.
		with open(log_file_path, "rt") as log_file:
			# Loop through our lines, looking for a specific line of text.
			#found = False
			for log_file_line in log_file:
				# Look for our key string.
				# It might be better to look for a substring rather than the whole string? Need to benchmark.
				if log_file_line == self.t_dipole_section_header:
					# We found our header.
					
					# The next line should be the table header, so we can skip it.
					log_file.readline()
					
					# Now loop through, splitting on white space.
					transition_dipoles = []
					for log_file_line in log_file:
						# Split on whitespace.
						split_line = log_file_line.split()
						
						try:
							# Add into our list.
							transition_dipoles.append({
								'state_level': int(split_line[0]),
								'origin_coords': (0.0, 0.0, 0.0),
								'vector_coords': (
									self.au_to_debye(float(split_line[1])),
									self.au_to_debye(float(split_line[2])),
									self.au_to_debye(float(split_line[3]))
								)
								})
						except (ValueError, IndexError):
							# No more data.
							break
					
					# Add this set of transition dipoles to our big list.
					transition_dipole_groups.append(transition_dipoles)
		
		# Return the last t-dipole moment (unless we have none, in which case we get upset).
		if len(transition_dipole_groups) == 0:
			raise Result_unavailable_error("Transition dipole moment", "There is no transition dipole moment data in '{}'".format(log_file_path))
		else:
			return transition_dipole_groups[-1]
	
	@classmethod
	def au_to_debye(self, au):
		"""
		Convert a dipole moment in au to debye.
		"""
		return au * 2.541746473
	
	@classmethod
	def bohr_to_angstrom(self, bohr_distance):
		"""
		Convert a length in bohr to angstrom.
		"""
		return bohr_distance * 0.529177210903
		
			