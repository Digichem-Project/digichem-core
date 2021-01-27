# General imports.
from pathlib import Path

# Silico imports.
from silico.result.dipole_moment import Dipole_moment
from silico.file.cube import Fchk_to_cube
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
	
	@classmethod
	def list_from_parser(self, parser):
		"""
		Get a list of TDMs from an output file parser.
		
		:param parser: An output file parser.
		:return: A list of TDM objects. An empty list will be returned if no TDM data is available.
		"""
		try:
			return [
				self(state_level = index +1,
				origin_coords = (0,0,0),
				vector_coords = (parser.au_to_debye(tdm[0]), parser.au_to_debye(tdm[1]), parser.au_to_debye(tdm[2])),
				atoms = parser.results.alignment) for index, tdm in enumerate(parser.data.etdips)]
		except AttributeError:
			return []
	
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
			cube_file = Fchk_to_cube.from_image_options(Path(output_dir, sub_dir_name, output_name + ".{}.cube".format(file_name)), cubegen_type = "MO", orbital = "HOMO", **kwargs)
			
		# Get our image.
		self._files['dipole_image'] = Dipole_image_maker.from_image_options(Path(output_dir, sub_dir_name, output_name + ".{}.jpg".format(file_name)), cube_file = cube_file, dipole_moment = self, **kwargs)
	
	def cleanup_intermediate_files(self):
		"""
		Remove any intermediate files that may have been created by this object.
		"""
		# Remove our cube file.
		super().cleanup_intermediate_files('cube_file')
	
		
			