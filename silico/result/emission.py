from silico.result.excited_states import Excited_state, Excited_state_list
from silico.image.excited_states import Excited_states_diagram_maker
from pathlib import Path
from silico.image.spectroscopy import Absorption_graph_maker
from silico.result import Result_object
from silico.exception.base import Silico_exception


class Relaxed_excited_state(Excited_state):
	"""
	Class for representing an emission energy from an excited state to a ground state.
	
	Note that while emission and absorption can be approximated as the reverse of each other, emission and absorption strictly have different energies.
	Excited states, as represented by the Excited_state class, are for absorption energies.
	This class is for emission energies, which are typically lower in energy (because the excited state has relaxed). 
	"""
	
	def __init__(self, ground_state_result, excited_state_result, transition_type, excited_state = None):
		"""
		Constructor for Relaxed_excited_state objects.
		
		:param ground_state_result: A Result_set object representing the ground state.
		:param excited_state_result: A Result_set object representing the excited state.
		:param excited_state: An optional Excited_state object. This is required (for example) in time dependent DFT where the total energy of excited_state_result is the ground state energy (at the excited state geometry).
		:param transition_type:  A string describing the type of transition, either 'adiabatic' (GS and ES relaxed) or 'vertical' (ES relaxed, GS @ ES geom).
		"""
		# We don't call the Excited_state constructor (yet) because it handles energy differently.
		Result_object.__init__(self)
		self.ground_state_result = ground_state_result
		self.excited_state_result = excited_state_result
		self.excited_state = excited_state
			
		if self.excited_state is not None:
			# If we have an excited state we can inherit certain properties from it.
			self.level = self.excited_state.level
			self.multiplicity_level = excited_state.multiplicity_level
			
			self.oscillator_strength = self.excited_state.oscillator_strength
		else:
			# For now we assume this is the lowest possible excited state (may change in future).
			self.level = 1
			self.multiplicity_level = 1
			
			# We don't have a concept of oscillator strength (yet?).
			self.oscillator_strength = None
			
		self.transition_type = transition_type
		
	@classmethod
	def from_results(self, main_result, *, ground_state_result = None, excited_state_result, transition_type, excited_state = None):
		"""
		A more intelligent constructor for Relaxed_excited_state object.
		
		:param main_result: A Result_set object containing main calculation results.
		:param ground_state_result: A Result_set object representing the ground state.
		:param excited_state_result: A Result_set object representing the excited state.
		:param transition_type:  A string describing the type of transition, either 'adiabatic' (GS and ES relaxed) or 'vertical' (ES relaxed, GS @ ES geom).
		:param excited_state: An optional Excited_state object. This is required (for example) in time dependent DFT where the total energy of excited_state_result is the ground state energy (at the excited state geometry). If excited_state is not an Excited_state object (and isn't None), it is passed as criteria to Excited_state_list.get_state(). 
		"""
		# Decide on what to use for our ground state
		if ground_state_result is None:
			# If our excited results has TD excited states and we are a vertical transition, then it can be our ground state.
			if len(excited_state_result.excited_states) > 0 and transition_type == "vertical":
				ground_state_result = excited_state_result
			elif transition_type == "adiabatic":
				# Check we are an Opt, and complain if not.
				if 'Optimisation' in main_result.metadata.calculations:
					# No TD excited states, so we'll use the main result object.
					ground_state_result = main_result
				else:
					# No explicit ground given and out main result is not an Opt, so it probably isn't suitable.			
					raise Silico_exception("Unable to determine ground state in adiabatic emission; no explicit ground state given and this calculation is not an optimisation")
			else:
				if 'Single Point' in main_result.metadata.calculations:
					ground_state_result = main_result
				else:
					# Vertical transition via the unrestricted triplet method; we could (in theory) be the ground state, but might also be optimised ground (which would give adiabatic, not vertical).
					# Sadly, there's no real way of knowing for sure (even if we are a single point, it could be a single point at the optimised geom).
					# The best we do is guess; if we are a singlepoint, then we are OK as the ground.
					# If we are not; then we'll stop here. The user can always explicitly set this object as the ground if it is correct.
					# TODO: Maybe just convert to a warning?
					raise Silico_exception("Unable to determine ground state in vertical emission; no explicit ground state given and this calculation is not a single point")
			
		# Now decide which excited state to use.
		excited_state = excited_state if excited_state is None or isinstance(excited_state, Excited_state) else excited_state_result.excited_states.get_state(excited_state)
		
		if excited_state is None and len(excited_state_result.excited_states) > 0:
			excited_state = excited_state_result.excited_states[0]
		
		return self(
			ground_state_result,
			excited_state_result,
			transition_type,
			excited_state)
	
	@property
	def excited_energy(self):
		"""
		The total energy of the excited state in this transition.
		"""
		# Start with our excited_state_result energy.
		excited_energy = self.excited_state_result.energy
		
		# Add the excited state energy (if we have it).
		if self.excited_state is not None:
			excited_energy += self.excited_state.energy
			
		return excited_energy
			
	@property
	def ground_energy(self):
		"""
		The total energy of the ground state in this transition.
		"""
		return self.ground_state_result.energy
	
	@property
	def energy(self):
		"""
		The energy of this transition (in eV).
		"""
		# Return the difference between ground and excited.
		return self.excited_energy - self.ground_energy
	
	@property
	def multiplicity(self):
		"""
		This is an alias of excited_multiplicity
		"""
		return self.excited_multiplicity
	
	@property
	def excited_multiplicity(self):
		"""
		The multiplicity (as a number) of the excited state in this emission transition.
		"""
		if self.excited_state is not None:
			return self.excited_state.multiplicity
		else:
			return self.excited_state_result.metadata.system_multiplicity
		
	@property
	def ground_multiplicity(self):
		"""
		The multiplicity (as a number) of the ground state in this emission transition.
		"""
		return self.ground_state_result.metadata.system_multiplicity
	
	@property
	def emission_type(self):
		"""
		The emission type (as a string), either fluorescence or phosphorescence.
		"""
		if self.ground_multiplicity == self.excited_multiplicity:
			return "fluorescence"
		else:
			return "phosphorescence"
		
		
	def set_file_options(self, output_dir, output_name, *, ground_state, **kwargs):
		"""
		Set the options that will be used to create images from this object.
		
		This method will also call set_file_options on the excited states contained in this list.
		
		:param output_dir: A pathlib Path object to the directory within which our files should be created.
		:param output_name: A string that will be used as the start of the file name of the files we create.
		"""
		# First our states diagram.
		self._files['excited_states_diagram'] = Excited_states_diagram_maker.from_image_options(
			Path(output_dir, output_name + ".{}_emission_states.png".format(self.transition_type)),
			excited_states = Excited_state_list([self]),
			ground_state = ground_state,
			**kwargs
		)
		
		# Now emission spectrum.
		self._files['simulated_emission_graph'] = Absorption_graph_maker.from_image_options(
			Path(output_dir, output_name + ".simulated_{}_emission_spectrum.png".format(self.transition_type)),
			excited_states = Excited_state_list([self]),
			**kwargs
			)

	@property
	def excited_states_diagram(self):
		return self.get_file('excited_states_diagram')
	
	@property
	def simulated_emission_graph(self):
		return self.get_file('simulated_emission_graph')
		