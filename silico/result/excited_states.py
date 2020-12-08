from itertools import filterfalse
from silico.result import Result_container
from silico.result.base import Result_object
import numpy
from silico.exception.base import Result_unavailable_error
from pathlib import Path
from silico.image.excited_states import Excited_states_diagram_maker
from silico.image.spectroscopy import Absorption_graph_maker
import colour
import math
# from logging import getLogger
# import silico

class Excited_state_list(Result_container):
	"""
	Class for representing a group of excited states.
	"""
	
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		
	@property
	def singlet_triplet_energy(self):
		"""
		Get ΔE st; the energy difference between the T1 and S1 excited states.
		
		The returned value is positive if S1 is higher in energy than T1 (the usual case), negative otherwise.
		
		:return: The singlet-triplet splitting energy (in eV).
		"""
		# Get our two states of interest.
		S1 = self.get_state("S(1)")
		T1 = self.get_state("T(1)")
		# Calculate their difference.))
		return float(S1) - float(T1)
		
	def get_state(self, criteria = None, *, state_symbol = None, level = None, mult_level = None):
		"""
		Retrieve a particular excited state from its state symbol (S(1), T(1) etc.).
		
		For multiplicities from 1 -> 4, single capital letters are used ('S', 'D', 'T' or 'Q' respectively). For higher multiplicities, the integer representation is used ('5', '6', '7' and so on).
		The multiplicity level follows in brackets, eg ('S(1)') refers to the first (lowest energy) singlet state.
		
		See Excited_state.state_symbol() for a full description of possible symbols.
		
		:raises ValueError: If the requested symbol could not be found in this list.
		:param criteria: Automatically determine which criteria to search by.
		:param state_symbol: The symbol (a string) to retrieve.
		:param level: The level (an int or string that looks like an int) to retrieve.
		:param mult_level: The multiplicity and multiplicity level (as a tuple, eg '(1,1)' for S(1), '(3,2)' for T(2)) to retrieve.
		:return: The requested excited state object.
		"""
		if criteria is not None:
			if isinstance(criteria, tuple):
				mult_level = criteria
			elif criteria.isdigit() or isinstance(criteria, int):
				level = int(criteria)
			else:
				state_symbol = criteria
		
		# Now get our search func.
		if state_symbol is not None:
			filter_func = lambda state: state.state_symbol != state_symbol
		elif level is not None:
			filter_func = lambda state: state.level != level
		elif mult_level is not None:
			filter_func = lambda state: state.multiplicity != mult_level[0] or state.multiplicity_level != mult_level[1]
		else:
			raise ValueError("Missing criteria to search by; specify one of 'criteria', 'state_symbol', 'level' or 'mult_level'")
		
		try:
			return next(filterfalse(filter_func, self))
		except Exception:
			raise Result_unavailable_error("Excited state", "could not find excited state with symbol = '{}'".format(state_symbol))
		
	
	def set_file_options(self, output_dir, output_name, *, ground_state, **kwargs):
		"""
		Set the options that will be used to create images from this object.
		
		This method will also call set_file_options on the excited states contained in this list.
		
		:param output_dir: A pathlib Path object to the directory within which our files should be created.
		:param output_name: A string that will be used as the start of the file name of the files we create.
		"""
		# First our states diagram.
		self._files['excited_states_diagram'] = Excited_states_diagram_maker.from_image_options(
			Path(output_dir, output_name + ".excited_states.png"),
			excited_states = self,
			ground_state = ground_state,
			**kwargs
		)
		
		# Then our simulated absorption graph.
		self._files['simulated_absorption_graph'] = Absorption_graph_maker.from_image_options(
			Path(output_dir, output_name + ".simulated_absorption_spectrum.png"),
			excited_states = self,
			**kwargs
		)
		
		# Then set our excited states.
		for state in self:
			state.set_file_options(output_dir, output_name, **kwargs)
			
	@property
	def excited_states_diagram(self):
		return self.get_file('excited_states_diagram')
	
	@property
	def simulated_absorption_graph(self):
		return self.get_file('simulated_absorption_graph')
	
	def group(self):
		"""
		Group the excited states in this list by multiplicity.
		
		:return: A dictionary of grouped excited states. Each key will correspond to a multiplicity (1, 2, 3 etc) and each value will be an Excited_state_list of states with that multiplicity. 
		"""
		# Dictionary of grouped states.
		grouped_excited_states = {}
		
		for excited_state in self:
			# Group each ES by multiplicity.
			try:
				# Try and append to our list.
				grouped_excited_states[excited_state.multiplicity].append(excited_state)
			except KeyError:
				# No list exists yet.
				grouped_excited_states[excited_state.multiplicity] = type(self)([excited_state])
				
		# Return.
		return grouped_excited_states
		
	@classmethod
	def from_parser(self, parser):
		"""
		Create an Excited_state_list object from an output file parser.
		
		:param parser: An output file parser.
		:return: The populated Excited_state_list object.
		"""
		return self(Excited_state.list_from_parser(parser))
			

class Excited_state_transition(Result_object):
	"""
	Class that represents a transition that contributes to a particular excited state.
	"""
	
	def __init__(self, level, starting_mo, ending_mo, coefficient):
		"""
		Constructor for excited state transitions.
		
		:param level: The 'level' if this transition. The most significant (highest probability) transition has level 1.
		:param starting_mo: The Molecular_orbital object from which this transition begins.
		:param ending_mo: The Molecular_orbital object to which the transition ends.
		:param coefficient: The coefficient of this orbital. Square to get the probability.
		"""
		self.level = level
		self.starting_mo = starting_mo
		self.ending_mo = ending_mo
		self.coefficient = coefficient
		
	@property
	def probability(self):
		"""
		The probability of this transition.
		"""
		return self.coefficient **2
		
	@classmethod
	def list_from_parser(self, parser):
		"""
		Create a list of excited state transitions from an output file parser.
		
		:param parser: An output file parser.
		:param alpha_mo_list: A Molecular_orbital_list object of the MOs of this system.
		:param beta_mo_list: A Molecular_orbital_list object of the beta MOs of this system. This can be left as null for restricted calcs.
		:return: A list of Excited_state_transition objects.
		"""
		try:
			# Create a tuple of our MOs (helps us later).
			MOs = (parser.results.molecular_orbitals, parser.results.beta_orbitals)
			
			# We'll first create an intermediate list of keyword dicts which we'll then sort.
			data_list = [ 
				{'starting_mo': MOs[starting_mo_AB][starting_mo_index], 'ending_mo': MOs[ending_mo_AB][ending_mo_index], 'coefficient': coefficient}
				for (starting_mo_index, starting_mo_AB), (ending_mo_index, ending_mo_AB), coefficient
				in parser.data.etsecs
			]
			
			# Sort by probability/coefficient.
			data_list.sort(key=lambda keywords: keywords['coefficient'], reverse=True)
			
			# Now get a list objects.
			return [self(index+1, keywords['starting_mo'], keywords['ending_mo'], keywords['coefficient']) for index, keywords in enumerate(data_list)]
			
		except IndexError:
			# Probably because one (or both) of our given mo_lists is empty (or too short).
			raise TypeError("Unable to construct excited state transition; transition is to/from an orbital that is not available")
		except AttributeError:
			# No data.
			return []

class Energy_state(Result_object):
	"""
	Class for representing different energy states of the same system.
	"""
	# Some constants.
	speed_of_light = 299792458 # m/s
	plancks_constant = 6.62607004e-34
	electron_volt = 1.602176634e-19 # J
	
	# Colour categories.
	colors = [
		{"max": 400, "name": "UV"},
		{"max": 420, "name": "Violet"},
		{"max": 470, "name": "Blue"},
		{"max": 505, "name": "Cyan"},
		{"max": 555, "name": "Green"},
		{"max": 595, "name": "Yellow"},
		{"max": 625, "name": "Orange"},
		{"max": 740, "name": "Red"},
		{"max": float("inf"), "name": "IR"}	
	]
	
	def __init__(self, level, multiplicity, multiplicity_level, energy):
		"""
		Constructor for Energy_state objects.
		
		:param level: The ordered (by energy) index of this energy state, where 0 is the GS and 1 is the lowest excited state.
		:param multiplicity_level: the ordered (by energy) index of this energy state in terms of states that share the same multiplicity.
		:param multiplicity: The multiplicity of this state as a number (eg, 1 for singlet, 2 for doublet). Fractional multiplicities are also accepted.
		:param energy: The energy of this state in eV. Whether this value is absolute or relative to another state depends on the implementing class.
		"""
		self.level = level
		self.multiplicity = multiplicity
		self.multiplicity_level = multiplicity_level
		self.energy = energy
	
	def __float__(self):
		"""
		Float of this class.
		"""
		return self.energy
	
	@classmethod
	def multiplicity_number_to_string(self, multiplicity):
		if multiplicity == 1:
			return "singlet"
		elif multiplicity == 2:
			return "doublet"
		elif multiplicity == 3:
			return "triplet"
		elif multiplicity == 4:
			return"quartet"
		elif multiplicity % 1 == 0:
			# Multiplicity is an integer, so return as a stringy whole number.
			return str(int(multiplicity))
		else:
			return str(multiplicity)
		
	
	@property
	def multiplicity_string(self):
		return self.multiplicity_number_to_string(self.multiplicity)
	
	@property
	def multiplicity_symbol(self):
		# Get a shorthand symbol if we can.
		if  self.multiplicity == 1:
			return "S"
		elif  self.multiplicity == 2:
			return "D"
		elif  self.multiplicity == 3:
			return "T"
		elif  self.multiplicity == 4:
			return "Q"
		elif self.multiplicity % 1 == 0:
			# Multiplicity is an integer, so return as a stringy whole number.
			return str(int(self.multiplicity))
		else:
			return str(self.multiplicity)
	
	@property
	def state_symbol(self):
		"""
		A short hand notation to identify this excited state.
		
		If the multiplicity is well defined (singlet, doublet, triplet etc), the symbol starts with an appropriate letter (S, D, T etc), otherwise a numeric multiplicity is used. The symbol ends with an integer in brackets, indicating the excited state's level.
		eg, S1 is the first singlet excited state, S2 is the second and so on.
		"""
		return "{}({})".format(self.multiplicity_symbol, self.multiplicity_level)

class Excited_state(Energy_state):
	"""
	Class for representing an excited state.
	"""
	
	def __init__(self, level, multiplicity_level, symmetry, energy, oscillator_strength, transitions, transition_dipole_moment = None):
		"""
		Constructor for excited state objects.
		
		:param level: The 'level' of this excited state (essentially an index), where the lowest state has a level of 1, increasing by 1 for each higher state.
		:param multiplicity_level: The 'level' of this state within the excited states that have the same multiplicity.
		:param symmetry: The symmetry of this excited state; a complicated term that contains our multiplicity among other things.
		:param energy: The energy of this excited state (in eV).
		:param oscillator_strength: The oscillator strength of this transition.
		:param transitions: The singly excited transitions which make up this transition.
		:param transition_dipole_moment: Optional transition dipole of the excited state.
		"""
		super().__init__(level, self.get_multiplicity_from_symmetry(symmetry), multiplicity_level, energy)
		self.symmetry = symmetry
		# Oscillator strength is a dimensionless float that describes the probability of the transition from the reference state (normally ground) to this excited state.
		self.oscillator_strength = oscillator_strength
		# The transitions which contribute to this state.
		self.transitions = transitions
		
		# If we were given a TDM, set its excited state to ourself.
		if transition_dipole_moment is not None:
			transition_dipole_moment.set_excited_state(self)
		self.transition_dipole_moment = transition_dipole_moment
		
	@classmethod
	def get_multiplicity_from_symmetry(self, symmetry_string):
		"""
		Get the multiplicity of of an excited state from its symmetry.
		
		See multiplicity() for a more detailed description of the format of multiplicity strings.
		:param symmetry_string: The symmetry string.
		:return: The multiplicity as a number (possibly an int, possibly a float).
		"""
		# Split the symmetry string on the dash (-) character.
		multiplicity_string = (symmetry_string.split('-', 1))[0]
		# Convert to a number if we have a string.
		if multiplicity_string == "Singlet":
			return 1
		elif multiplicity_string == "Doublet":
			return 2
		elif multiplicity_string == "Triplet":
			return 3
		elif multiplicity_string == "Quartet":
			return 4
		else:
			# Try and cast to float.
			# Float multiplicities don't make sense in the real world, but are possible results from calculations.
			return float(multiplicity_string)
	
	@property
	def wavelength(self):
		"""
		The wavelength that corresponds to the energy of this excited state (in nm).
		
		:raises Result_unavailable_error: If the energy of this state is 0.
		"""
		try:
		# λ = (c * h) / e
		#return ((self.speed_of_light * self.plancks_constant) / (self.energy * self.electron_volt)) * 1000000000
			return self.energy_to_wavelength(self.energy)
		except FloatingPointError:
			# Our energy is zero.
			raise Result_unavailable_error('excited state wavelength', "excited state '{}' energy is 0 eV".format(self.state_symbol))
		
	@property
	def color(self):
		"""
		The 'color' that corresponds to the energy of this excited state (as a string).
		"""
		# Go through our color presets and see which one we match.
		for color in self.colors:
			if self.wavelength < color["max"]:
				# Good match.
				return color["name"]
		return "???"
	
	@property
	def CIE_XYZ(self):
		"""
		The CIE XYZ tristimulus values of the 'color' that corresponds to the energy of this excited state (as a numpy array).
		"""
		try:
			return colour.wavelength_to_XYZ(self.wavelength)
		except ValueError:
			# Wavelength is out of our colour range, so we can't see it.
			return numpy.zeros(3)
		
	@property
	def CIE_xy(self):
		"""
		The CIE xy chromaticity coordinates of the 'color' that corresponds to the energy of this excited state.
		"""
		try:
			return colour.XYZ_to_xy(self.CIE_XYZ)
		except Exception:
			return numpy.zeros(2)
		
				
	@property
	def rgb(self):
		"""
		The RGB values of the 'color' that corresponds to the energy of this excited state (as a list of [r, g, b] from 0 -> 255).
		"""
		return self.xyz_to_rgb(self.CIE_XYZ)
		
	@classmethod
	def xyz_to_rgb(self, XYZ):
		rgb =  [numpy.clip(clr, 0, math.inf) for clr in colour.XYZ_to_sRGB(XYZ)]
		
		# Now we normalise if one of our values exceeds 1.
		if max(rgb) > 1:
			rgb = [clr / max(rgb) for clr in rgb]
			
		# Now convert to 0 -> 255 and return.
		return [int(clr * 255) for clr in rgb]
		
	
	def set_file_options(self, output_dir, output_name, **kwargs):
		"""
		Set the options that will be used to create images from this object.
		
		:param output_dir: A pathlib Path object to the directory within which our files should be created.
		:param output_name: A string that will be used as the start of the file name of the files we create.
		"""
		# Set our TDP if we have one.
		if self.transition_dipole_moment is not None:
			self.transition_dipole_moment.set_file_options(output_dir, output_name, **kwargs)
	
	
	@classmethod
	def wavelength_to_energy(self, emission_wavelength):
		"""
		Convert an emission wavelength (in nm) to energy (in eV).
		"""
		#TODO: No reason for this method to have 'emission' in its name (it's very old).
		# e = (c * h) / λ
		return ((self.speed_of_light * self.plancks_constant) / (emission_wavelength / 1000000000)) / self.electron_volt
	
	@classmethod
	def energy_to_wavelength(self, energy):
		"""
		Convert an energy (in eV) to wavelength (in nm).
		"""
		# λ = (c * h) / e
		return ((self.speed_of_light * self.plancks_constant) / (energy * self.electron_volt)) * 1000000000
	
	@classmethod
	def wavenumbers_to_energy(self, wavenumbers):
		"""
		Convert wavenumbers (in cm-1) to energy (in eV).
		"""
		return self.wavelength_to_energy((1 / wavenumbers) * 10000000)
	
		
	@classmethod
	def list_from_parser(self, parser):
		"""
		Create a list of Excited_state objects from an output file parser.
				
		:param parser: An output file parser.
		:return: A list of Excited_state objects in the same order as given in cclib.
		"""
		try:
			# List of excited states.
			excited_states = []
			# Dictionary of each of our multiplicities (singlets, triplets etc) so we know our multiplicities.
			multiplicities = {}
			
			# Assemble cclib's various arrays into a single list.
			excited_states_data = list(zip(parser.data.etsyms, parser.data.etenergies, parser.data.etoscs))	
			
			# Loop through our data.
			for index, (symmetry, energy, oscillator_strength) in enumerate(excited_states_data):
				# Get our multiplicity.
				mult = Excited_state.get_multiplicity_from_symmetry(symmetry)
				# Try and see if we've already got some states of this mult.
				mult_level = 1
				try:
					# Add one to the level.
					multiplicities[mult] += 1
					mult_level = multiplicities[mult]
				except KeyError:
					# We haven't seen this mult before.
					multiplicities[mult] = 1
					
				# Relevant transition dipole moments.
				try:
					tdm = parser.results.transition_dipole_moments[index]
				except IndexError:
					tdm = None
					
				# Get and append our object.
				excited_states.append(
					self(
						level = index +1,
						multiplicity_level = mult_level,
						symmetry = symmetry,
						energy = self.wavenumbers_to_energy(energy),
						oscillator_strength = oscillator_strength,
						transitions = Excited_state_transition.list_from_parser(parser),
						transition_dipole_moment = tdm
					)
				)
			
			# All done, return our list.
			return excited_states
		except AttributeError:
			return []

	