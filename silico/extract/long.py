# Extractors for tabulating data from a single result
# In contrast to the normal table formats (dict, CSV and table), long tables are not used for comparing different result files, but rather for tabulating data from one result file (such as orbitals, atoms, excited states etc).

from silico.extract import Result_extractor
from collections import OrderedDict
from silico.extract import Result_extractor_group
import silico.extract
from silico.misc import Layered_dict
from silico.exception.base import Result_unavailable_error, Silico_exception
from silico.result.angle import Angle
from silico.result.spectroscopy import Absorption_emission_graph,\
	Spectroscopy_graph


class Long_tabular_group_extractor(Result_extractor_group):
	"""
	Abstract group extractor for long tables (good for reading orbitals, atoms etc in text files and terminals).
	"""
	
	@classmethod
	def get_default_extractors(self, **kwargs):
		"""
		Get a list of default extractor objects that can be used to convert a result set to dict format.
		
		:param **kwargs: Keyword args that will be passed as is to each extractor class to construct.
		"""
		return [
			Atoms_long_extractor(),
			SCF_long_extractor(),
			MP_long_extractor(),
			CC_long_extractor(),
			Orbitals_long_extractor(),
			Beta_long_extractor(),
			Vibrations_long_extractor(),
			Excited_state_long_extractor(),
			Excited_state_transitions_long_extractor(),
			TDM_long_extractor()
		]
	
	@classmethod
	def tabulate(self, fieldnames, table_data):
		"""
		Convert some raw table data into our desired format.
		
		Inheriting classes should write their own implementation.
		
		:param fieldnames: List of fieldnames/table headers.
		:param table_data: List of table rows.
		"""
		raise NotImplementedError()
	
	def join(self, extracted_results):
		"""
		Join together results from multiple extractor classes.
		
		:param extracted_results: A list of results extracted by this extractor group (one for each extractor). The format of each extracted result depends on the extractor group class.
		:return: A joined representation of the results. The format will depend on the extractor group class.
		"""
		# We can make use of Layered_dict again to tabulate for us.
		text_tables = []
		# Tabulating is much easier here than for summary tables (because each row has the same keys).
		for extracted_result in extracted_results:
			if extracted_result != None:
				fieldnames, table_data = Layered_dict.tabulate(extracted_result)
				
				# Now turn into text with our class' tabulate() function of choice.
				text_tables.append(self.tabulate(fieldnames, table_data))
		
		# Join multiple sections by empty newline.
		return "\n".join(text_tables)
	
	def join_results(self, extracted_results):
		"""
		Method called to combine a list of extracted results from multiple result sets.
		"""
		# Just use our parents join() (which joins with newlines).
		return super().join(extracted_results)
	

class Long_table_extractor(Result_extractor):
	"""
	Top-level class for long table extractors.
	"""
	# Almost by definition, long extractors don't allow criteria (because they are for tabulating all data of a type).
	ALLOW_CRITERIA = False
	LIST_NAME = None
	
	def __init__(self, *args, **kwargs):
		"""
		Constructor for long tables.
		
		The constructor signature can differ for each extractor; positional args are used to specify criteria (filters which instruct the extractor what data to retrieve. Some extractors take no criteria, some optionally take criteria and some require criteria. The meaning of keyword args depends entirely on the implementing class.
		:param *args: Possibly optional criteria. Each criterion should be a (possibly length 1) tuple, the format depends on the inheriting class.
		:param **kwargs: Possibly optional keyword arguments. See the inheriting class for meaning.
		"""
		super().__init__(*args, **kwargs)
		#self.fieldnames = []

	def _extract_item(self, index, item):
		"""
		Helper method called to convert one item in the list we convert to an OrderedDict.
		
		Inheriting classes should write their own implementation.
		"""
		raise NotImplementedError()
	
	def _get_list(self, result):
		"""
		Method called to get the list over which this object will extract.
		
		This default implementation returns the attribute of result named by self.LIST_NAME
		"""
		return getattr(result, self.LIST_NAME)

	def _extract(self, result):
		"""
		Extract a result set.
		
		This implementation calls _extract_item() once for every item in the list returned by _get_list().
		"""
		# First, get our list.
		lst = self._get_list(result)
		
		# Get upset if the list is empty.
		if len(lst) == 0:
			raise Result_unavailable_error(self.LIST_NAME, "there are no items")
		
		# Extract using _extract_item().
		return [self._extract_item(index, item) for index, item in enumerate(lst)]
		
class Atoms_long_extractor(Long_table_extractor):
	"""
	Long table extractor for atoms.
	"""
	CLASS_HANDLE = silico.extract.GEOM_CLASS_HANDLE
	LIST_NAME = "alignment"
	
	def _extract_item(self, index, atom):
		"""
		"""
		return OrderedDict({
			"Atom No.": index +1,
			"Element": str(atom.element),
			"Proton No.": atom.element.number,
			"Isotope mass /gmol-1": atom.mass,
			"X coord": atom.coords[0],
			"Y coord": atom.coords[1],
			"Z coord": atom.coords[2]
			
		})
		
class _Orbitals_long_extractor(Long_table_extractor):
	"""
	Abstract class for fetching orbitals in long table format. See Orbitals_dict_extractor and Beta_dict_extractor for concrete classes. 
	"""
	LIST_NAME = "molecular_orbitals"
	
	def _extract_item(self, index, orbital):
		"""
		Convert an orbital object to an OrderedDict.
		"""
		return OrderedDict({
			'Label': orbital.label,
			'Level': orbital.level,
			'HOMO difference': orbital.HOMO_difference,
			'Symmetry': orbital.symmetry,
			'Energy /eV': orbital.energy
		})
		
	def _extract(self, result):
		"""
		Extract a result set.
		"""
		orbital_list = super()._extract(result)
		# We sort so we decrease in energy, might seem a bit weird but helps in comparing HOMO/LUMO levels.
		orbital_list.sort(key=lambda orb: orb['Level'], reverse=True)
		return orbital_list
		
class Orbitals_long_extractor(_Orbitals_long_extractor):
	"""
	Long table extractor for (alpha) orbitals.
	"""
	CLASS_HANDLE = silico.extract.ORBITALS_CLASS_HANDLE
	
class Beta_long_extractor(_Orbitals_long_extractor):
	"""
	Long table extractor for (alpha) orbitals.
	"""
	LIST_NAME = "beta_orbitals"
	CLASS_HANDLE = silico.extract.BETA_CLASS_HANDLE
	
class Vibrations_long_extractor(Long_table_extractor):
	"""
	Long table extractor for vibrations.
	"""
	CLASS_HANDLE = silico.extract.VIBRATIONS_CLASS_HANDLE
	LIST_NAME = "vibrations"
	
	def _extract_item(self, index, vibration):
		"""
		Convert a vibration object to an OrderedDict.
		"""
		return OrderedDict({
			'Level': vibration.level,
			'Symmetry': vibration.symmetry,
			'Frequency /cm-1': vibration.frequency,
			'Intensity /km mol-1': vibration.intensity
		})

class _Energy_long_extractor(Long_table_extractor):
	"""
	Abstract class for energy in long format.
	"""
	
	def _extract_item(self, index, energy):
		"""
		Convert an energy to an OrderedDict.
		"""
		return OrderedDict({
			"Step No.": index +1,
			'Energy /eV': energy
		})

		
class SCF_long_extractor(_Energy_long_extractor):
	"""
	Long table extractor for SCF energy.
	"""
	CLASS_HANDLE = silico.extract.SCF_CLASS_HANDLE
	LIST_NAME = 'SCF_energies'
		
class MP_long_extractor(_Energy_long_extractor):
	"""
	Long table extractor for MP energy.
	"""
	CLASS_HANDLE = silico.extract.MP_CLASS_HANDLE
	LIST_NAME = 'MP_energies'
	
class CC_long_extractor(_Energy_long_extractor):
	"""
	Long table extractor for CC energy.
	"""
	CLASS_HANDLE = silico.extract.CC_CLASS_HANDLE
	LIST_NAME = 'CC_energies'
	
class Excited_state_long_extractor(Long_table_extractor):
	"""
	Long table extractor for excited states.
	"""
	CLASS_HANDLE = silico.extract.EXCITED_STATE_CLASS_HANDLE
	LIST_NAME = "excited_states"
	
	def _extract_item(self, index, excited_state):
		"""
		Convert an excited_state to an OrderedDict.
		"""
		# First get general excited state stuff.
		data = Layered_dict({
			"Level": excited_state.level,
			"Symbol": excited_state.state_symbol,
			"Symmetry": excited_state.symmetry,
			"Energy /eV": excited_state.energy,
			"Wavelength /nm": excited_state.wavelength,
			"Colour": excited_state.color,
			"Oscillator strength": excited_state.oscillator_strength
		})
		
		# Now add data for each transition.
		# Layered_dict will take care of maintaining the order for us.
		for transition in excited_state.transitions:
			data.append(OrderedDict({
				"Transition {} orbitals:".format(transition.level): "{} -> {}".format(transition.starting_mo.label, transition.ending_mo.label),
				"Transition {} probability:".format(transition.level): transition.probability
			}))

# 		data.update({
# 			'Transitions (probabilty)': "\n".join(["{} -> {} ({:0.2f})".format(transition.starting_mo.label, transition.ending_mo.label, transition.probability) for transition in excited_state.transitions])
# 		})
	
		return data
	
class TDM_long_extractor(Long_table_extractor):
	"""
	Long table extractor for transition dipole moments.
	"""
	CLASS_HANDLE = silico.extract.TDM_CLASS_HANDLE
	# We don't actually use LIST_NAME for this class, but we still use its name for an error message.
	LIST_NAME = "TDM"
	
	def _get_list(self, result):
		"""
		Method called to get the list over which this object will extract.
		"""
		return [excited_state.transition_dipole_moment for excited_state in result.excited_states if excited_state.transition_dipole_moment is not None]
	
	def _extract_item(self, index, dipole_moment):
		"""
		Convert an excited_state to an OrderedDict.
		"""
		angle_symbol = Angle.units_to_pretty_units(Angle._default_angle_units)
		return OrderedDict({
			"Excited state level": dipole_moment.excited_state.level,
			"Excited state symbol": dipole_moment.excited_state.state_symbol,
			"Origin X coord /D": dipole_moment.origin_coords[0],
			"Origin Y coord /D": dipole_moment.origin_coords[1],
			"Origin Z coord /D": dipole_moment.origin_coords[2],
			"Vector X coord /D": dipole_moment.vector_coords[0],
			"Vector Y coord /D": dipole_moment.vector_coords[1],
			"Vector Z coord /D": dipole_moment.vector_coords[2],
			"X-axis angle /{}".format(angle_symbol): float(Angle(dipole_moment.X_axis_angle, output_units = Angle._default_angle_units)),
			"XY-plane angle /{}".format(angle_symbol): float(Angle(dipole_moment.XY_plane_angle, output_units = Angle._default_angle_units)), 
			"Total /D": dipole_moment.total
		})
		
class Excited_state_transitions_long_extractor(Long_table_extractor):
	"""
	Long table extractor for ES transitions.
	"""
	CLASS_HANDLE = silico.extract.EXCITED_STATE_TRANSITIONS_CLASS_HANDLE
	LIST_NAME = "transitions"
	
	def _get_list(self, result):
		"""
		Method called to get the list over which this object will extract.
		"""
		return [(excited_state, transition) for excited_state in result.excited_states for transition in excited_state.transitions]
	
	def _extract_item(self, index, transition_tuple):
		"""
		Convert an excited_state to an OrderedDict.
		"""
		excited_state = transition_tuple[0]
		transition = transition_tuple[1]
		
		
		return OrderedDict({
			"Excited state level": excited_state.level,
			"Excited state symbol": excited_state.state_symbol,
			"Starting orbital": transition.starting_mo.label,
			"Ending orbital": transition.ending_mo.label,
			"Probability": transition.probability
		})
		
class Absorption_spectrum_long_extractor(Long_table_extractor):
	"""
	Extractor for simulated absorption graphs
	"""
	CLASS_HANDLE = ["ABS", "absorption", "absorption_graph"]
	LIST_NAME = "absorptions"
	
	def _get_list(self, result):
		"""
		Method called to get the list over which this object will extract.
		"""
		try:
			# First, lets get our spectrum!
			spectrum = Absorption_emission_graph.from_excited_states(result.excited_states, use_jacobian = self.config['absorption_spectrum']['use_jacobian'])
			
			return [(energy, intensity) for energy, intensity in spectrum.plot_cumulative_gaussian(self.config['absorption_spectrum']['fwhm'], self.config['absorption_spectrum']['gaussian_resolution'], self.config['absorption_spectrum']['gaussian_cutoff'])]
		except Silico_exception:
			# No values to plot.
			return []
	
	def _extract_item(self, index, coord):
		"""
		Convert a coordinate to an OrderedDict.
		"""
		return OrderedDict({'Wavelength /nm': coord[0], 'Intensity': coord[1]})
	
class Absorption_energy_spectrum_long_extractor(Long_table_extractor):
	"""
	Extractor for simulated absorption graphs
	"""
	CLASS_HANDLE = ["ABSE", "absorption_energy", "absorption_energy_graph"]
	LIST_NAME = "absorptions"
	
	def _get_list(self, result):
		"""
		Method called to get the list over which this object will extract.
		"""
		try:
			# First, lets get our spectrum!
			spectrum = Spectroscopy_graph([(excited_state.energy, excited_state.oscillator_strength) for excited_state in result.excited_states])
			
			return [(energy, intensity) for energy, intensity in spectrum.plot_cumulative_gaussian(self.config['absorption_spectrum']['fwhm'], self.config['absorption_spectrum']['gaussian_resolution'], self.config['absorption_spectrum']['gaussian_cutoff'])]
		except Silico_exception:
			# No values to plot.
			return []
	
	def _extract_item(self, index, coord):
		"""
		Convert a coordinate to an OrderedDict.
		"""
		return OrderedDict({'Energy /eV': coord[0], 'Oscillator Strength': coord[1]})

		
class IR_spectrum_long_extractor(Long_table_extractor):
	"""
	Extractor for simulated IR spectra
	"""
	CLASS_HANDLE = ["IR", "infrared"]
	LIST_NAME = "vibrations"
	
	def _get_list(self, result):
		"""
		Method called to get the list over which this object will extract.
		"""
		try:
			# First, lets get our spectrum!
			spectrum = Spectroscopy_graph.from_vibrations(result.vibrations)
			
			return [(energy, intensity) for energy, intensity in spectrum.plot_cumulative_gaussian(self.config['IR_spectrum']['fwhm'], self.config['IR_spectrum']['gaussian_resolution'], self.config['IR_spectrum']['gaussian_cutoff'])]
		except Silico_exception:
			# No values to plot.
			return []
	
	def _extract_item(self, index, coord):
		"""
		Convert a coordinate to an OrderedDict.
		"""
		return OrderedDict({'Frequency /cm-1': coord[0], 'Intensity /km mol-1': coord[1]})
	