# Summary extractors.
# These extractor classes form the basis for several other extractors, they extract a summary view of a particular result (or subsection of a result).
# Summary extractors alone are sufficient in most use cases, but they are not great for repetitive result data (atoms, vibrations and excited states). 

from _collections import OrderedDict
from silico.exception.base import Result_unavailable_error
from silico.extract import Result_extractor
import silico.extract
from silico.result.angle import Angle
from silico.extract.base import Result_extractor_group
from silico.misc import Layered_dict
from silico import misc

class Summary_group_extractor(Result_extractor_group):
	"""
	Abstract top-level class for extracting particular sections from a number of Result_set objects in summary format.	
	"""
	
	def __init__(self, *extractor_objects, ignore = False, **kwargs):
		"""
		Constructor for Summary_group_extractor objects.
		"""		
		# Check to see if there any DEF extractors, and replace them if there are.
		replaced_extractors = []
		for extractor_object in extractor_objects:
			if isinstance(extractor_object, Default_summary_extractor):
				# Replaced with defaults.
				replaced_extractors.extend(self.get_default_extractors())
				
				# And ignore
				ignore = True
			else:
				# Not default, just add.
				replaced_extractors.append(extractor_object)
		extractor_objects = replaced_extractors
		
		# Check to make sure we always have either name or metadata.
		if len(extractor_objects) != 0 and not True in [silico.extract.METADATA_CLASS_HANDLE == extractor.CLASS_HANDLE for extractor in extractor_objects]:
			extractor_objects = list(extractor_objects)
			extractor_objects.insert(0, Name_summary_extractor(**kwargs))
		
		super().__init__(*extractor_objects, ignore = ignore, **kwargs)
		self.fieldnames = []
		
	@classmethod
	def get_default_extractors(self, **kwargs):
		"""
		Get a list of default extractor objects that can be used to convert a result set to summary format.
		"""
		return [
			Metadata_summary_extractor(**kwargs),
			Vibrations_summary_extractor(**kwargs),
			Geometry_summary_extractor(**kwargs),
			SCF_summary_extractor(**kwargs),
			MP_summary_extractor(**kwargs),
			CC_summary_extractor(**kwargs),
			Orbitals_summary_extractor("+0", **kwargs),
			Orbitals_summary_extractor("+1", **kwargs),
			Orbitals_summary_extractor(**kwargs),
			Beta_summary_extractor("+0", **kwargs),
			Beta_summary_extractor("+1", **kwargs),
			Beta_summary_extractor(**kwargs),
			PDM_summary_extractor(**kwargs),
			TDM_summary_extractor(**kwargs),
			Excited_states_summary_extractor(**kwargs),
			Excited_states_summary_extractor((1,1), **kwargs),
			Excited_state_transitions_summary_extractor((1,1), 1, **kwargs),
			Excited_states_summary_extractor((3,1), **kwargs),
			Excited_state_transitions_summary_extractor((3,1), 1, **kwargs),
			SOC_summary_extractor("S(0)", "T(1)", **kwargs),
			SOC_summary_extractor("S(1)", "T(1)", **kwargs)
		]
	
	def join(self, extracted_results):
		"""
		Join together results from multiple extractor classes.
		
		:param extracted_results: A list of results extracted by this extractor group (one for each extractor). The format of each extracted result depends on the extractor group class.
		:return: A joined representation of the results. The format will depend on the extractor group class.
		"""
		return Layered_dict(*extracted_results)
		
	
	def join_results(self, extracted_results):
		"""
		Method called to combine a list of extracted results from multiple result sets.
		"""
		# Get our table data from Layered_dict.
		self.fieldnames, result_rows = Layered_dict.tabulate(extracted_results)
		
		# Return the table data.
		return result_rows

class Summary_extractor(Result_extractor):
	"""
	Abstract, top-level class for all summary type extractor classes.
	
	Summary extractors return Layered_dict objects.
	
	See the extractors in text, CSV or table for concrete implementations.
	"""
	
	def _extract(self, result, **kwargs):
		"""
		Extract data from a Result_set in summary format.
		
		This default implementation does nothing, inheriting classes should write their own.
		The returned data is expected to be an OrderedDict or Layered_dict.
		
		Each key in the dict should be a 2-membered tuple of the form (section_name, property_name).
		"""
		raise NotImplementedError()
	
class Default_summary_extractor(Summary_extractor):
	"""
	Default extractor.
	
	This dummy extractor is replaced with a list of default extractors.
	"""
	
	CLASS_HANDLE = ["DEF", "default", "NORM", "normal"]

class Name_summary_extractor(Summary_extractor):
	"""
	Summary extractor for result name.
	
	This dummy extractor is used to ensure the result name is always included no matter what the user chooses.
	"""
	
	CLASS_HANDLE = ["NAME"]
	
	def _extract(self, result):
		"""
		Convert a Result_set into an OrderedDict object.
		"""
		return OrderedDict({
			'Name': result.safe_get('metadata', 'name')
		})

class Metadata_summary_extractor(Summary_extractor):
	"""
	Summary extractor for metadata.
	"""
	
	CLASS_HANDLE = silico.extract.METADATA_CLASS_HANDLE
	
	def _extract(self, result):
		"""
		Convert a Result_set into an OrderedDict object.
		"""
		return OrderedDict({
			'Name': result.safe_get('metadata', 'name'),
			'Date': misc.date_to_string(result.safe_get('metadata', 'date')) if result.safe_get('metadata', 'date') is not None else None,
			'Duration': misc.timedelta_to_string(result.safe_get('metadata', 'duration')) if result.safe_get('metadata', 'duration') is not None else None,
			'Package': result.safe_get('metadata', 'package'),
			'Package version': result.safe_get('metadata', 'package_version'),
			'Calculations': result.safe_get('metadata', 'calculations_string'),
			'Method': result.safe_get('metadata', 'calc_methods_string'),
			'Functional': result.safe_get('metadata', 'calc_functional'),
			'Basis set': result.safe_get('metadata', 'calc_basis_set'),
			'Success': result.safe_get('metadata', 'calc_success'),
			'Optimisation converged': result.safe_get('metadata', 'optimisation_converged'),
			'Calculation temperature /K': result.safe_get('metadata', 'calc_temperature'),
			'Calculation pressure /atm': result.safe_get('metadata', 'calc_pressure'),
			#'Formula': result.safe_get('alignment', 'formula_string'),
			'Charge': result.safe_get('alignment', 'charge'),
			'Multiplicity': result.safe_get('metadata', 'system_multiplicity')
		})
		
class Geometry_summary_extractor(Summary_extractor):
	"""
	Summary extractor for geometry (atoms etc).
	"""
	
	CLASS_HANDLE = silico.extract.GEOM_CLASS_HANDLE
	
	def _extract(self, result):
		"""
		Convert a Result_set into an OrderedDict object.
		"""
		return OrderedDict({
			'Formula': result.alignment.formula_string,
			'Exact mass /gmol-1': result.safe_get('alignment', 'mass'),
			'Molar mass /gmol-1': result.alignment.molar_mass,
			'No. atoms': len(result.alignment),
			'Alignment method': result.alignment.CLASS_HANDLE[0],
			'X extension /Å': result.alignment.X_length,
			'Y extension /Å': result.alignment.Y_length,
			'Z extension /Å': result.alignment.Z_length,
			'Linearity ratio': result.alignment.get_linear_ratio(),
			'Planarity ratio': result.alignment.get_planar_ratio()
		})

class _Orbitals_summary_extractor(Summary_extractor):
	"""
	Abstract class for fetching orbitals in summary format. See Orbitals_summary_extractor and Beta_summary_extractor for concrete classes.
	"""
	ORBITAL_TYPE = "molecular_orbitals"
	ALLOW_CRITERIA = True
	
	def _extract_with_criteria(self, orbital_id, *, result):
		"""
		Convert a Result_set into an OrderedDict object.
		
		:param orbital_id: A string or int describing the orbital to extract. See Molecular_orbital.get_orbital().
		"""
		# First, get the orbitals we were asked for.
		orbital = getattr(result, self.ORBITAL_TYPE).get_orbital(orbital_id)
		
		# Now convert to dict.
		return OrderedDict({
			'{} energy /eV'.format(orbital.label): orbital.energy
		})
		
	def _extract(self, result):
		"""
		Convert a Result_set into an OrderedDict object.
		"""
		spin_type = " ({})".format(getattr(result, self.ORBITAL_TYPE).spin_type) if getattr(result, self.ORBITAL_TYPE).spin_type != "none" else ""
		return OrderedDict({
			#'HOMO{} energy /eV'.format(spin_type): getattr(result, self.ORBITAL_TYPE).HOMO_energy,
			#'LUMO{} energy /eV'.format(spin_type): getattr(result, self.ORBITAL_TYPE).LUMO_energy,
			'HOMO/LUMO{} energy /eV'.format(spin_type): getattr(result, self.ORBITAL_TYPE).HOMO_LUMO_energy,
			'No. virtual orbitals{}'.format(spin_type): len([orbital for orbital in getattr(result, self.ORBITAL_TYPE) if orbital.HOMO_difference > 0]),
			'No. occupied orbitals{}'.format(spin_type): len([orbital for orbital in getattr(result, self.ORBITAL_TYPE) if orbital.HOMO_difference <= 0])
		})
		
class Orbitals_summary_extractor(_Orbitals_summary_extractor):
	"""
	Summary extractor for (alpha) orbitals.
	"""
	CLASS_HANDLE = silico.extract.ORBITALS_CLASS_HANDLE
	
class Beta_summary_extractor(_Orbitals_summary_extractor):
	"""
	Summary extractor for (alpha) orbitals.
	"""
	ORBITAL_TYPE = "beta_orbitals"
	CLASS_HANDLE = silico.extract.BETA_CLASS_HANDLE
	
class Vibrations_summary_extractor(Summary_extractor):
	"""
	Summary extractor for vibrations.
	"""
	ALLOW_CRITERIA = True
	CLASS_HANDLE = silico.extract.VIBRATIONS_CLASS_HANDLE
	
	def _extract_with_criteria(self, vibration_index, *, result):
		"""
		Convert a Result_set into an OrderedDict object.
		"""
		try:
			vibration = result.vibrations[int(vibration_index) -1]
		except IndexError:
			raise Result_unavailable_error("vibrational frequency", "no vibration with index {}".format(vibration_index))
		
		return OrderedDict({
			'Vibration {} symmetry'.format(vibration.level): vibration.symmetry,
			'Vibration {} frequency /cm-1'.format(vibration.level): vibration.frequency,
			'Vibration {} intensity /km mol-1'.format(vibration.level): vibration.intensity
		})
		
	def _extract(self, result):
		"""
		Convert a Result_set into an OrderedDict object.
		"""
		if len(result.vibrations) == 0:
			raise Result_unavailable_error("vibrational frequencies", "there are no vibrational frequencies")
		
		return OrderedDict({
			'No. vibrations': len(result.vibrations),
			'No. negative frequencies': len(result.vibrations.negative_frequencies)
		})
		
class Excited_states_summary_extractor(Summary_extractor):
	"""
	dict extractor for ES.
	"""
	ALLOW_CRITERIA = True
	CLASS_HANDLE = silico.extract.EXCITED_STATE_CLASS_HANDLE
	
	def _extract_with_criteria(self, excited_state_id, *, result):
		"""
		Convert a Result_set into an OrderedDict object.
		
		:param excited_state_id: A string or int describing the excited state to get the dipole of. See Excited_state.get_state().
		"""
		es = result.excited_states.get_state(excited_state_id)
		
		return OrderedDict({
			'{} symmetry'.format(es.state_symbol): es.symmetry,
			'{} energy /eV'.format(es.state_symbol): es.energy,
			'{} wavelength /nm'.format(es.state_symbol): es.wavelength,
			'{} colour'.format(es.state_symbol): es.color,
			'{} CIE X'.format(es.state_symbol): es.CIE_xy[0],
			'{} CIE Y'.format(es.state_symbol): es.CIE_xy[1],
			'{} oscillator strength'.format(es.state_symbol): es.oscillator_strength
		})
		
	def _extract(self, result):
		"""
		Convert a Result_set into an OrderedDict object.
		"""
		if len(result.excited_states) == 0:
			raise Result_unavailable_error("excited states", "there are no excited states")
		
		return OrderedDict({
			'ΔEst /eV': result.safe_get('excited_states', 'singlet_triplet_energy'),
			'No. excited states': len(result.excited_states)
		})
		
		
class SOC_summary_extractor(Summary_extractor):
	"""
	dict extractor for ES.
	"""
	ALLOW_CRITERIA = True
	CLASS_HANDLE = silico.extract.SOC_CLASS_HANDLE
	
	def _extract_with_criteria(self, state1, state2, *, result):
		"""
		Convert a Result_set into an OrderedDict object.
		
		"""
		soc = result.spin_orbit_coupling.between(state1, state2)
		
		return OrderedDict({
			'<{}|Hso|{}> /cm-1'.format(soc.singlet_state.state_symbol, soc.triplet_state.state_symbol): soc.wavenumbers,
			'<{}|λ|{}> /cm-1'.format(soc.singlet_state.state_symbol, soc.triplet_state.state_symbol): soc.mixing_coefficient,
		})
		
# 	def _extract(self, result):
# 		"""
# 		Convert a Result_set into an OrderedDict object.
# 		"""
# 		if len(result.spin_orbit_coupling) == 0:
# 			raise Result_unavailable_error("SOC", "there is no spin-orbit coupling0")
# 		
# 		return OrderedDict({
# 			'ΔEst /eV': result.safe_get('excited_states', 'singlet_triplet_energy'),
# 			'No. excited states': len(result.excited_states)
# 		})
		
	
class Excited_state_transitions_summary_extractor(Summary_extractor):
	"""
	Summary extractor for ES transitions.
	"""
	ALLOW_CRITERIA = True
	FORCE_CRITERIA = True
	CLASS_HANDLE = silico.extract.EXCITED_STATE_TRANSITIONS_CLASS_HANDLE
	
	def _extract_transition(self, excited_state, transition = None):
		"""
		Extractor helper function, extracts from a given state and transition.
		"""
		# If we weren't given a transition, call ourself again, once for each possible transition, and combine into a single Layered_dict.
		if transition is None:
			return Layered_dict().update(*[self._extract_transition(excited_state, transition) for transition in excited_state.transitions])
		
		# Normal execution, extract the transition.
		return Layered_dict({
			"{} transition {} orbitals".format(excited_state.state_symbol, transition.level): "{} -> {}".format(transition.starting_mo.label, transition.ending_mo.label),
			"{} transition {} probability".format(excited_state.state_symbol, transition.level): transition.probability
		})
	
	def _extract_with_criteria(self, excited_state_id, transition_index = None, *, result):
		"""
		Convert a Result_set into an OrderedDict object.
		"""	
		# First get out excited state.
		excited_state = result.excited_states.get_state(excited_state_id)
		
		# Get our transition index, watching out for user error.
		try:
			transition_index = int(transition_index)-1 if transition_index is not None else None
		except Exception:
			raise TypeError("the excited state transition identifier '{}' is not valid".format(transition_index))
		
		# Now get the transition object; None is fine here because it is supported by _extract_transition().
		try:
			transition = excited_state.transitions[transition_index] if transition_index is not None else None
		except IndexError:
			raise Result_unavailable_error("excited state transition", "there is no transition no. '{}' for excited state '{}'".format(transition_index+1, excited_state_id))
		
		# Now extract proper.
		return self._extract_transition(excited_state, transition)
			
	
class _Dipole_summary_extractor(Summary_extractor):
	"""
	Abstract class for extracting dipole moments to summary format.
	"""
	ALLOW_CRITERIA = True
	
	def _extract_with_criteria(self, dipole_moment):
		"""
		Convert a Result_set into an OrderedDict object.
		
		This abstract class method takes the dipole_moment to extract as its argument.
		"""
		if dipole_moment is None:
			raise Result_unavailable_error("dipole moment", "there is no dipole of the requested type")
		
		angle_symbol = Angle.units_to_pretty_units(Angle._default_angle_units)
		return OrderedDict({
			'{} Origin X coord /D'.format(dipole_moment.name): dipole_moment.origin_coords[0],
			'{} Origin Y coord /D'.format(dipole_moment.name): dipole_moment.origin_coords[1],
			'{} Origin Z coord /D'.format(dipole_moment.name): dipole_moment.origin_coords[2],
			'{} Vector X coord /D'.format(dipole_moment.name): dipole_moment.vector_coords[0],
			'{} Vector Y coord /D'.format(dipole_moment.name): dipole_moment.vector_coords[1],
			'{} Vector Z coord /D'.format(dipole_moment.name): dipole_moment.vector_coords[2],
			'{} x-axis angle /{}'.format(dipole_moment.name, angle_symbol): float(Angle(dipole_moment.X_axis_angle, output_units = Angle._default_angle_units)),
			'{} xy-plane angle /{}'.format(dipole_moment.name, angle_symbol): float(Angle(dipole_moment.XY_plane_angle, output_units = Angle._default_angle_units)), 
			'{} Total /D'.format(dipole_moment.name): dipole_moment.total
		})

class TDM_summary_extractor(_Dipole_summary_extractor):
	"""
	Summary extractor for transition dipole moments.
	"""
	
	CLASS_HANDLE = silico.extract.TDM_CLASS_HANDLE
	
	def _extract_with_criteria(self, excited_state_id, *, result):
		"""
		Convert a Result_set into an OrderedDict object.
		
		:param excited_state_id: A string or int describing the excited state to get the dipole of. See Excited_state.get_state().
		"""
		# Get the dipole we've been asked for.
		# Criteria works the same as for excited states (because TDMs are stored under their respective ES).
		dipole = result.excited_states.get_state(excited_state_id).transition_dipole_moment
		return super()._extract_with_criteria(dipole_moment = dipole)
	
	def _extract(self, result):
		"""
		Convert a Result_set into an OrderedDict object.
		"""
		return super()._extract_with_criteria(dipole_moment = result.transition_dipole_moment)
	
class PDM_summary_extractor(_Dipole_summary_extractor):
	"""
	Summary extractor for permanent dipole moments.
	"""
	ALLOW_CRITERIA = False
	CLASS_HANDLE = silico.extract.PDM_CLASS_HANDLE
	
	def _extract(self, result):
		"""
		Convert a Result_set into an OrderedDict object.
		"""
		return super()._extract_with_criteria(result.dipole_moment)
	
class _Energy_summary_extractor(Summary_extractor):
	"""
	Abstract class for Energy_list data to summary format.
	"""
	ALLOW_CRITERIA = True
	ENERGY_TYPE = ''
	
	def _extract_with_criteria(self, energy_step, *, result):
		"""
		Convert a Result_set into an OrderedDict object.
		"""
		energies = getattr(result, self.ENERGY_TYPE)
		energy_step = int(energy_step)
		return OrderedDict({
			'{} step {} energy /eV'.format(energies.energy_type, energy_step +1): energies[energy_step]
		})
			
	def _extract(self, result):
		"""
		Convert a Result_set into an OrderedDict object.
		"""
		energies = getattr(result, self.ENERGY_TYPE)
		return OrderedDict({
			'No. of {} steps'.format(energies.energy_type): len(energies),
			'{} energy /eV'.format(energies.energy_type): energies.final,
			'{} energy /kJmol-1'.format(energies.energy_type): energies.eV_to_kJmol(energies.final)
		})
		
class SCF_summary_extractor(_Energy_summary_extractor):
	"""
	Summary extractor for SCF energy.
	"""
	CLASS_HANDLE = silico.extract.SCF_CLASS_HANDLE
	ENERGY_TYPE = 'SCF_energies'
		
class MP_summary_extractor(_Energy_summary_extractor):
	"""
	Summary extractor for MP energy.
	"""
	CLASS_HANDLE = silico.extract.MP_CLASS_HANDLE
	ENERGY_TYPE = 'MP_energies'
	
class CC_summary_extractor(_Energy_summary_extractor):
	"""
	Summary extractor for CC energy.
	"""
	CLASS_HANDLE = silico.extract.CC_CLASS_HANDLE
	ENERGY_TYPE = 'CC_energies'
		