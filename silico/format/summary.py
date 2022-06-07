# Summary formats.
# These format classes form the basis for several other formats, they extract a summary view of a particular result (or subsection of a result).
# Summary formats alone are sufficient in most use cases, but they are not great for repetitive result data (atoms, vibrations and excited states). 

# General imports.
from _collections import OrderedDict

# Silico imports.
from silico.exception.base import Result_unavailable_error
from silico.format import Result_format
import silico.format
from silico.result.angle import Angle
from silico.format.base import Result_format_group
from silico.misc import Layered_dict
from silico import misc


class Summary_group_format(Result_format_group):
    """
    Abstract top-level class for extracting particular sections from a number of Result_set objects in summary format.    
    """
    
    def __init__(self, *format_objects, ignore = False, **kwargs):
        """
        Constructor for Summary_group_format objects.
        """        
        # Check to see if there any DEF formats, and replace them if there are.
        replaced_formats = []
        for format_object in format_objects:
            if isinstance(format_object, Default_summary_format):
                # Replaced with defaults.
                replaced_formats.extend(self.get_default_formats())
                
                # And ignore
                ignore = True
            else:
                # Not default, just add.
                replaced_formats.append(format_object)
        format_objects = replaced_formats
        
        # Check to make sure we always have either name or metadata.
        if len(format_objects) != 0 and not True in [silico.format.METADATA_CLASS_HANDLE == format.CLASS_HANDLE for format in format_objects]:
            format_objects = list(format_objects)
            format_objects.insert(0, Name_summary_format(**kwargs))
        
        super().__init__(*format_objects, ignore = ignore, **kwargs)
        self.fieldnames = []
        
    @classmethod
    def get_default_formats(self, **kwargs):
        """
        Get a list of default format objects that can be used to convert a result set to summary format.
        """
        return [
            Metadata_summary_format(**kwargs),
            Vibrations_summary_format(**kwargs),
            Geometry_summary_format(**kwargs),
            SCF_summary_format(**kwargs),
            MP_summary_format(**kwargs),
            CC_summary_format(**kwargs),
            Orbitals_summary_format("+0", **kwargs),
            Orbitals_summary_format("+1", **kwargs),
            Orbitals_summary_format(**kwargs),
            Beta_summary_format("+0", **kwargs),
            Beta_summary_format("+1", **kwargs),
            Beta_summary_format(**kwargs),
            PDM_summary_format(**kwargs),
            TDM_summary_format(**kwargs),
            Excited_states_summary_format(**kwargs),
            Excited_states_summary_format((1,1), **kwargs),
            Excited_state_transitions_summary_format((1,1), 1, **kwargs),
            Excited_states_summary_format((3,1), **kwargs),
            Excited_state_transitions_summary_format((3,1), 1, **kwargs),
            SOC_summary_format("S(0)", "T(1)", **kwargs),
            SOC_summary_format("S(1)", "T(1)", **kwargs)
        ]
    
    def join(self, extracted_results):
        """
        Join together results from multiple format classes.
        
        :param extracted_results: A list of results extracted by this format group (one for each format). The format of each extracted result depends on the format group class.
        :return: A joined representation of the results. The format will depend on the format group class.
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

class Summary_format(Result_format):
    """
    Abstract, top-level class for all summary type format classes.
    
    Summary formats return Layered_dict objects.
    
    See the formats in text, CSV or table for concrete implementations.
    """
    
    def _extract(self, result, **kwargs):
        """
        Extract data from a Result_set in summary format.
        
        This default implementation does nothing, inheriting classes should write their own.
        The returned data is expected to be an OrderedDict or Layered_dict.
        
        Each key in the dict should be a 2-membered tuple of the form (section_name, property_name).
        """
        raise NotImplementedError()
    
class Default_summary_format(Summary_format):
    """
    Default format.
    
    This dummy format is replaced with a list of default formats.
    """
    
    CLASS_HANDLE = ["DEF", "default", "NORM", "normal"]

class Name_summary_format(Summary_format):
    """
    Summary format for result name.
    
    This dummy format is used to ensure the result name is always included no matter what the user chooses.
    """
    
    CLASS_HANDLE = ["NAME"]
    
    def _extract(self, result):
        """
        Convert a Result_set into an OrderedDict object.
        """
        return OrderedDict({
            'Name': result.safe_get('metadata', 'name')
        })

class Metadata_summary_format(Summary_format):
    """
    Summary format for metadata.
    """
    
    CLASS_HANDLE = silico.format.METADATA_CLASS_HANDLE
    
    def _extract(self, result):
        """
        Convert a Result_set into an OrderedDict object.
        """
        return OrderedDict({
            'Name': result.safe_get('metadata', 'name'),
            'Date': misc.date_to_string(result.safe_get('metadata', 'date')) if result.safe_get('metadata', 'date') is not None else None,
            'Duration /s': result.safe_get('metadata', 'duration').total_seconds() if result.safe_get('metadata', 'duration') is not None else None,
            'Package': result.safe_get('metadata', 'package'),
            'Package version': result.safe_get('metadata', 'package_version'),
            'Calculations': result.safe_get('metadata', 'calculations_string'),
            'Method': result.safe_get('metadata', 'methods_string'),
            'Functional': result.safe_get('metadata', 'functional'),
            'Basis set': result.safe_get('metadata', 'basis_set'),
            'Success': result.safe_get('metadata', 'success'),
            'Optimisation converged': result.safe_get('metadata', 'optimisation_converged'),
            'Calculation temperature /K': result.safe_get('metadata', 'temperature'),
            'Calculation pressure /atm': result.safe_get('metadata', 'pressure'),
            'Charge': result.safe_get('alignment', 'charge'),
            'Multiplicity': result.safe_get('metadata', 'multiplicity')
        })
        
class Geometry_summary_format(Summary_format):
    """
    Summary format for geometry (atoms etc).
    """
    
    CLASS_HANDLE = silico.format.GEOM_CLASS_HANDLE
    
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

class _Orbitals_summary_format(Summary_format):
    """
    Abstract class for fetching orbitals in summary format. See Orbitals_summary_format and Beta_summary_format for concrete classes.
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
            'HOMO-LUMO{} energy /eV'.format(spin_type): getattr(result, self.ORBITAL_TYPE).HOMO_LUMO_energy,
            'No. virtual orbitals{}'.format(spin_type): len([orbital for orbital in getattr(result, self.ORBITAL_TYPE) if orbital.HOMO_difference > 0]),
            'No. occupied orbitals{}'.format(spin_type): len([orbital for orbital in getattr(result, self.ORBITAL_TYPE) if orbital.HOMO_difference <= 0])
        })
        
class Orbitals_summary_format(_Orbitals_summary_format):
    """
    Summary format for (alpha) orbitals.
    """
    CLASS_HANDLE = silico.format.ORBITALS_CLASS_HANDLE
    
class Beta_summary_format(_Orbitals_summary_format):
    """
    Summary format for (alpha) orbitals.
    """
    ORBITAL_TYPE = "beta_orbitals"
    CLASS_HANDLE = silico.format.BETA_CLASS_HANDLE
    
class Vibrations_summary_format(Summary_format):
    """
    Summary format for vibrations.
    """
    ALLOW_CRITERIA = True
    CLASS_HANDLE = silico.format.VIBRATIONS_CLASS_HANDLE
    
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
        
class Excited_states_summary_format(Summary_format):
    """
    dict format for ES.
    """
    ALLOW_CRITERIA = True
    CLASS_HANDLE = silico.format.EXCITED_STATE_CLASS_HANDLE
    
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
            'ΔE(ST) /eV': result.safe_get('excited_states', 'singlet_triplet_energy'),
            'No. excited states': len(result.excited_states)
        })
        
class Emission_summary_format(Summary_format):
    
    ALLOW_CRITERIA = True
    
    def _extract_with_criteria(self, multiplicity, *, result):
        """
        Convert a Result_set into an OrderedDict object.
        
        """
        multiplicity = float(multiplicity)
        
        try:
            emission = getattr(result, self.EMISSION_ATTR)[multiplicity]
        except KeyError:
            raise Result_unavailable_error("adiabatic emission", "there are no emission with multiplicity '{}'".format(multiplicity))
            
        title = emission.transition_type.capitalize() + " " + emission.state_symbol + " " + "Emission"
        emission_rate = emission.safe_get('emission_rate')
        
        return OrderedDict({
            '{} energy /eV'.format(title): emission.energy,
            '{} wavelength /nm'.format(title): emission.wavelength,
            '{} colour'.format(title): emission.color,
            '{} oscillator strength'.format(title): emission.oscillator_strength,
            '{} rate /s-1'.format(title): emission_rate,
        })
        
class Adiabatic_emission_format(Emission_summary_format):
    
    CLASS_HANDLE = silico.format.ADIABATIC_EMISSION_CLASS_HANDLE
    EMISSION_ATTR = 'adiabatic_emission'
    
class Vertical_emission_format(Emission_summary_format):
    
    CLASS_HANDLE = silico.format.VERTICAL_EMISSION_CLASS_HANDLE
    EMISSION_ATTR = 'vertical_emission'
        
class SOC_summary_format(Summary_format):
    """
    dict format for ES.
    """
    ALLOW_CRITERIA = True
    CLASS_HANDLE = silico.format.SOC_CLASS_HANDLE
    
    def _extract_with_criteria(self, state1, state2, *, result):
        """
        Convert a Result_set into an OrderedDict object.
        
        """
        soc = result.spin_orbit_coupling.between(state1, state2)
        
        return OrderedDict({
            '<{}|Hso|{}> /cm-1'.format(soc.singlet_state.state_symbol, soc.triplet_state.state_symbol): soc.wavenumbers,
            '<{}|λ|{}> /cm-1'.format(soc.singlet_state.state_symbol, soc.triplet_state.state_symbol): soc.mixing_coefficient,
        })
        
    
class Excited_state_transitions_summary_format(Summary_format):
    """
    Summary format for ES transitions.
    """
    ALLOW_CRITERIA = True
    FORCE_CRITERIA = True
    CLASS_HANDLE = silico.format.EXCITED_STATE_TRANSITIONS_CLASS_HANDLE
    
    def _extract_transition(self, excited_state, transition = None):
        """
        Format helper function, extracts from a given state and transition.
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
            
    
class _Dipole_summary_format(Summary_format):
    """
    Abstract class for extracting dipole moments to summary format.
    """
    ALLOW_CRITERIA = True
    
    def _format_dipole(self, dipole_moment):
        """
        Convert a Result_set into an OrderedDict object.
        
        This abstract class method takes the dipole_moment to extract as its argument.
        """
        if dipole_moment is None:
            raise Result_unavailable_error("dipole moment", "there is no dipole of the requested type")
        
        angle_symbol = Angle.units_to_pretty_units(Angle._default_angle_units)
        return OrderedDict({
            '{} Vector X coord /D'.format(dipole_moment.name): dipole_moment.vector_coords[0],
            '{} Vector Y coord /D'.format(dipole_moment.name): dipole_moment.vector_coords[1],
            '{} Vector Z coord /D'.format(dipole_moment.name): dipole_moment.vector_coords[2],
            '{} x-axis angle /{}'.format(dipole_moment.name, angle_symbol): float(Angle(dipole_moment.X_axis_angle, output_units = Angle._default_angle_units)),
            '{} xy-plane angle /{}'.format(dipole_moment.name, angle_symbol): float(Angle(dipole_moment.XY_plane_angle, output_units = Angle._default_angle_units)), 
            '{} Total /D'.format(dipole_moment.name): dipole_moment.total
        })

class TDM_summary_format(_Dipole_summary_format):
    """
    Summary format for transition dipole moments.
    """
    
    CLASS_HANDLE = silico.format.TDM_CLASS_HANDLE
    
    def _format_dipole(self, dipole_moment):
        if dipole_moment is None:
            raise Result_unavailable_error("dipole moment", "there is no dipole of the requested type")
        
        angle_symbol = Angle.units_to_pretty_units(Angle._default_angle_units)
        
        return OrderedDict({
            # Electric
            "Electric {} vector x coord /D".format(dipole_moment.name): dipole_moment.electric.vector_coords[0],
            "Electric {} vector y coord /D".format(dipole_moment.name): dipole_moment.electric.vector_coords[1],
            "Electric {} vector z coord /D".format(dipole_moment.name): dipole_moment.electric.vector_coords[2], 
            "Electric {} /D".format(dipole_moment.name): dipole_moment.electric.total,
            "Electric {} x-axis angle /{}".format(dipole_moment.name, angle_symbol): float(Angle(dipole_moment.electric.X_axis_angle, output_units = Angle._default_angle_units)),
            "Electric {} xy-plane angle /{}".format(dipole_moment.name, angle_symbol): float(Angle(dipole_moment.electric.XY_plane_angle, output_units = Angle._default_angle_units)),
            # Magnetic
            "Magnetic {} vector x coord /au".format(dipole_moment.name): dipole_moment.magnetic.vector_coords[0],
            "Magnetic {} vector y coord /au".format(dipole_moment.name): dipole_moment.magnetic.vector_coords[1],
            "Magnetic {} vector z coord /au".format(dipole_moment.name): dipole_moment.magnetic.vector_coords[2], 
            "Magnetic {} /au".format(dipole_moment.name): dipole_moment.magnetic.total,
            "Magnetic {} x-axis angle /{}".format(dipole_moment.name, angle_symbol): float(Angle(dipole_moment.magnetic.X_axis_angle, output_units = Angle._default_angle_units)),
            "Magnetic {} xy-plane angle /{}".format(dipole_moment.name, angle_symbol): float(Angle(dipole_moment.magnetic.XY_plane_angle, output_units = Angle._default_angle_units)),
            # Contrast
            "Electric {} /esu⋅cm".format(dipole_moment.name): dipole_moment.electric.gaussian_cgs,
            "Magnetic {} /erg⋅G-1".format(dipole_moment.name): dipole_moment.magnetic.gaussian_cgs,
            "Electric & magnetic {} angle /{}".format(dipole_moment.name, angle_symbol): dipole_moment.angle().angle,
            "Electric & magnetic {} cos(angle)".format(dipole_moment.name): dipole_moment.cos_angle(),
            "Electric & magnetic {} dissymmetry factor".format(dipole_moment.name): dipole_moment.g_value,
        })
    
    def _extract_with_criteria(self, excited_state_id, *, result):
        """
        Convert a Result_set into an OrderedDict object.
        
        :param excited_state_id: A string or int describing the excited state to get the dipole of. See Excited_state.get_state().
        """
        # Get the dipole we've been asked for.
        # Criteria works the same as for excited states (because TDMs are stored under their respective ES).
        dipole_moment = result.excited_states.get_state(excited_state_id).transition_dipole_moment
        return self._format_dipole(dipole_moment)
    
    def _extract(self, result):
        """
        Convert a Result_set into an OrderedDict object.
        """
        return self._format_dipole(dipole_moment = result.transition_dipole_moment)
    
class PDM_summary_format(_Dipole_summary_format):
    """
    Summary format for permanent dipole moments.
    """
    ALLOW_CRITERIA = False
    CLASS_HANDLE = silico.format.PDM_CLASS_HANDLE
    
    def _extract(self, result):
        """
        Convert a Result_set into an OrderedDict object.
        """
        return super()._format_dipole(result.dipole_moment)
    
class _Energy_summary_format(Summary_format):
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
        
class SCF_summary_format(_Energy_summary_format):
    """
    Summary format for SCF energy.
    """
    CLASS_HANDLE = silico.format.SCF_CLASS_HANDLE
    ENERGY_TYPE = 'SCF_energies'
        
class MP_summary_format(_Energy_summary_format):
    """
    Summary format for MP energy.
    """
    CLASS_HANDLE = silico.format.MP_CLASS_HANDLE
    ENERGY_TYPE = 'MP_energies'
    
class CC_summary_format(_Energy_summary_format):
    """
    Summary format for CC energy.
    """
    CLASS_HANDLE = silico.format.CC_CLASS_HANDLE
    ENERGY_TYPE = 'CC_energies'
        