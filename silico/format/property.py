# Formats for tabulating data from a single result
# In contrast to the normal table formats (dict, CSV and table), property tables are not used for comparing different result files, but rather for tabulating data from one result file (such as orbitals, atoms, excited states etc).

from silico.format import Result_format
from collections import OrderedDict
from silico.format import Result_format_group
import silico.format
from silico.misc import Layered_dict
from silico.exception.base import Result_unavailable_error, Silico_exception
from silico.result.angle import Angle
from silico.result.spectroscopy import Absorption_emission_graph,\
    Spectroscopy_graph


class Tabular_property_group_format(Result_format_group):
    """
    Abstract group format for property tables (good for reading orbitals, atoms etc in text files and terminals).
    """
    
    @classmethod
    def get_default_formats(self, **kwargs):
        """
        Get a list of default format objects that can be used to convert a result set to dict format.
        
        :param **kwargs: Keyword args that will be passed as is to each format class to construct.
        """
        return [
            Atoms_property_format(**kwargs),
            SCF_property_format(**kwargs),
            MP_property_format(**kwargs),
            CC_property_format(**kwargs),
            Orbitals_property_format(**kwargs),
            Beta_property_format(**kwargs),
            Vibrations_property_format(**kwargs),
            Excited_state_property_format(**kwargs),
            Excited_state_transitions_property_format(**kwargs),
            SOC_property_format(**kwargs),
            TDM_property_format(**kwargs)
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
        Join together results from multiple format classes.
        
        :param extracted_results: A list of results extracted by this format group (one for each format). The format of each extracted result depends on the format group class.
        :return: A joined representation of the results. The format will depend on the format group class.
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
    

class Property_table_format(Result_format):
    """
    Top-level class for property table formats.
    """
    # Almost by definition, property formats don't allow criteria (because they are for tabulating all data of a type).
    ALLOW_CRITERIA = False
    LIST_NAME = None
    
    def __init__(self, *args, **kwargs):
        """
        Constructor for property tables.
        
        The constructor signature can differ for each format; positional args are used to specify criteria (filters which instruct the format what data to retrieve. Some formats take no criteria, some optionally take criteria and some require criteria. The meaning of keyword args depends entirely on the implementing class.
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
        
class Atoms_property_format(Property_table_format):
    """
    Property table format for atoms.
    """
    CLASS_HANDLE = silico.format.GEOM_CLASS_HANDLE
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
        
class _Orbitals_property_format(Property_table_format):
    """
    Abstract class for fetching orbitals in property table format. See Orbitals_dict_format and Beta_dict_format for concrete classes. 
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
        # We sort so we decrease in energy, might seem a bit weird but helps in comparing HOMO-LUMO levels.
        orbital_list.sort(key=lambda orb: orb['Level'], reverse=True)
        return orbital_list
        
class Orbitals_property_format(_Orbitals_property_format):
    """
    Property table format for (alpha) orbitals.
    """
    CLASS_HANDLE = silico.format.ORBITALS_CLASS_HANDLE
    
class Beta_property_format(_Orbitals_property_format):
    """
    Property table format for (alpha) orbitals.
    """
    LIST_NAME = "beta_orbitals"
    CLASS_HANDLE = silico.format.BETA_CLASS_HANDLE
    
class Vibrations_property_format(Property_table_format):
    """
    Property table format for vibrations.
    """
    CLASS_HANDLE = silico.format.VIBRATIONS_CLASS_HANDLE
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

class _Energy_property_format(Property_table_format):
    """
    Abstract class for energy in property format.
    """
    
    def _extract_item(self, index, energy):
        """
        Convert an energy to an OrderedDict.
        """
        return OrderedDict({
            "Step No.": index +1,
            'Energy /eV': energy
        })

        
class SCF_property_format(_Energy_property_format):
    """
    Property table format for SCF energy.
    """
    CLASS_HANDLE = silico.format.SCF_CLASS_HANDLE
    LIST_NAME = 'SCF_energies'
        
class MP_property_format(_Energy_property_format):
    """
    Property table format for MP energy.
    """
    CLASS_HANDLE = silico.format.MP_CLASS_HANDLE
    LIST_NAME = 'MP_energies'
    
class CC_property_format(_Energy_property_format):
    """
    Property table format for CC energy.
    """
    CLASS_HANDLE = silico.format.CC_CLASS_HANDLE
    LIST_NAME = 'CC_energies'
    
class Excited_state_property_format(Property_table_format):
    """
    Property table format for excited states.
    """
    CLASS_HANDLE = silico.format.EXCITED_STATE_CLASS_HANDLE
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
            "CIE X": excited_state.CIE_xy[0],
            "CIE Y": excited_state.CIE_xy[1],
            "Oscillator strength": excited_state.oscillator_strength
        })
        
        # Now add data for each transition.
        # Layered_dict will take care of maintaining the order for us.
        for transition in excited_state.transitions:
            data.append(OrderedDict({
                "Transition {} orbitals:".format(transition.level): "{} -> {}".format(transition.starting_mo.label, transition.ending_mo.label),
                "Transition {} probability:".format(transition.level): transition.probability
            }))
        return data
    
class SOC_property_format(Property_table_format):
    """
    Property table format for SOC.
    """
    CLASS_HANDLE = silico.format.SOC_CLASS_HANDLE
    LIST_NAME = "spin_orbit_coupling"
    
    def _extract_item(self, index, soc):
        """
        """
        return OrderedDict({
            "Singlet": soc.singlet_state.state_symbol,
            "Triplet": soc.triplet_state.state_symbol,
            "+1 /cm-1": soc.positive_one,
            "0 /cm-1": soc.zero,
            "-1 /cm-1": soc.negative_one,
            "RSS /cm -1": soc.wavenumbers,
            "Hso /eV": soc.energy,
            "ΔE /eV": soc.splitting_energy,
            "λ": soc.mixing_coefficient,
            
        })
    
    
class TDM_property_format(Property_table_format):
    """
    Property table format for transition dipole moments.
    """
    CLASS_HANDLE = silico.format.TDM_CLASS_HANDLE
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
            # Meta
            "Excited state level": dipole_moment.electric.excited_state.level,
            "Excited state symbol": dipole_moment.electric.excited_state.state_symbol,
            # Electric
            "Electric dipole moment vector x coord /D": dipole_moment.electric.vector_coords[0],
            "Electric dipole moment vector y coord /D": dipole_moment.electric.vector_coords[1],
            "Electric dipole moment vector z coord /D": dipole_moment.electric.vector_coords[2], 
            "Electric dipole moment /D": dipole_moment.electric.total,
            "Electric dipole moment x-axis angle /{}".format(angle_symbol): float(Angle(dipole_moment.electric.X_axis_angle, output_units = Angle._default_angle_units)),
            "Electric dipole moment xy-plane angle /{}".format(angle_symbol): float(Angle(dipole_moment.electric.XY_plane_angle, output_units = Angle._default_angle_units)),
            # Magnetic
            "Magnetic dipole moment vector x coord /au": dipole_moment.magnetic.vector_coords[0],
            "Magnetic dipole moment vector y coord /au": dipole_moment.magnetic.vector_coords[1],
            "Magnetic dipole moment vector z coord /au": dipole_moment.magnetic.vector_coords[2], 
            "Magnetic dipole moment /au": dipole_moment.magnetic.total,
            "Magnetic dipole moment x-axis angle /{}".format(angle_symbol): float(Angle(dipole_moment.magnetic.X_axis_angle, output_units = Angle._default_angle_units)),
            "Magnetic dipole moment xy-plane angle /{}".format(angle_symbol): float(Angle(dipole_moment.magnetic.XY_plane_angle, output_units = Angle._default_angle_units)),
            # Contrast
            "Electric dipole moment /esu⋅cm": dipole_moment.electric.gaussian_cgs,
            "Magnetic dipole moment /erg⋅G-1": dipole_moment.magnetic.gaussian_cgs,
            "Electric & magnetic dipole moment angle /{}".format(angle_symbol): dipole_moment.angle().angle,
            "Electric & magnetic dipole moment cos(angle)": dipole_moment.cos_angle(),
            "Electric & magnetic dipole moment dissymmetry factor": dipole_moment.g_value,
        })
        
class Excited_state_transitions_property_format(Property_table_format):
    """
    Property table format for ES transitions.
    """
    CLASS_HANDLE = silico.format.EXCITED_STATE_TRANSITIONS_CLASS_HANDLE
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
        
class Absorption_spectrum_property_format(Property_table_format):
    """
    Format for simulated absorption graphs
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
    
class Absorption_energy_spectrum_property_format(Property_table_format):
    """
    Format for simulated absorption graphs
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

        
class IR_spectrum_property_format(Property_table_format):
    """
    Format for simulated IR spectra
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
    