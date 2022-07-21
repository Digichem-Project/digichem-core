# General imports.
from itertools import zip_longest
import warnings

# Silico imports.
from silico.result import Result_container
from silico.result import Result_object
from silico.result import Floatable_mixin, Unmergeable_container_mixin


class Vibrations_list(Result_container, Unmergeable_container_mixin):
    """
    Class for representing a group of molecular vibrations.
    """
    
    MERGE_WARNING = "Attempting to merge lists of vibrations that are not identical; non-equivalent vibrations will be ignored"
        
    @property
    def negative_frequencies(self):
        """
        Get a Vibrations_list object of the vibrations in this list that have negative frequencies (they are imaginary).
        """
        warnings.warn("negative_frequencies is deprecated, use negative instead")
        return self.negative
    
    @property
    def negative(self):
        """
        Get a Vibrations_list object of the vibrations in this list that have negative frequencies (they are imaginary).
        """
        return type(self)([vibration for vibration in self if vibration.frequency < 0])
    
    @classmethod
    def from_parser(self, parser):
        """
        Get an Vibrations_list object from an output file parser.
        
        :param parser: An output file parser.
        :return: A Vibrations_list object. The list will be empty if no vibration frequency data is available.
        """
        return self(Vibration.list_from_parser(parser))

    
class Vibration(Result_object, Floatable_mixin):
    """
    Class for representing vibrational frequencies.
    """
    
    def __init__(self, level, frequency, intensity, symmetry = None):
        """
        Vibration object constructor.
        
        :param level: The numerical index of the vibration.
        :param frequency: Frequency of the vibration (cm-1).
        :param intensity: Intensity of the vibration in IR (km mol-1)
        :param symmetry: Symmetry term of the vibration (string).
        """
        self.level = level
        self.symmetry = symmetry
        self.frequency = frequency
        self.intensity = intensity
        
    def __float__(self):
        """
        A float representation of this object.
        """
        return float(self.frequency)
        
    @classmethod
    def list_from_parser(self, parser):
        """
        Create a list of Vibration objects from an output file parser.
        
        :param parser: An output file parser.
        :return: A list of Vibration objects. The list will be empty if no frequency data is available.
        """
        try:
            return [self(index+1, frequency, intensity, symmetry) for index, (symmetry, frequency, intensity) in enumerate(zip_longest(getattr(parser.data, 'vibsyms', []), parser.data.vibfreqs, parser.data.vibirs, fillvalue = None))]
        except AttributeError:
            return []