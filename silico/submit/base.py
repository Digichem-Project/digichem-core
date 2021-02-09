import os
import numpy

from silico.config.configurable import Configurable
from silico.exception.configurable import Configurable_exception

class Configurable_target(Configurable):
    """
    Top level class for user-configurable submit targets (Calculation types, program types, method types etc.)
    """
    
    def get_children(self, target_list):
        """
        Get all the configurable target objects from a list that contain this object as a parent.
        """
        return [child for child in target_list if child.compatible_parent(self)]
    
    def compatible_parent(self, parent):
        """
        Determine whether a configurable is allowed to be a parent of this configurable.
        
        """
        return len(set(parent.NAMES).intersection(self.parents)) != 0
        
    def validate_parent(self, parent):
        """
        Determine whether a configurable is allowed to be a parent of this configurable, raising an exception if not.
        
        :raises Configurable_exception: If the given parent is not valid.
        :return: None.
        """
        if not self.compatible_parent(parent):
            raise Configurable_exception(self, "{} '{}' is not compatible with {} '{}'".format(parent.TYPE, parent.NAME, self.TYPE, self.NAME))
        
    @classmethod
    def get_available_CPUs(self):
        """
        Determine how many CPUs there are available to this process.
        
        It is generally not a good idea to rely on this function to determine how many CPUs to use for your calculation, not least because the environment that silico runs in is often very different to that the calculation is performed in. This function mainly exists as a fall-back.
        """
        try:
            return len(os.sched_getaffinity(0))
            # The python docs mention that sched_getaffinity() is not available on all systems, but it is not clear what exception will be raised in such a case. These two seem most likely.
        except (AttributeError, NotImplementedError):
            return os.cpu_count()
        
    
class Memory():
    """
    Class for storing memory-amounts.
    """
    
    # The units that we know about.
    UNITS = {
        'TB': 1000000000000,
        'GB': 1000000000,
        'MB': 1000000,
        'KB': 1000,
         'B': 1
        }
    
    # When outputting units, whether to separate the number and unit with a space.
    SPACE_UNIT = False
    
    def __init__(self, value = None, print_decimal = False, round_decimal = True):
        """
        """
        self.value = None
        self.auto = value
        self.print_decimal = print_decimal
        self.round_decimal = round_decimal
    
    @property
    def auto(self):
        """
        Get the amount of memory, as a string, using an automatically determined suffix.
        """
        # First, get an ordered list from our known units.
        ordered_units = sorted(self.UNITS.items(), key = lambda item: item[1])
        
        # Now see where our number best fits (thanks numpy).
        suffix_index = numpy.digitize((abs(self.value),), [ordered_unit[1] for ordered_unit in ordered_units])[0]
        
        # Wrap if we're out of bounds and convert to proper index.
        if suffix_index == 0:
            suffix_index += 0
        elif suffix_index == len(ordered_units):
            suffix_index -= 1
        else:
            suffix_index -= 1
        
        # Now check to see if we have a decimal part which is non zero.
        while not self.print_decimal and suffix_index > -1 and (self.value / ordered_units[suffix_index][1]) % 1 != 0:
            # Got a decimal part, use the next smallest unit.
            suffix_index -= 1
            
        if suffix_index < 0:
            # We ran out of units and couldn't remove our fraction.
            if self.round_decimal:
                # We're allowed to round our decimal away.
                suffix_index = 0
            else:
                raise ValueError("Unable to represent memory value '{}'B as non decimal".format(self.value))
        
        # Now get our value in the correct units.
        value = self.value / ordered_units[suffix_index][1]
        if not self.print_decimal:
            value = int(value)
        
        return "{}{}{}".format(value, " " if self.SPACE_UNIT else "", ordered_units[suffix_index][0])
            
        
    @auto.setter
    def auto(self, value):
        """
        Set the amount of memory, automatically determining the unit suffix (KB, MB, GB etc).
        """
        if isinstance(value, str):
            # First, determine the suffix (if any).
            suffix = None
            for unit, amount in  sorted(self.UNITS.items(), key = lambda item: item[1], reverse = True):
                if value.lower().endswith(unit.lower()):
                    suffix = unit
                    # This is a bit of a hack because all suffixes contain 'B' at the end.
                    break
                    
            # Now remove the suffix (if there is one), convert to float and multiply by the value of suffix (to get bytes).
            if suffix is not None:
                value = float(value[:-len(suffix)]) * self.UNITS[suffix]
        
        # Convert to int (because not sure it makes sense to represent fractions of bytes) and store.
        self.value = int(value)
    
    @property    
    def MiB(self):
        """
        The value of this memory in MiB.
        """
        return round(self.value/1048576, None)
    
    @property    
    def MB(self):
        """
        The value of this memory in MB.
        """
        return round(self.value/1000000, None)
        
    def __float__(self):
        """
        Floatify this memory amount (returning the number of bytes).
        """
        return float(self.value)
    
    def __int__(self):
        """
        Intify this memory amount (returning the number of bytes).
        
        In most cases, the int and float equivalents should be the same (because fractional bytes are not supported).
        """
        return int(self.value)
    
    def __str__(self):
        """
        Stringify this memory amount (returning the number of bytes) with an appropriate suffix.
        """
        return self.auto
        