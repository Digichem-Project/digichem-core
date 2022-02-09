# General imports.
import os
import numpy
import yaml

# Silico imports.
from silico.config.configurable.option import Option
from silico.config.configurable.base import Configurable_class_target
import silico.logging
from silico.exception.configurable import Configurable_loader_exception


class Method_target(Configurable_class_target):
    """
    Top level class for user-configurable method targets (Calculation types, program types, destination types etc.)
    """
    
    hidden = Option(help = "If True, this method will not appear in lists (but can still be specified by the user). Useful for methods that should not be used naively.", type = bool, default = False, no_edit = True)
    warning = Option(help = "A warning message to display when this method is chosen.", default = None, type = str, no_edit = True)
            
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
    

def parse_method_from_file(file_name, method_library):
    """
    Parse a method from a given file name.
    
    The file is expected to be in yaml format and contain three definitions: destination, program and calculation,
    each of which corresponds to one of the parts of a method. Each definition can either be a string, which
    identifies an existing definition from the method_library, or a dict, which contains a new definition which
    will be parsed.
    
    :param file_name: Path to a file to parse.
    :param method_library: A loaded library fo method definitions.
    :returns: A tuple of (destination, program, calculation)
    """
    # First, parse the file with yaml.
    with open(file_name, "rt") as method_file:
        raw = yaml.safe_load(method_file)
    
    # Split into our three parts.
    try:
        raw_method = {'destination': raw.pop('destination'), 'program': raw.pop('program'), 'calculation': raw.pop('calculation')}
        
    except Exception as e:
        raise Exception("Unable to parse method file '{}'; missing one or more of 'destination', 'program' or 'calculation'".format(file_name)) from None
    
    # Give a warning if there's anything extra.
    if len(raw) > 0:
        silico.logging.get_logger().warning("Ignoring additional data found in method file '{}':\n{}".format(file_name, raw))
        
    method = {}
    
    # Now for each of the three parts, either fetch a definition from our config library, or else get a new definition from the options given.
    for method_part_name, raw_method_part in raw_method.items():
        if not isinstance(raw_method_part, dict):
            # The method has been given as a string; fetch the definition from our library.
            method_part = getattr(method_library, method_part_name + "s").resolve(raw_method_part)
            
        else:
            # The method has been given as a dict. Create a new definition.
            # First, get the class we're looking for.
            try:
                cls = Method_target.from_class_handle(method_part_name).from_class_handle(raw_method_part['class_name'])
            
            except KeyError:
                # No classname specified.
                raise Exception("Cannot parse '{}' definition from file '{}'; missing class_name option".format(method_part_name, file_name)) from None
            
            except ValueError:
                # Couldn't find given class.
                # This might not be the best exception
                raise Configurable_loader_exception(method_part, method_part_name, file_name, "class_name '{}' is not recognised".format(method_part['class_name']))
            
            # Now make from it (this will also validate).
            method_part = cls(file_name = file_name, validate = True, **raw_method_part)
            
            # And finalize so we can make real objects from it.
            method_part.finalize()
            
        method[method_part_name] = method_part
        
    return (method['destination'], method['program'], method['calculation'])    