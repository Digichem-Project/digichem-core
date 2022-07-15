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


def parse_method_from_file(file_name, method_library, resolve = True):
    """
    Parse a method from a given file name.
    
    The file is expected to be in yaml format and contain three definitions: destination, program and calculation,
    each of which corresponds to one of the parts of a method. Each definition can either be a string, which
    identifies an existing definition from the method_library, or a dict, which contains a new definition which
    will be parsed.
    
    :param file_name: Path to a file to parse.
    :param method_library: A loaded library fo method definitions.
    :param resolve: Whether to resolve method parts that are IDs to a method in the library. If False, the method ID will be returned instead.
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
            if resolve:
                try:
                    method_part = getattr(method_library, method_part_name + "s").resolve(raw_method_part)
                
                except Exception as e:
                    raise Exception("Cannot resolve {} ID '{}' from file '{}'".format(method_part_name, raw_method_part, file_name)) from e
                
            else:
                # Just save the id.
                method_part = raw_method_part
            
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
            method_part = cls(file_name = file_name, validate_now = True, **raw_method_part)
            
            # And finalize so we can make real objects from it.
            method_part.finalize()
            
        method[method_part_name] = method_part
        
    return (method['destination'], method['program'], method['calculation'])    