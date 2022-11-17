# Exceptions relating to Configurable objects.
import textwrap

from silico.exception.base import Silico_exception
from silico.config.configurable.util import getopt


class Configurable_exception(Silico_exception):
    
    def __init__(self, configurable, reason):
        """
        Constructor for Configurable_exception objects.
        
        :param configurable: The Configurable object in which the error occurred.
        :param reason: The reason why the error occurred.
        """
        self.configurable = configurable
        self.reason = reason


    def __str__(self, *args, **kwargs):
        """
        Get a string description of this error.
        """
        hierarchy = self.get_file_hierarchy()
        
        message = "Error in '{}'; {}".format(self.configurable.description, self.reason)
        
        if hierarchy is not None:
            message += "; loaded from file(s): {}".format(hierarchy)
        
        return message
    
    def get_file_hierarchy(self):
        """
        Get the file hierarchy from which this configurable was loaded.
        
        Each configurable can be loaded from more than one file, because later files can overwrite certain options given in earlier files etc.
        """
        try:
            hierarchy = [(loader.TAG, loader.file_name) for loader in self.configurable.loader_list if loader.TAG is not None]
            
        except AttributeError:
            # We have no file list.
            return None
        
        if len(hierarchy) == 0:
            # No hierarchy, return a single file name instead.
            return self.configurable.file_name
        
        elif len(hierarchy) == 1:
            return "{}: {}".format(hierarchy[0][0], hierarchy[0][1])
        
        else:
            return "\n" + textwrap.indent("\n".join("{}: {}".format(tag, file) for tag, file in hierarchy), "\t")


class Configurable_option_exception(Configurable_exception):
    
    def __init__(self, configurable, option, reason):
        """
        Constructor for Configurable_option_exception objects.
        """
        self.option = option
        super().__init__(configurable, "error in option '{}'; {}".format(option.full_name, reason))
        

class Disallowed_choice_exception(Configurable_option_exception):
    """
    Exception raised when an value has been set for a configurable option that is not one of the allowed choices.
    """
    
    def __init__(self, configurable, option, value):
        """
        Constructor.
        
        :param configurable: The configurable where the error has occured.
        :param option: The option with the disallowed choice.
        :param value: The disallowed value.
        """
        super().__init__(configurable, option, "value '{}' ({}) is not one of the allowed choices: {}".format(value, type(value), ", ".join(str(choice) for choice in option.choices)))


class Missing_option_exception(Configurable_option_exception, AttributeError):
    """
    Exception raised when a required option is not set in a Configurable object.
    """
    
    def __init__(self, configurable, option):
        """
        Constructor for Missing_option_exceptions.
        
        :param configurable: The Configurable object in which the error occurred.
        :param option: The option with the missing value.
        """
        super().__init__(configurable, option, "a value is required but has not been set")
        AttributeError.__init__(self)
    

class Configurable_loader_exception(Silico_exception):
    """
    Exceptions raised when reading and parsing configurable files.
    """
    
    def __init__(self, config, TYPE, file_name, reason = None):
        """
        :param config: The config dict loaded.
        :param TYPE: The TYPE of the configurable.
        :param file_name: The file name (if any) from which the config was loaded.
        :param reason: The reason the exception was raised.
        """
        self.config = config
        self.TYPE = TYPE
        self.file_name = file_name
        self.reason = reason
        
    def __str__(self):
        """
        Get a string description of this error.
        """
        message = "Error loading configurable"
        
        if getopt(self.config, "link", "tag", default = None) is not None:
            message += " '{}'".format(self.config['link']['tag'])
            
        message += " of type '{}'".format(self.TYPE)
        
        if self.file_name is not None:
            message += " from file '{}'".format(self.file_name)
            
        if self.reason is not None:
            message += "; {}".format(self.reason)
        
        return message
    
class Short_tag_path_error(Silico_exception):
    """
    A given TAG path was too short.
    """
    
    def __init__(self, tag_path):
        """
        Constructor for Short_tag_path_error
        
        :param tag_path: A list of TAG names that is too short.
        """
        self.tag_path = tag_path
        
    def __str__(self):
        """
        Get a string description of this error.
        """
        return "could not resolve TAG list '{}', too few TAG names given".format(self.tag_path)
    
class Long_tag_path_error(Silico_exception):
    """
    A given TAG path was too long.
    """
    
    def __init__(self, path, remaining_tags):
        """
        Constructor for Long_tag_path_error
        
        :param path: The loader path that has been built so far.
        :param remaining_tags: The tags that remain to be resolved.
        """
        self.path = path
        self.remaining_tags = remaining_tags
        
    def __str__(self):
        """
        Get a string description of this error.
        """
        return "could not resolve remaining items in TAG list '{}', already reached the single loader '{}'".format(self.remaining_tags, self.path[-1].TAG)
    
class Unresolvable_tag_path_error(Silico_exception):
    """
    A given TAG path is not resolvable (because it could refer to more than one configurable).
    """ 
    
    def __init__(self, tag, possible_loaders):
        """
        Constructor for Unresolvable_tag_path_error.
        
        :param tag: The duplicated TAG name.
        :param possible_paths: A list of configurable loaders that could match the given TAG.
        """
        self.tag = tag
        self.possible_loaders = possible_loaders

    def __str__(self):
        """
        Get a string description of this error.
        """
        msg = "Could not uniquely identify configurable by TAG: '{}'; there are multiple possible configurables that match:\n".format(self.tag)
        
        matching = ""
        for loader_path in self.possible_loaders:
            matching += " : ".join([loader.TAG for loader in loader_path if loader.TAG is not None]) + "\n"
        
        msg += textwrap.indent(matching, "\t")
            
        return msg

        