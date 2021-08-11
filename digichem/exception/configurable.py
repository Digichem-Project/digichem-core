# Exceptions relating to Configurable objects.

from silico.exception import Silico_exception


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
        
        message = "Error in configurable '{}'".format(self.configurable.description)
        
        if hierarchy is not None:
            message += " from file(s): {}".format(hierarchy)
            
        message += ";\n{}".format(self.reason)
        return message
    
    def get_file_hierarchy(self):
        """
        Get the file hierarchy from which this configurable was loaded.
        
        Each configurable can be loaded from more than one file, because later files can overwrite certain options given in earlier files etc.
        """
        hierarchy = self.configurable.file_hierarchy
        
        if len(hierarchy) == 0:
            return None
        
        elif len(hierarchy) == 1:
            return "{}: {}".format(hierarchy[0][0], hierarchy[0][1])
        
        else:
            return "\n" + "\n".join("{}: {}".format(tag, file) for tag, file in hierarchy)


class Configurable_option_exception(Configurable_exception):
    
    def __init__(self, configurable, option, reason):
        """
        Constructor for Configurable_option_exception objects.
        """
        super().__init__(configurable, "error in option '{}'; {}".format(option.full_name, reason))


class Configurable_class_exception(Configurable_exception):
    """
    Exceptions occurring on the class of a Configurable (not the object).
    """
    
    @property
    def configurable_desc(self):
        """
        A string that describes the Configurable class that threw this error.
        """
        return self.configurable.__name__
    
    def __str__(self):
        """
        Get a string description of this error.
        """
        return "Error in '{}'; {}".format(self.configurable_desc, self.reason)
    

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
        
        if self.config.get("TAG", None) is not None:
            message += " '{}'".format(self.config['TAG'])
            
        message += " of type '{}'".format(self.TYPE)
        
        if self.file_name is not None:
            message += " from file '{}'".format(self.file_name)
            
        if self.reason is not None:
            message += "; {}".format(self.reason)
        
        return message


class Missing_option_exception(AttributeError, Configurable_option_exception):
    """
    Exception raised when a required option is not set in a Configurable object.
    """
    
    def __init__(self, configurable, option_name):
        """
        Constructor for Missing_option_exceptions.
        
        :param configurable: The Configurable object in which the error occurred.
        :param option: The name of the option which is missing.
        """
        super().__init__()
        self.configurable = configurable
        self.option_name = option_name
        
    def __str__(self):
        """
        Get a string description of this error.
        """
        return "Error in configurable '{}' from file '{}'; {}".format(self.configurable.description, self.configurable.file_name, "option '{}' is required but is not set".format(self.option_name))
        