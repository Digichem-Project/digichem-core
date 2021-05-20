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
        return "Error in configurable '{}' from file '{}'; {}".format(self.configurable.description, self.configurable.file_name, self.reason)


class Configurable_option_exception(Configurable_exception):
    
    def __init__(self, configurable, option, reason):
        """
        Constructor for Configurable_option_exception objects.
        """
        super().__init__(configurable, "error in option '{}'; {}".format(option.name, reason))

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
        