def defres(value):
    """
    Evaluate a function argument that might be set to a default.
    
    If value is a Default object, the 'default' attribute of that object is returned. Otherwise, value is returned.
    """
    if isinstance(value, Default):
        return value.default
    
    else:
        return value

class Default():
    """
    A class used to signify that an argument given to a function is the default.
    """
    
    def __init__(self, default):
        """
        Constructor for Default objects.
        
        :param default: The real default value.
        """
        self.default = default