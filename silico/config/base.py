import deepmerge

class Config(dict):
    """
    Class that represents a collection of config options.
    """
    
    def __init__(self, *args, FILE_NAME = None, **kwargs):
        """
        Constructor for Configurable objects.
        
        Configurable objects are enhanced dicts, so this constructor takes the same args as dict does, with the following exceptions:
        
        :param FILE_NAME: The name of the config file from which this Configurable was loaded.
        """
        super().__init__(*args, **kwargs)
        self.FILE_NAME = FILE_NAME
                
    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            raise KeyError("Unknown/missing config option '{}'".format(key))

    @classmethod
    def merge_dict(self, new, old):
        """
        Recursively merge two dictionaries
         
        Any keys specified in new will overwrite those specified in old.
         
        :param new: A new dictionary to merge into the old.
        :param old: An old dictionary to be overwritten by new.
        """
        return deepmerge.always_merger.merge(old, new)
    
    def merge(self, new):
        """
        Recursively merge another dictionary into this config object.
        
        :param new: A new Config to merge into the old (this Config).
        """
        return self.merge_dict(new, self)
            