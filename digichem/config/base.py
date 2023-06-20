import deepmerge
import yaml
from deepmerge import Merger

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
        # Taken from deepmerge docs: https://deepmerge.readthedocs.io/en/latest/
        merger = Merger(
            # pass in a list of tuple, with the
            # strategies you are looking to apply
            # to each type.
            [
                (list, ["override"]),
                (dict, ["merge"]),
                (set, ["union"])
            ],
            # next, choose the fallback strategies,
            # applied to all other types:
            ["override"],
            # finally, choose the strategies in
            # the case where the types conflict:
            ["override"]
        )
        return merger.merge(old, new)
    
    def merge(self, new):
        """
        Recursively merge another dictionary into this config object.
        
        :param new: A new Config to merge into the old (this Config).
        """
        return self.merge_dict(new, self)


class Auto_type():
    """
    A class that when called, converts a string into a more appropriate type automatically.
    
    The rules for deciding which type to convert to are the same as for parsing yaml using pyyaml, as that is the module that is relied on for conversion.
    """
    
    def __new__(cls, value):
        return yaml.safe_load(value)
        