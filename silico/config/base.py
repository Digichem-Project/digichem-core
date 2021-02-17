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
    def merge_dict(self, new, old, *, none_to_old = False):
        """
        Recursively merge two dictionaries
        
        Any keys specified in new will overwrite those specified in old.
        
        :param new: A new dictionary to merge into the old.
        :param old: An old dictionary to be overwritten by new.
        :param none_to_old: If True, values of None in the new dictionary will be ignored (and so the value from old will be used). If False (the default), None will be set in the returned dictionary.
        """
        # This is inspired by https://stackoverflow.com/questions/823196/yaml-merge-in-python.
        # Only attempt to merge if both are dicts (otherwise we just return new straight away).
        if isinstance(new, dict) and isinstance(old, dict):
            for key, value in new.items():
                if key not in old:
                    # The option is new, so we can just set without danger.
                    old[key] = value
                else:
                    # Both old and new have the same key, so we need to merge.
                    old[key] = self.merge_dict(value, old[key], none_to_old = none_to_old)
            return old
        # This removed behaviour was intended for the main config file (it allows you to explicitly reset an option to the one in the default.yaml file) but it clashes with the submit config files (where setting None is designed to unset an option set by a PARENT)...
        elif new is None and none_to_old:
            return old
        else:
            return new
    
    def merge(self, new, *, none_to_old = False):
        """
        Recursively merge two Config objects.
        
        Any options specified in new will overwrite those specified in this Config.
        
        You probably want to deepcopy new and old before calling this method.
        
        :param new: A new Config to merge into the old (this Config).
        :param none_to_old: If True, values of None in the new dictionary will be ignored (and so the value from old will be used). If False (the default), None will be set in the returned dictionary.
        """
        return self.merge_dict(new, self, none_to_old = none_to_old)
            