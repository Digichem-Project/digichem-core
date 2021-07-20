from pathlib import Path

class Splitter():
    """
    Class for splitting configurable files into multiple configurable files.
    """
    
    def __init__(self, *parents):
        """
        Constructor for Splitter objects.
        
        :param parent: The 'parent' configurables to split. Each parent can either be a directory or one of the files within that directory.
        """
        self.parents = [Path(parent) for parent in parents]
        
    @classmethod
    def find_containing_directory(self, hint):
        """
        Find the first containing directory from a given file path hint.
        
        The first containing directory is either the hint itself, if the hint is a directory, or the containing directory of hint.
        
        :raises FileNotFoundError: If hint does not exist.
        :param hint: A file path to the find the first containing directory of.
        """
        hint = hint.resolve(strict = True)
        if hint.is_dir():
            return hint
        else:
            return hint.parent
        
    def split(self, option, values):
        """
        Split the configurables based on an option and a list of values.
        
        Each configurable will have a one new child created for each value in values.
        
        :param option: The name of an option to split upon.
        :param values: A list of values to split upon; each of newly created splits will have one of values set.
        """
        
    

