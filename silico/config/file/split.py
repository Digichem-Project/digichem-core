from silico.misc.directory import copytree

from pathlib import Path
import tempfile
import yaml

class Splitter():
    """
    Class for splitting configurable files into multiple configurable files.
    """
    
    def __init__(self, parent):
        """
        Constructor for Splitter objects.
        
        :param parent: The 'parent' configurable to split. Each parent can either be a directory or one of the files within that directory.
        """
        self.parent = self.find_containing_directory(Path(parent))
        self.tmp_dir = None
        
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
            raise NotADirectoryError("{} is not a directory".format(hint))
            #return hint.parent
        
    def split(self, option, values):
        """
        Split the configurables based on an option and a list of values.
        
        Each configurable will have a one new child created for each value in values.
        
        :param option: The name of an option to split upon.
        :param values: A list of values to split upon; each of newly created splits will have one of values set.
        """
        # First, move any directories (existing children) in our parent directory to a temporary location.
        # Create a temporary location.
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            
            # Move each child.
            for existing_child in self.parent.iterdir():
                if existing_child.is_dir():
                    existing_child.rename(Path(temp_dir, existing_child.name))
                    
            # Create a new subdirectory for each value we've been given.
            value_number = 1
            for value in values:
                # Decide on the name of the new child, combining the value and the number.
                new_child = Path(self.parent, "{:02} {}".format(value_number, value))
                
                # Make the dir.
                new_child.mkdir()
                
                # Now make a yaml file with the option and save in the newly made child.
                with open(Path(new_child, "01.yaml"), "wt") as child_file:
                    yaml.dump({option: value}, child_file)
                    
                # Finally, copy the original children of parent to the new child.
                for old_child in temp_dir.iterdir():
                    if old_child.is_dir():
                        copytree(old_child, Path(new_child, old_child.name))
                
                value_number += 1
                
            
            
            
        

