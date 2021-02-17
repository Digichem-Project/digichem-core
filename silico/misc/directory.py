from pathlib import Path
import os
import shutil

# Based on https://stackoverflow.com/questions/431684/how-do-i-change-the-working-directory-in-python/24176022#24176022
class cd:
    """
    Context manager for temporarily changing the working directory.
    """
    
    def __init__(self, directory):
        """
        Constructor for cd.
        
        :param directory: Path to the directory to change to.
        """
        self.new_directory = Path(directory)
        self.old_directory = None
        
    def __enter__(self):
        # Save our current working directory.
        self.old_directory = os.getcwd()
        
        # And change to the new directory.
        os.chdir(str(self.new_directory))
        
    def __exit__(self, exc_type, exc_value, traceback):
        # Restore our old directory.
        os.chdir(self.old_directory)
        

def copytree(src, dst, symlinks = False, ignore = None, copy_function = shutil.copy2):
    """
    Fixed implementation of shutil.copytree that doesn't arbitrarily fail if src exists.
    
    Adapted from https://stackoverflow.com/questions/1868714/how-do-i-copy-an-entire-directory-of-files-into-an-existing-directory-using-pyth
    """
    if not os.path.exists(dst):
        os.makedirs(dst)
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            copytree(s, d, symlinks = symlinks, ignore = ignore, copy_function = copy_function)
        else:
            copy_function(s, d)