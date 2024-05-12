import collections
import math
from pathlib import Path
import tempfile
import os
import shutil
import sys
import warnings

from digichem.datas import get_resource

def expand_path(pth):
    """
    Expand variables (both $VAR and '~') in a string.
    
    This function is similar to calling both os.path.expanduser() and os.path.expandvar(),
    but with some additional functionality for relative paths (which are interpreted relative
    to the silico data dir).
    """
    pth = str(pth)
    
    new_pth = pth.replace("$SILICO", str(get_resource("data")))
    if new_pth != pth:
        warnings.warn("the '$SILICO' magic variable is deprecated, use '$DIGICHEM' instead", DeprecationWarning)

    pth = new_pth
    pth = pth.replace("$DIGICHEM", str(get_resource("data")))
    pth = os.path.expanduser(pth)
    pth = os.path.expandvars(pth)
    return pth

def tail(file, num_lines = 20):
    """
    Return the last n lines of a file.
    """
    last_lines = collections.deque(maxlen = num_lines)
    
    # We'll read through the file from the top.
    # This is probably inefficient for huge files but is easy to implement and I don't care just now.
    for line in file:
        last_lines.append(line)
        
    return list(last_lines)


def smkdir(dir_name, max_iter = math.inf):
    """
    Safe mkdir.
    
    Will attempt to create the given directory. If a file with the given name already exists, a directory with the same name except with an appropriate number appended will be created instead.
    
    If max_iter == 1, this method acts just like a normal mkdir().
    
    :param dir_name: The path to a directory to try and create.
    :param max_iter: The maximum number of iterations to try and create the directory. Default is no max.
    :returns: The path to the directory actually created.
    """
    # First try and make our base directory.
    counter = 1
    dir_name = str(dir_name)
    directory = None
    while True:
        try:
            directory = Path(dir_name + " {}".format(str(counter).zfill(2)) if counter != 1 else dir_name)
            directory.mkdir(parents = True)
            break
        except FileExistsError:
            if counter < max_iter:
                counter +=1
            else:
                raise
            
    return directory


# TODO: Would be nice to implement this as a context manager that acts like a file (better compatibility with yaml.dump() for example).
def atomic_write(file, data):
    """
    Atomically write to a file.
    
    :param file: The name of file to write to.
    :param data: Data to write to the file.
    """
    file = Path(file)
    # Open a temp file in the same directory as the real file (so we're likely to be on the same file system).
    with tempfile.NamedTemporaryFile("wt", dir = file.parent, delete = False) as temp_write:
        # Write data to the temp file.
        temp_write.write(data)
        # Force it to disk.
        temp_write.flush()
        os.fsync(temp_write)
        # Rename the temp file over the real file (overwriting it).
        os.rename(temp_write.name, file)


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
        

def copytree(src, dst, symlinks = False, ignore = None, copy_function = shutil.copy):
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
            
def copy_contents(src, dst, copy_function = shutil.copy):
    """
    Copy the top-level contents of one directory to another.
    """
    Path(dst).mkdir(parents = True, exist_ok = True)
    for child in Path(src).iterdir():
        if child.is_file():
            copy_function(child, Path(dst, child.name))
            
            
def rmtree(path, *args, **kwargs):
    """
    Fixed implementation of shutil.rmtree that doesn't arbitrarily require the first path to be a directory.
    """
    try:
        os.remove(path)
        
    except IsADirectoryError:
        shutil.rmtree(path, *args, **kwargs)


class Multi_file_wrapper():
    """
    A class that can be used to transparently handle both 'normal' files and stdin/stdout.
    """
    
    def __init__(self, file, mode = "r", *args, **kwargs):
        # If the filename is the special symbol '-' (a dash), then 'open' stdout as appropriate).
        if file == '-':
            if 'w' in mode or 'a' in mode:
                # We are writing.
                if 'b' in mode:
                    # We are a binary format, use a different stdout that accepts binary output.
                    object.__setattr__(self, 'file', sys.stdout.buffer)
                    #self.file = sys.stdout.buffer
                else:
                    object.__setattr__(self, 'file', sys.stdout)
                    #self.file = sys.stdout
            elif '+' in mode:
                # We are updating (not allowed).
                raise ValueError("Invalid mode specified '{}', updating is not valid for stdin/stdout".format(mode))
            else:
                # We are reading.
                if 'b' in mode:
                    object.__setattr__(self, 'file', sys.stdin.buffer)
                    #self.file = sys.stdin.buffer
                else:
                    object.__setattr__(self, 'file', sys.stdin)
                    #self.file = sys.stdin
            # We don't close stdin/stdout.
            object.__setattr__(self, 'should_close_file', False)
            #self.should_close_file = False
        else:
            # Normal file.
            
            object.__setattr__(self, 'file', open(file, mode, *args, **kwargs))
            object.__setattr__(self, 'should_close_file', True)
            #self.file = open(file, mode, *args, **kwargs)
            #self.should_close_file = True
            
    def __getattr__(self, name):
        """
        Magic method so we can pretend to be a real file.
        """
        return getattr(self.file, name)
    
    def __setattr__(self, name, value):
        """
        Magic method so we can pretend to be a real file.
        """
        return setattr(self.file, name, value)
        
    def __delattr__(self, name):
        """
        Magic method so we can pretend to be a real file.
        """
        if name in self.__dict__:
            object.__delattr__(self, name)
        else:
            return delattr(self.file, name)
            
            
    def close(self):
        if self.should_close_file:
            self.file.close()
        else:
            self.file.flush()
            
    def __enter__(self):
        """
        Magic function for the 'with' keyword.
        """
        return self
        
    def __exit__(self, etype, value, traceback):
        """
        Magic function for the 'with' keyword, called at the end of the block.
        """
        # Close our file.
        self.close()
        