import collections
import math
from pathlib import Path
import tempfile
import os

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
        
        
        