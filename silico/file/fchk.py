# General imports.
import subprocess
import os

# Silico imports.
from silico.file import File_converter
from silico.exception.base import File_maker_exception
import silico.file.types as file_types
import silico.logging
from silico.submit import Memory

class Chk_to_fchk(File_converter):
    """
    Class for creating Gaussian fchk files from Gaussian chk files.
    """
        
    # 'Path' to the formchk executable.
    formchk_executable = "formchk"
    
    # Text description of our input file type, used for error messages etc.
    #input_file_type = "chk"
    input_file_type = file_types.gaussian_chk_file
    # Text description of our output file type, used for error messages etc.
    #output_file_type = "fchk"
    output_file_type = file_types.gaussian_fchk_file
    
    def __init__(self, *args, chk_file = None, fchk_file = None, memory = None, **kwargs):
        """
        Constructor for Chk_to_fchk objects.
        
        See Image_maker for a full signature.
        
        :param output: The filename/path to the fchk file (this path doesn't need to point to a real file yet; we will use this path to write to).
        :param chk_file: Optional chk_file to use to generate this fchk file.
        :param fchk_file: An optional file path to an existing fchk file to use. If this is given (and points to an actual file), then a new fchk will not be made and this file will be used instead.
        :param memory: The amount of memory for formchk to use.
        """
        super().__init__(*args, input_file = chk_file, existing_file = fchk_file, **kwargs)
        memory = memory if memory is not None else "3 GB"
        self.memory = Memory(memory)
        
    def make_files(self):
        """
        Make the files referenced by this object.
        """
        # The signature we'll use to call formchk.
        signature = [
            "{}".format(self.formchk_executable),
            self.input_file.name,
            str(self.output.absolute())
        ]
        
        # Sadly (but not unexpectedly), the g09 version of formchk (and maybe other versions too) will crash if their are certain characters (at least '(' and ')') in the input directory name. The final part of the input name appears to be fine, as does everything in the output dir.
        # This is easily fixed by making the output dir absolute and cd'ing into the input directory.
                    
        try:
            formchk_proc =  subprocess.run(
                signature,
                # Capture both stdout and stderr.
                stdout = subprocess.PIPE,
                stderr = subprocess.STDOUT,
                universal_newlines = True,
                cwd = str(self.input_file.parent),
                env = dict(os.environ, GAUSS_MEMDEF = str(self.memory))
                )
        except FileNotFoundError:
            raise File_maker_exception(self, "Could not locate formchk executable '{}'".format(self.formchk_executable))
        
        # If something went wrong, dump output.
        if formchk_proc.returncode != 0:
            # An error occured.
            raise File_maker_exception(self, "Formchk did not exit successfully:\n{}".format(formchk_proc.stdout))
        else:
            # Everything appeared to go ok.
            # Dump formchk output if we're in debug.
            silico.logging.get_logger().debug(formchk_proc.stdout)
            