from pathlib import Path
from datetime import datetime
from silico.misc.io import smkdir
import math

class Silico_directory():
    """
    Top level class for Directory helper classes.
    
    Sadly, we can't subclass pathlib Paths yet.
    """
    
    # Characters that are considered unsafe and are replaced with underscores.
    UNSAFE_CHARS = ["/", "\\"]

    def __init__(self, *args, **kwargs):
        """
         Constructor for Calculation_directory objects.
         """
        
        self.path = Path(*args, **kwargs)
        
    def __str__(self):
        return str(self.path)
    
    @classmethod
    def safe_name(self, file_name, unsafe_chars = None):
        """
        Get a safe version of a file name.
        
        Note that both '/' and '\' are considered unsafe, so supplying paths to this function is probably not what you want to do.
        
        :param file_name: The unsafe file_name.
        :param unsafe_chars: A list of characters to replace. If None is supplied, the default list is used.
        """
        if unsafe_chars is None:
            unsafe_chars = self.UNSAFE_CHARS
        
        safe_str = file_name
        for unsafe_char in self.UNSAFE_CHARS:
            safe_str = safe_str.replace(unsafe_char, "_")
        return safe_str

class Molecule_directory(Silico_directory):
    """
    Class that represents directory holding several calculations for a single molecule.
    """        
        
    @classmethod
    def from_calculation(self, calculation):
        """
        Create a Molecule_directory object from a Calculation_target.
        """
        return self(calculation.output, calculation.molecule_name)
    
    def create_structure(self):
        """
        Create the directory structure of this Molecule_directory object.
        
        :return: True normally, False if the structure already exists.
        """
        try:
            self.path.mkdir(parents = True)
        except FileExistsError:
            return False
        return True
        
    def get_calculation_dirs(self):
        """
        Get a list of Calculation_directory objects that are in this molecule dir.
        """
        return [Calculation_directory(calc_path) for calc_path in self.path.iterdir()]

class Calculation_directory(Silico_directory):
    """
    Class that represents the directory hierarchy where we submit calculations to.
    """
    
    def __init__(self, molecule_directory, program, name, *, directory_time = None, create = False, program_sub_folder = False, prepend_program_name = True, append_program_name = False):
        """
        Constructor for calculation directory objects.
        
        :param molecule_directory: Molecule_directory object in which the calculation is to take place.
        :param program: Name of the program carrying out the calculation. This can be None if program_sub_folder, prepend_program_name and append_program_name are all False.
        :param name: Name of the calculation.
        :param create: Whether to create the directory structure now.
        :param program_sub_folder: Whether to include a subdirectory for the program name.
        :param prepend_program_name: Whether to prepend the program name to the calculation directory.
        :param append_program_name: Whether to append the program name to the calculation directory.
        """
        self.molecule_directory = molecule_directory
        self.program = program
        self.name = name
        if prepend_program_name:
            self.name = "{} {}".format(self.program, self.name)
        if append_program_name:
            self.name = "{} {}".format(self.name, self.program)
        self.program_sub_folder = program_sub_folder
        # Believe this is deprecated...
        self.directory_time = directory_time if directory_time is not None else datetime.now()
        if create:
            self.molecule_directory.create_structure()
            self.create_structure(True)
            
    @classmethod
    def from_calculation(self, calculation, create = False):
        """
        Create a Calculation_directory object from a Calculation_target object.
        """
        return self(
            Molecule_directory.from_calculation(calculation),
            calculation.program.name,
            self.safe_name(calculation.name),
            create = create,
            program_sub_folder = calculation.structure['program_sub_folder'],
            prepend_program_name = calculation.structure['prepend_program_name'],
            append_program_name = calculation.structure['append_program_name']
        )        
    
    @property
    def path(self):
        """
        Pathlib Path object to this calculation dir.
        """
        if self.program_sub_folder:
            return Path(str(self.molecule_directory), self.program, self.name)
        else:
            return Path(str(self.molecule_directory), self.name)
        
    @path.setter
    def path(self, value):
        """
        """
        self.name = value.name
    
    @property
    def input_directory(self):
        """
        Full path to the calculation input directory.
        """
        return Path(str(self) + "/Input")
    
    @property
    def prep_directory(self):
        """
        Full path to the calculation prep directory.
        
        This directory is similar to the input directory; it is used for temporary input and setup.
        """
        return Path(str(self) + "/Setup")
    
    @property
    def output_directory(self):
        """
        Full path to the calculation output directory.
        """
        return Path(str(self) + "/Output")
    
    @property
    def result_directory(self):
        """
        Full path to the results directory.
        """
        return Path(str(self) + "/Results")
    
    @property
    def scratch_directory(self):
        """
        Full path to the scratch directory.
        
        Note that this directory is not the scratch directory as written to by quantum chemistry programs, but rather where we attempt to save the scratch in cases of error etc.
        """
        return Path(str(self) + "/Scratch")
    
    @property
    def flag_directory(self):
        """
        Full path to the flags directory.
        
        Flags are empty text files where the name of the file conveys information about the state of the calculation (STARTED, SUCCESS, ERROR etc)
        """
        return Path(str(self) + "/Flags")
    
    @property
    def log_directory(self):
        """
        Full path to the logging directory.
        """
        return Path(str(self) + "/Logs")
        
        
    def set_flag(self, flag, safe = True):
        """
        Set a flag file.
        
        Possible flags can be found in silico.submit.flag.
        
        :param flag: The flag to set (a Flag enum member). If the flag is already set; this is a noop.
        :param safe: If True, this method will not raise exceptions.
        """
        try:
            Path(self.flag_directory, flag.name).touch()
        except Exception:
            if safe:
                # This is ok.
                pass
            else:
                raise
        
    def get_flag(self, flag, safe = True):
        """
        Determine whether a flag file has been set or not.
        
        Possible flags can be found in silico.submit.flag.
        
        :param flag: The flag to check (a Flag enum member).
        :param safe: If True, this method will not raise exceptions.
        :return: True if the flag has been set, False otherwise.
        """
        try:
            return Path(self.flag_directory, flag.name).exists()
        except Exception:
            if safe:
                # This is ok.
                pass
            else:
                raise
    
    def del_flag(self, flag):
        """
        Unset a flag file. Possible flags can be found in silico.submit.flag.
        
        Possible flags can be found in silico.submit.flag.
        
        :param flag: The flag to unset (a Flag enum member). If the flag is not already set; this is a noop.
        """
        try:
            Path(self.flag_directory, flag.name).unlink()
        except FileNotFoundError:
            # This is ok.
            pass
    
    @property
    def report_directory(self):
        """
        Full path to the calculation report directory.
        """
        return Path(str(self) + "/Report")
    
    @property
    def log_file(self):
        """
        Full path to the file to which info should be logged.
        """
        return Path(self.log_directory, "silico.out")
    
    def create_structure(self, adjust = False):
        """
        Create the directory structure of this Calculation_directory object.
        
        :return True normally, False if the directory structure has already been made.
        """
        retval = 0
        
        self.name = smkdir(self.path, 1 if not adjust else math.inf).name
        
        sub_dirs = [self.input_directory, self.output_directory, self.flag_directory, self.log_directory]
        
        for sub_dir in sub_dirs:
            try:
                sub_dir.mkdir(parents = True)
            except FileExistsError:
                retval += 1
        return retval != len(sub_dirs)
        
        
        