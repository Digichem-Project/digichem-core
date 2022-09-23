import silico.file.types as file_types
from silico.file.base import File_converter, Dummy_file_maker
from silico.submit.destination.local import Series
from silico.submit.calculation.gaussian import make_NTO_calc
from silico.submit.memory import Memory
import tempfile
from silico.input.chk import Chk_input
from silico.exception.base import File_maker_exception
import shutil


class Chk_to_NTO_chk(File_converter):
    """
    A class for converting an existing gaussian checkpoint file to a new chk file that contains certain natural transition orbitals (NTOs).
    """
    
    # Text description of our input file type, used for error messages etc.
    #input_file_type = "chk"
    input_file_type = file_types.gaussian_chk_file
    # Text description of our output file type, used for error messages etc.
    #output_file_type = "fchk"
    output_file_type = file_types.gaussian_NTO_chk_file
    
    def __init__(self, *args, chk_file, calc_t, prog_t, silico_options, **kwargs):
        """
        Constructor for Chk_to_fchk objects.
        
        :param output: Path to the file where the new chk file will be written.
        :param chk_file: Path to an existing Gaussian .chk file from which NTOs will be extracted.
        :param calc_t: The calculation template that will be run to calculate the new density.
        :param prog_t: The program template that will be run to calculate the new density.
        :param silico_options: Global silico options.
        """
        super().__init__(*args, input_file = chk_file, **kwargs)
        
        self.silico_options = silico_options
        
        # Save our calculation, program and destination templates.
        # We use an in series destination so we will block while the calc runs.
        self.destination_t = Series(
            name = "NTO chk"
        )
        self.destination_t.finalize()
            
        # Given calc program.
        self.prog_t = prog_t
        
        # Given calc.
        self.calc_t = calc_t
        
    @classmethod
    def from_options(self, output, *args, chk_file, transition, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        
        :param output: Path to the file where the new chk file will be written.
        :param chk_file: Path to an existing Gaussian .chk file from which NTOs will be extracted.
        :param transition: The integer number of the transition to calculate NTOs for.
        :param options: Global silico options.
        """
        # First, get our program.
        prog_t = options['report']['gaussian']['program']
        
        # Give up if no program available.
        if prog_t is None:
            return Dummy_file_maker(output, "No program definition is available (set the report: gaussian: program: option)")
        
        # Next, generate our calculation.
        calc_t = make_NTO_calc(
            name = chk_file,
            memory = Memory(options['report']['gaussian']['memory']),
            num_CPUs = options['report']['gaussian']['num_CPUs'],
            transition = transition,
            scratch_path = options['report']['gaussian']['scratch_path']
        )
        
        # And continue.
        return self(
            output,
            *args,
            chk_file = chk_file,
            calc_t = calc_t.inner_cls,
            prog_t = prog_t.inner_cls,
            silico_options = options,
            **kwargs
        )
        
    @classmethod
    def from_calculation(self, *args, calculation, chk_file, transition, options, **kwargs):
        """
        Create a Turbomole cube maker from an existing Turbomole calculation.
        """
        calc_t = calculation.NTO_calc(transition = transition)
        
        return self(
            *args,
            chk_file = chk_file,
            calc_t = calc_t.inner_cls,
            prog_t = type(calculation.program),
            silico_options = options,
            **kwargs
        )
        
    def make_files(self):
        """
        Make the files referenced by this object.
        """
        # Link.
        destination = self.destination_t()
        prog = self.prog_t(destination)
        calc = self.calc_t(prog, global_silico_options = self.silico_options)
        
        # We'll write our calc to a tempdir.
        with tempfile.TemporaryDirectory() as tempdir:
            #calc.prepare_from_directory(tempdir, self.input_file, molecule_name = "Anadens", additional_files = [(self.first_density, self.first_density_file_name), (self.second_density, self.second_density_file_name)])
            calc.prepare(tempdir, Chk_input(self.input_file))
            
            # Go.
            try:
                calc.submit()
            except Exception as e:
                raise File_maker_exception(self, "Failed to make Gaussian NTO chk file") from e
            
            # Move the (hopefully) created output file to our real destination.
            try:
                src = prog.chk_file_path
                dst = self.output
                shutil.move(src, dst, copy_function = shutil.copy)
                
            except FileNotFoundError as e:
                raise File_maker_exception(self, "The requested chk file could not be found, perhaps the calculation was setup incorrectly?") from e
            
            