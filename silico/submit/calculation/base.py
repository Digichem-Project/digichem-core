import getpass
from pathlib import Path
import copy

from silico.submit.flag import Flag
from silico.exception.base import Submission_error
from silico.submit.memory import Memory
from silico.config.configurable.option import Option
from silico.config.configurable.options import Options
from silico.input import si_from_file
from silico.submit.base import Method_target
from silico.input.directory import Calculation_directory_input
from silico.submit.basis import BSE_basis_set
from silico.exception.configurable import Configurable_exception
from silico.submit.translate import Basis_set, Functional


class Calculation_target(Method_target):
    """
    Abstract top-level class for calculation targets.
    """
    # Top level Configurable for calculations.
    CLASS_HANDLE = ("calculation",)
    meta = Options(
        TYPE = Option(default = "calculation")
    )
    
            
    @classmethod
    def safe_name(self, file_name):
        """
        Get a filename safe for Gaussian (and other programs).
        
        What constitutes a safe name from Gaussian's point of view is not entirely clear, to play it safe we'll only allow alpha-numeric characters, dots and underscores.
        
        :param file_name: The file name to make safe, note that a path (containing path separators) will not maintain its structure after a call to safe_name().
        :return: The safe path name.
        """
        # Adapted from https://stackoverflow.com/questions/7406102/create-sane-safe-filename-from-any-unsafe-string
        safe_chars = "._"
        return "".join([char if char.isalnum() or char in safe_chars else "_" for char in file_name])
    
    def expand(self, calculations):
        """
        Expand this calculation target if it represents multiple real calcs.
        
        Objects should return a number (possibly 0) of calculation targets that will actually be submitted to; these targets need not be the object itself.
                
        :return: A list of ready-to-go calculation targets.
        """
        raise NotImplementedError("ABC Calculation_target cannot be expanded")

    @classmethod
    def link(self, methods, *, global_silico_options, prepare_only = False):
        """
        Prepare a number of Calculation_target objects for submission by creating an ordered, linked list.
        
        :param methods: A list of 3-membered tuples (destination, program, calculation) that are to be prepared.
        :param global_silico_options: A dict like object of default options to use for the calculations that will be linked. Note that the options given here can be overridden by specific options given to each calculation...
        :return: The method (a 3-membered tuple of (destination, program, calculation)) that is to be submitted first.
        """
        first = None
        previous = None
        
        # These objects are class templates.
        for destination_t, program_t, calculation_t in methods:
            # Expand calculation (because the 'calculation' could actually be a meta-calc representing multiple real calcs).
            for expanded_calculation_t in calculation_t.expand(global_silico_options.calculations):
                
                # Init the destination, prog and calc.
                # This also links the three together.
                destination = destination_t()
                prog = program_t(destination)
                calc = expanded_calculation_t(prog, global_silico_options = global_silico_options, prepare_only = prepare_only)
                
                # TODO: Need to ensure that the given calculation and program objects are compatible.
                # This is tricky because of circular imports,
                # probably means we should refactor our organisation of programs and calculations
                
                # If the calc was part of a series, set the series name.
                if "Series" in calculation_t.CLASS_HANDLE:
                    calc.series_name = calculation_t.combined_report_name
                
                # Keep track of the first.
                if first is None:
                    first = calc
                                    
                # Link (doubly)
                if previous is not None:
                    previous.next = calc
                    calc.previous = previous
                    
                previous = calc
        
        # Return the first calculation in the chain.
        return first


class AI_calculation_mixin():
    """
    Abstract mixin class for calculation types that are ab-initio (from first principles).
    """
    
    electron = Options(help = "Options for controlling electron occupancy (how many electrons there are per molecule, and how they fill orbitals).",
        multiplicity = Option(help = "Forcibly set the molecule multiplicity. Leave blank to use the multiplicity given in the coordinate file.", default = None, type = int),
        charge = Option(help = "Forcibly set the molecule charge. Leave blank to use the charge given in the coordinate file.", default = None, type = float),
        unrestricted = Option(help = "Whether to perform an unrestricted calculation", type = bool, default = False)
    )
    
    basis_set = Options(help = "The basis set to use.",
        internal = Option(help = "The name of a basis set native to the relevant calculation program.", type = Basis_set, exclude = "exchange"),
        exchange = Option(help = "The definition of a (number of) basis sets to use from the Basis Set Exchange (BSE), in the format 'basis set name': 'applicable elements' (for example: '6-31G(d,p)': '1,3-4,B-F').", type = BSE_basis_set, dump_func = lambda option, configurable, value: dict(value), exclude = "internal", edit_vtype = "dict", default = lambda option, configurable: BSE_basis_set())
    )
    
    @property
    def charge(self):
        """
        The molecule/system charge that we'll actually be using in the calculation.
        
        Unlike the electron['charge'] attribute, this property will translate None to the actual charge to be used.
        """
        return int(self.electron['charge'] if self.electron['charge'] is not None else self.input_coords.implicit_charge)
    
    @property
    def multiplicity(self):
        """
        The molecule/system multiplicity that we'll actually be using in the calculation.
        
        Unlike the electron['multiplicity'] attribute, this property will translate None to the actual multiplicity to be used.
        """
        return int(self.electron['multiplicity'] if self.electron['multiplicity'] is not None else self.input_coords.implicit_multiplicity)
    
    @property
    def unpaired_electrons(self):
        """
        The number of unpaired electrons at this given multiplicity.
        """
        return self.multiplicity -1
    
    @property
    def basis_set_name(self):
        """
        A descriptive label of the basis set being used for the calculation.
        """
        if self.basis_set['builtin'] is not None:
            return self.basis_set['builtin']
        
        elif self.basis_set['external'] is not None:
            return str(self.basis_set['external'])
        
        else:
            return None


#########################
# Validation Functions. #
#########################
def validate_method(option, owning_obj, value):
    """Function for validating method choices."""
    # First check exactly one method has been selected.
    found = None
    for method_type in ("hf", "dft", "mp", "cc"):
        if value[method_type]['calc']:
            if found is not None:
                # A method has already been chosen.
                raise Configurable_exception(owning_obj, "The '{}' method cannot be selected at the same time as the '{}' method".format(found, method_type))
            
            else:
                found = method_type
                
    if found is None:
        # No method given.
        raise Configurable_exception(owning_obj, "No calculation method (HF, DFT, MP, CC etc) has been chosen.")
    
    return True

def validate_dft(option, owning_obj, value):
    """Function for validating DFT choices."""
    if value['calc'] and value['functional'] is None:
        # DFT requested but no functional.
        raise Configurable_exception(owning_obj, "No DFT functional has been chosen {}: {}.".format(value['calc'], value['functional']))
    
    return True

def validate_properties(option, owning_obj, value):
    """Function for validating requested calculation properties."""
    if not any([value[prop_type]['calc'] for prop_type in ("sp", "opt", "freq", "es")]):
        raise Configurable_exception(owning_obj, "No calculation properties (SP, Opt, Freq, ES etc) have been chosen.")
    
    return True

def validate_solvent(option, owning_obj, value):
    """Function for validating solvent model."""
    if value['calc'] and value['solvent'] is None:
        raise Configurable_exception(owning_obj, "No solvent has been chosen.")
    
    return True


class Concrete_calculation(Calculation_target):
    """
    Top-level class for real calculations.
    """
    
    CLASS_HANDLE = ()
    DIRECTORY_CALCULATION = False
    
    # A list of strings describing the expected input file types (file extensions) for calculations of this class. The first item of this list will be passed to obabel via the -o flag. 
    INPUT_FILE_TYPES = []
    
    # Configurable options.
    performance = Options(
        memory = Option(help = "The amount of memory to use for the calculation. Typical memory suffixes, such as B (byte), KB (kilobyte), MB (megabyte), and GB (byte) are accepted.", required = True, type = Memory),
        num_cpu = Option("num_cpu", help = "An integer specifying the number of CPUs to use for the calculation.", default = 1, type = int)
    )
    
    # Method options, choses the type of calculation (DFT, MP2, CC etc).
    method = Options(help = "Options for controlling the computational method used. Only one method should be chosen at a time.", validate = validate_method,
        # Here, we use the 'calc' sub option to turn on or off each option.
        # This may seem unnecessary, but it allows the other options for each method to be set without turning them on,
        # so that they can be inherited for child options.
        hf = Options(help = "Options for the HF method.",
            calc = Option(help = "Whether to use the Hartree-Fock method", type = bool, default = False),
        ),
        dft = Options(help = "Options for the density functional theory (DFT) method.", validate = validate_dft,
            calc = Option(help = "Whether to use the DFT method.", type = bool, default = False),
            functional = Option(help = "The DFT functional to use.", type = Functional, default = None),
            dispersion = Option(help = "The empirical dispersion method to use, note that not all methods are suitable with all functionals.", choices = ("PFD", "GD2", "GD3", "GD3BJ", None), default = None),
            # TODO: Look at supporting grid for all calculations.
            grid = Option(help = "The size of the numerical integration grid.", default = None)
        ),
        mp = Options(help = "Options for the Møller–Plesset (MP) method.",
            calc = Option(help = "Whether to use the MP method.", type = bool, default = False),
            level = Option(help = "The MP order to use.", choices = ("MP2", "MP3", "MP4", "MP4(DQ)", "MP4(SDQ)", "MP5"), default = "MP2")
        ),
        cc = Options(help = "Options for the Coupled-Cluster (CC) method.",
            calc = Option(help = "Whether to use the CC method.", type = bool, default = False),
            level = Option(help = "The CC level to use.", type = str, choices = ("CCD", "CCSD", "CCSD(T)"), default = "CCD")
        )
    )
    
    scf = Options(help = "Options for controlling the self-consistent field (SCF) method.",
        # Turbomole does not seem to offer different SCF procedures, where as Gaussian (and others?) do,
        # so options here will be program specfic.
        method = Option(help = "The SCF method or procedure to use.", default = None),
        iterations = Option(help = "The maximum number of SCF iterations to permit before aborting. If not specified, program defaults will be used.", type = int, default = None),
        convergence = Option(help = "The SCF convergence criteria, expressed as 10^n (where n is the value of this option).", type = int, default = None),
        damping = Options(help = "Options for damping SCF oscillations.",
            calc = Option(help = "Whether to enable SCF damping", type = bool, default = False),
            # Other options here are largely specific to each calc program.
        )
    )
    
    # Properties options, selects what to calculate (excited states, frequencies etc).
    properties = Options(help = "Options for controlling which properties to calculate.", validate = validate_properties,
        sp = Options(help = "Options for calculating single-point energies.",
            calc = Option(help = "Whether to calculate single-point energies.", type = bool, default = False)
        ),
        opt = Options(help = "Options for calculating optimised geometries.",
            calc = Option(help = "Whether to calculate optimised geometries." ,type = bool, default = False),
            iterations = Option(help = "The maximum number of optimisation iterations to permit before aborting. If not specified, program defaults will be used.", type = int, default = None),
        ),
        freq = Options(help = "Options for calculating vibrational frequencies.",
            calc = Option(help = "Whether to calculate vibrational frequencies.", type = bool, default = False),
            numerical = Option(help = "Whether to calculate vibrational frequencies numerically (as opposed to analytically). Numerical calcualtion is generally inferior, but is usually more broadly available. If analytical frequencies are not available, most programs will switch to numerical automatically.", type = bool, default = False)
        ),
        # TODO: Implement.
#         nmr = Options(help = "Options for calculating nuclear-magnetic resonance (NMR) shielding tensors and magnetic susceptibilities.",
#             calc = Option(help = "Whether to calculate NMR.", type = bool, default = False),
#         ),
        es = Options(help = "Options for calculating excited states.",
            calc = Option(help = "Whether to calculate excited states.", type = bool, default = False),
            method = Option(help = "The method to use to calculate excited states.", choices = ("CIS", "CIS(D)", "TD-DFT", "TDA", "EOMCCSD"), default = "CIS"),
            multiplicity = Option(help = "The multiplicity of the excited states to calculate. This option is only meaningful for singlet ground state molecules.", choices = ("Singlet", "Triplet", "50-50"), default = "Singlet"),
            num_states = Option(help = "The number of excited states to calculate. If 50-50 is given as the multiplicity, this is the number of each multiplicity to calculate.", type = int, default = 10),
            state_of_interest = Option(help = "The principal excited state of interest (SOS). The exact meaning and use of the SOS depends on the calculation being performed, but it is used, for example, to select which excited state to optimise (for ES and Opt calcs).", type = int, default = 1),
        )
    )
    
    # Options for controlling simulated solvation.
    # solvation might be a better name for this option?
    solution = Options(help = "Options for using a simulated solvent environment.", validate = validate_solvent,
        calc = Option(help = "Whether to use a solution model.", type = bool, default = False),
        model = Option(help = "The solvent model to use.", choices = ("PCM", "CPCM", "Dipole", "IPCM", "SCIPCM", "SMD"), default = "PCM"),
        solvent = Option(help = "The name of the solvent to use.", default = None)
    )
    
    scratch_options = Options(help = "Options that control the use of the scratch directory",
        use_scratch = Option(help = "Whether to use a scratch directory. False will disable the scratch directory, and is not recommended", default = True, type = bool),
        path = Option(help = "Path to the top of the scratch directory. For each calculation, a sub-directory will be created under this path. Note that some calculation programs (Gaussian16 at least) cannot handle whitespace in this path.", default = "/scratch", type = str),
        use_username = Option(help = "Whether to create a subdirectory for each user", default = True, type = bool),
        keep = Option(help = "Whether to copy any leftover files from the scratch directory once the calculation has completed successfully", default = False, type = bool),
        rescue = Option(help = "Whether to copy files from the scratch directory if the calculation fails or stops unexpectedly", default = True, type = bool),
        force_delete = Option(help = "Whether to always delete the scratch directory at the end of the calculation, even if output files could not be safely copied", default = False, type = bool),
        all_output = Option(help = "Whether to output all files to the scratch directory. If False, only intermediate files will be written to scratch (.log will be written directly to output directory for example)", default = False, type = bool)
    )
    
    structure = Options(help = "Options that control the calculation folder structure",
        program_sub_folder = Option(help = "Whether to use a separate subfolder for each calculation program", default = False, type = bool),
        prepend_program_name = Option(help = "Whether to prepend the calculation program name to the calculation folder", default = True, type = bool),
        append_program_name = Option(help = "Whether to append the calculation program name to the calculation folder", default = False, type = bool),
    )
    
    post = Options(help = "Options that control post processing of the calculation results.",
        write_summary = Option(help = "Whether to write Silico summary text files to the 'Results' folder at the end of the calculation", default = True, type = bool),
        write_report = Option(help = "Whether to write a Silico PDF report to the 'Report' folder at the end of the calculation", default = True, type = bool)
    )
    
    custom_silico_options = Option(
        "silico_options",
        help = "Silico options to overwrite for this calculation",
        default = {},
        type = dict
    )
        
    @property
    def mem_per_cpu(self):
        """
        Get the amount of memory to assign (per CPU).
        """
        return Memory(float(self.performance['memory']) / self.num_cpu)
    
    @property
    def num_cpu(self):
        """
        The number of CPUs to use for the calculation, this (roughly ?) translates to the number of worker processes that will be used.
        
        This property will translate 'auto' into a number of CPUs, use _num_CPUs if this is not desirable.
        """
        if self.performance['num_cpu'] == "auto":
            return self.get_available_CPUs()
        else:
            return self.performance['num_cpu']
    
    @num_cpu.setter
    def num_cpu(self, value):
        """
        Set the number of CPUs to use for the calculation. In addition to an exact integer amount, the string "auto" can also be supplied, in which case all available CPUs will be used.
        """
        # Set.
        self._num_CPUs = value
    
    def expand(self, calculations):
        """
        Expand this calculation target if it represents multiple real calcs.
        
        Objects should return a number (possibly 0) of calculation targets that will actually be submitted to; these targets need not be the object itself.
        
        This default implementation simply returns the same object.
        
        :param calculations: A Configurable_list of other calculations that this calculation could expand into.
        :return: A list of ready-to-go calculation targets.
        """
        return (self,)
            
    @property
    def is_unrestricted(self):
        """
        Determine whether this calculation will be using an unrestricted method.
        """
        # If we've been told to be unrestricted, use that.
        if self.electron['unrestricted']:
            return True
        
        # Otherwise, we can guess based on our multiplicity.
        return self.multiplicity == 1
    
    
    ############################
    # Class creation mechanism #
    ############################
    
    class _actual(Calculation_target._actual):
        """
        Inner class for calculations.
        """
        
        def __init__(self, program, *, global_silico_options, prepare_only = False):
            """
            Constructor for calculation objects.
            
            :param program: A Program_target_actual object to submit to.
            :param global_silico_options: A dict like object of options to use for this calculation. Note that the options given here can be overridden by specific options given to this calculation...
            :param prepare_only: Whether to only prepare this calculation and not perform it. 
            """    
            # Next is a linked list of Calculation_targets. We will call submit() on next once this object has been submitted.
            self.next = None
            self.previous = None
            self.output = None
            self.input_coords = None
            self.program = program
            self.global_silico_options = global_silico_options
            self.program.calculation = self
            # If this calculation was chosen as part of a series (meta-calc), this is the name of that series.
            # If this calc was chosen as an individual, this will be None.
            self.series_name = None
            self.prepare_only = prepare_only

        @property
        def silico_options(self):
            """
            Get the specific silico options to this calculation.
            
            This property is a speed-hack; it combines global and specific silico_options the first time it is called and caches the result.
            """
            try:
                return self._silico_options
            except AttributeError:
                # First time.
                # Deep copy silico_options (because we're going to merge it).
                self._silico_options = copy.deepcopy(self.global_silico_options)
                
                # Merge silico_options with the global options.
                self._silico_options.deep_merge(self.custom_silico_options)
                            
                return self._silico_options
        
        @property
        def scratch_directory(self):
            """
            Path to the scratch directory to use for this calculation. This will return None if no scratch is to be used.
            """
            # Return None if we're not using scratch.
            if not self.scratch_options['use_scratch']:
                return None
            
            # Start with the root path.
            directory = self.scratch_options['path']
            
            # Add username if we've been asked to.
            if self.scratch_options['use_username']:
                directory += "/" + getpass.getuser()
            
            # Make a path and return.
            return Path(directory, self.program.destination.unique_name)
        
        @property
        def molecule_name(self):
            """
            The name of the molecule under study.
            
            This name is 'safe' for Gaussian and other sensitive programs.
            """
            return self.safe_name(self.input_coords.implicit_name)
        
        @property
        def descriptive_name(self):
            """
            Get a name that describes the calculation and file together.
            """
            return "{} {}".format(self.molecule_name, self.meta['name'])
            
        def prepare(self, output, input_coords, additional_files = None):
            """
            Prepare this calculation for submission.
            
            :param output: Path to a directory to perform the calculation in.
            :param input_coords: A Silico_coords object containing the coordinates on which the calculation will be performed.
            :param additional_files: An optional list of paths to additional files required for the calculation. These files will be copied into the output directory without modification. In addition, each 'path' can be a tuple where the first item is the path to an additional file and the second is a new file name (relative to the calculation dir) to copy to.
            """
            # Because we normally run the program in a different environment to where we are currently, we need to load all input files we need into memory so they can be pickled.
            self.output = output
            self.input_coords = input_coords
            # TODO: Improve handling of additional files.
            # In general, this implementation is fine, but it will fail if the pickled file is moved to another machine before resuming (for example), because the files will be left behind.
            # We don't currently support this functionality, but one day we might...
            additional_files = additional_files if additional_files is not None else []
            self.additional_files = []
            for additional_file in additional_files:
                # Files that have not been given an explicit destination name will have one auto set.
                if not isinstance(additional_file, tuple):
                    self.additional_files.append((additional_file, Path(additional_file).name))
                
                else:
                    self.additional_files.append(additional_file)
                
            
        def prepare_from_calculation(self, output, calculation, additional_files = None, directory = False):
            """
            Alternative submission constructor that gets the input coordinates from a previously finished calc.
            
            :param output: Path to a directory to perform the calculation in.
            :param calculation: A previously submitted calculation.
            :param additional_files: An optional list of paths to additional files required for the calculation. These files will be copied into the output directory without modification. In addition, each 'path' can be a tuple where the first item is the path to an additional file and the second is a new file name (relative to the calculation dir) to copy to.
            :param directory: If True, the entire Output directory of the given calculation is used as input, useful for restarting previous calculations etc. If False (the default), only the output coords of the given calculation are used.
            """
            if not directory:
                self.prepare_from_file(
                    output,
                    calculation.program.next_coords,
                    input_format = calculation.OUTPUT_COORD_TYPE,
                    molecule_name = calculation.molecule_name,
                    molecule_charge = calculation.input_coords.charge,
                    molecule_multiplicity = calculation.input_coords.multiplicity,
                    # Don't gen3D or add H (we want to use output coordinates exactly as-is).
                    gen3D = False,
                    additional_files = additional_files
                )
                
            else:
                self.prepare_from_directory(
                    output,
                    calculation.program.destination.calc_dir.output_directory,
                    molecule_name = calculation.molecule_name,
                    additional_files = additional_files
                )
            
        def prepare_from_file(self,
            output,
            input_file_path,
            *,
            input_format = None,
            gen3D = None,
            molecule_name = None,
            molecule_charge = None,
            molecule_multiplicity = None,
            additional_files = None):
            """
            Alternative submission constructor that first reads in an input file.
            
            :param output: Path to a directory to perform the calculation in.
            :param input_file_path: Path (string or pathlib.Path) to the file to submit. This file should contain input coordinates for the system under study.
            :param gen3D: If True (and convert is True or "auto") and the loaded molecule does not have 3D coordinates, obabel will be used to generate them.
            :param molecule_name: Optional molecule name to use for the new calculation. If None, a name will be chosen based on the given file.
            :param molecule_charge: optional molecule charge to use for the new calculation. If None, a charge will be taken from the given file or otherwise a default will be used.
            :param molecule_multiplicity: optional molecule multiplicity to use for the new calculation. If None, a multiplicity will be taken from the given file or otherwise a default will be used.
            :param additional_files: An optional list of paths to additional files required for the calculation. These files will be copied into the output directory without modification. In addition, each 'path' can be a tuple where the first item is the path to an additional file and the second is a new file name (relative to the calculation dir) to copy to.
            """    
            # First, try and convert our given input file to the universal silico input format.
            try:
                # Load file.
                input_coords = si_from_file(input_file_path, input_format, name = molecule_name, charge = molecule_charge, multiplicity = molecule_multiplicity, gen3D = gen3D)
            except Exception:
                raise Submission_error(self, "failed to prepare input file (input format may not be supported)", file_name = input_file_path)
            
            # Prep normally.
            self.prepare(output, input_coords, additional_files = additional_files)
            
        def prepare_from_directory(self, output, calculation_directory, *, molecule_name = None, additional_files = None):
            """
            Alternative submission constructor that uses an existing calculation directory as input.
            
            :param output: Path to a directory to perform the calculation in.
            :param calculation_directory: Path to an existing calculation directory to use as input.
            :param molecule_name: Optional molecule name to use for the new calculation. If None, a name will be chosen based on the given file.
            :param additional_files:  An optional list of paths to additional files required for the calculation. These files will be copied into the output directory without modification. In addition, each 'path' can be a tuple where the first item is the path to an additional file and the second is a new file name (relative to the calculation dir) to copy to.
            """
            # Get the input file.
            input_file = Calculation_directory_input(calculation_directory, molecule_name)
            
            # Prep normally.
            self.prepare(output, input_file, additional_files = additional_files)
            
        def submit(self):
            """
            Submit a calculation to this Calculation_target.
            """            
            # Watch for errors so we can set flags appropriately.
            try:                
                # Do submit.
                self.program.submit()
                
            except Exception:
                # Something went wrong, set the error flag.
                self.program.destination.calc_dir.set_flag(Flag.ERROR)
                self.program.destination.calc_dir.set_flag(Flag.DONE)
                raise
            
            # Can't put this in finally because it will catch Submission_paused...
            self.program.destination.calc_dir.set_flag(Flag.DONE)
            
            # We don't wrap our own submit_post() in flags because this method submits the next calculation; any errors that occur relate to the new submission, not this one which has just finished.        
            self.post()
            
        def post(self):
            """
            Step 4/4 of the submission process, this method is called after submission.
            
            If there are any more calculations in the chain; we will call submit() on the next here.
            """
            if self.next is not None and not self.prepare_only:
                # We have another calculation to do.
                # Submit our output file.
                try:                    
                    # Prepare.
                    #self.next.prepare_from_file(self.output, self.program.calc_output_file_path)
                    self.next.prepare_from_calculation(self.output, self)
                    
                    # Go.
                    self.next.submit()
                    
                except Exception:
                    raise Submission_error(self, "failed to submit to next calculation in chain")