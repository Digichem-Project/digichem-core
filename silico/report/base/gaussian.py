# General imports.
from pathlib import Path
from logging import getLogger

# Silico imports.
from silico.report import Report
from silico.file.fchk import Chk_to_fchk
from silico.file.cube import Fchk_to_spin_cube, Fchk_to_cube
from silico.result.result import Result_set
from silico.exception import Unknown_file_type_exception, Silico_exception
import silico.file.types as file_types
import silico
from silico.parser import get_parser
from silico.result.alignment.base import Alignment

class Gaussian_report(Report):
    """
    A specialised report object for processing Gaussian results.
    """
    
    def __init__(self, result, *, chk_file_path = None, fchk_file_path = None, options):
        """
        Constructor for Gaussian reports.
        
        One (or both) of either chk_file_path or fchk_file_path must be given in order for orbital/molecular images to be rendered.
        
        :param result: A result_set object.
        :param chk_file_path: Optional path to a chk_file.
        :param fchk_file_path: Optional path to an fchk_file.
        """
        self.chk_file_path = chk_file_path
        self.fchk_file_path = fchk_file_path
        
        # Continue in parent.
        super().__init__(result, options = options)
        
    @classmethod
    def from_calculation_files(
            self,
            *input_files,
            gaussian_log_file = None,
            chk_file_path = None,
            fchk_file_path = None,
            discover_additional_inputs = True,
            alignment_class_name = None,
            options,
            **kwargs):
        """
        A more intelligent constructor that can automatically determine file type(s).
        
        Keyword arguments can be used to specify files with a given type to avoid guesswork.
        
        Note that all files given should be from the same calculation, else bizarre behaviour may occur.
        
        This method is designed to be called at the start of programs and so is quite verbose.
        
        :param *input_files:  Calculation result files to be analysed, file type will be determined automatically (from file extension in most cases).
        :param gaussian_log_file: A Gaussian .log file to analyse.
        :param chk_file_path: A Gaussian .chk (checkpoint) file to analyse. This file will be safely ignored if a .fchk file is also specified.
        :param fchk_file_path: A Gaussian .fchk (formatted checkpoint) file to analyse.
        :param discover_additional_inputs: If True, the input directory will be searched for missing input files.
        :param alignment_class_name: Optional string matching the class handle of an alignment class to use for analysis. If None, this will be interpreted from options.
        :param options: A silico Config dictionary which contains various options that control the appearance of this report.
        """
        # A dictionary of known files types.
        files = {
            file_types.gaussian_log_file: gaussian_log_file,
            file_types.gaussian_chk_file: chk_file_path,
            file_types.gaussian_fchk_file: fchk_file_path
        }
        
        # Go through our generic input_files and determine their type.
        for input_file in input_files:
            # Check type.
            key_name = None
            if file_types.gaussian_log_file.check(input_file):
                key_name = file_types.gaussian_log_file
            elif file_types.gaussian_chk_file.check(input_file):
                key_name = file_types.gaussian_chk_file
            elif file_types.gaussian_fchk_file.check(input_file):
                key_name = file_types.gaussian_fchk_file
            else:
                # Get upset if we don't recognise the given file.
                raise Unknown_file_type_exception(input_file, expected = "calculation result file")
            
            if files[key_name] is not None:
                # Give a warning if multiple of the same file type have been given.
                getLogger(silico.logger_name).warning("Ignoring input file '{}'; a {} file has already been specified: {}".format(input_file, key_name, files[key_name]))
            else:
                # Save the new file.
                files[key_name] = input_file
        
        # If we're allowed, try and get any missing inputs.
        if discover_additional_inputs:
            # Get a list of files that we actually know about.        
            concrete_files = [input_file for input_file in files.values() if input_file is not None]
            
            # See which inputs are missing.
            for file_type in files:
                if files[file_type] is None:
                    # Missing, see if we can find one.
                    files[file_type] = self.find_additional_inputs(concrete_files, file_type)
                    
                    # Print a message if we found something.
                    if files[file_type] is not None:
                        getLogger(silico.logger_name).info("Found '{}' in input directory; using as {} file".format(files[file_type], file_type))
        
        # Check we have enough files to actually continue.
        if files[file_types.gaussian_log_file] is None:
            raise Silico_exception("Missing required file type '{}'".format(file_types.gaussian_log_file))
        
        # Get a result set.
        # First decide on which alignment class we're using.
        alignment_class = Alignment.from_class_handle(options['alignment'] if alignment_class_name is None else alignment_class_name)
        result = get_parser(files[file_types.gaussian_log_file]).process(alignment_class)
        
        # Load emission results if we're given file names instead of results.
        for emission in ['vertical_emission_ground_result', 'adiabatic_emission_ground_result', 'emission_excited_result']:
            if kwargs.get(emission, None) is not None and not isinstance(kwargs.get(emission, None), Result_set):
                # This emission 'result' is not a result (assume it is a path); try and load it.
                try:
                    #kwargs[emission] = Result_set.from_calculation_file(kwargs[emission], alignment_class_name = alignment_class_name)
                    kwargs[emission] = get_parser(kwargs[emission]).process(alignment_class)
                except Exception:
                    raise Silico_exception("Error loading emission result file '{}'".format(kwargs[emission]))
        
        # Add emission energies to result.
        result.add_emission(
            **kwargs
        )
                
        # Use our proper constructor.
        report = self(result,chk_file_path = files[file_types.gaussian_chk_file], fchk_file_path = files[file_types.gaussian_fchk_file], options = options)

        # All done.
        return report
    
    def setup_cubes(self, output_dir, output_name):
        """
        Setup the cube files which will be used to render images.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        # First, get our fchk file (from which cubes are made in Gaussian.
        self.fchk_file = Chk_to_fchk(
            Path(output_dir, output_name + ".fchk"),
            chk_file = self.chk_file_path,
            fchk_file = self.fchk_file_path,
        )
        
        ################
        # Spin density #
        ################
        self.cubes['spin_density'] = Fchk_to_spin_cube.from_options(
            Path(output_dir, "Spin Density", output_name + ".spin.cube"),
            fchk_file = self.fchk_file,
            cubegen_type = "Spin",
            orbital = "SCF",
            options = self.options
        )
        
        
        ############
        # Orbitals #
        ############
        # We need to set images for both alpha and beta orbitals (if we have them).
        for orbital_list in (self.result.molecular_orbitals, self.result.beta_orbitals):
            for orbital in orbital_list:
                # First, decide what type of orbital we need.
                if orbital.spin_type == "alpha":
                    cubegen_type = "AMO"
                elif orbital.spin_type == "beta":
                    cubegen_type = "BMO"
                else:
                    cubegen_type = "MO"
                
                # Save cube.
                self.cubes[orbital.label] = Fchk_to_cube.from_options(
                    Path(output_dir, orbital.label, output_name + ".{}.cube".format(orbital.label)),
                    cubegen_type = cubegen_type,
                    orbital = orbital.level,
                    options = self.options)
        
        
        #############
        # Structure #
        #############
        # If we have an orbital cube, we can just reuse this for our structure.
        if "HOMO" in self.cubes:
            self.cubes['structure'] = self.cubes['HOMO']
        elif "HOMO (alpha)" in self.cubes:
            self.cubes['structure'] = self.cubes['HOMO (alpha)']
        else:
            # No MO cubes available, create one for structure.
            # We'll just use the HOMO to get our cube, as it almost certainly should exist.
            self.cubes['structure'] = Fchk_to_cube.from_options(
                Path(output_dir, "Structure", output_name + ".structure.cube"),
                cubegen_type = "MO",
                orbital = "HOMO",
                options = self.options)
    
    def cleanup(self):
        """
        Remove any intermediate files that may have been created by this object.
        """
        # Delete our fchk file.
        self.fchk_file.delete(lazy = True)
        
        # Continue.
        super().cleanup()
      
        