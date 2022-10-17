# General imports.
from pathlib import Path
import shutil
from subprocess import TimeoutExpired, CalledProcessError
import os
import re
from mako.lookup import TemplateLookup
import subprocess
from itertools import chain
import warnings
import tempfile

# Silico imports.
from silico.submit.program.base import Program_target
from silico.config.configurable.option import Option
from silico.exception.base import Submission_error
from silico.misc.directory import copytree
from silico.parse import parse_calculation
import silico
from silico.input.directory import Calculation_directory_input
from silico.misc.io import tail


class Turbomole(Program_target):
    """
    Top level class for submitting calculations to Turbomole.
    """
    
    CLASS_HANDLE = ("Turbomole",) 

    define_executable = Option(help = "Name/path of the define executable", default = "define", type = str)
    mp2prep_executable = Option(help = "Name/path of the mp2prep executable", default = "mp2prep", type = str)
    root = Option(help = "Path to the directory where turbomole is installed", required = True, type = str)
    init_file = Option(help = "Path to the turbomole Config_turbo_env script which is run to set-up Turbomole", required = True, type = str)
    
    # Regex for matching the UFF section.
    UFF_SECTION = r"(\$uff\n(.|\n)+?)\n\$.+"
    
    ############################
    # Class creation mechanism #
    ############################
    
    class _actual(Program_target._actual):
        """
        Inner class for programs.
        """
        
        @property
        def environ(self):
            """
            Environmental variables for turbomole.
            """
            env = os.environ.copy()
            # Add the turbomole directory.
            env['TURBODIR'] = self.turbomole_root
            # Also extend path.
            
        
        @property
        def coord_file_path(self):
            """
            Path to the input coord file.
            """
            return Path(self.destination.calc_dir.input_directory, "coord")
        
        @property
        def control_file_path(self):
            """
            Path to the input control file.
            """
            return Path(self.destination.calc_dir.input_directory, "control")
        
        @property
        def define_output_path(self):
            """
            Path to the file in which define output is written.
            """
            return Path(self.destination.calc_dir.log_directory, "define.out")
        
        @property
        def turbomole_output_path(self):
            """
            Path to the file in which turbomole output is written.
            """
            return Path(self.destination.calc_dir.output_directory, self.calculation.molecule_name).with_suffix(".log")
        
        @property
        def calc_output_file_path(self):
            """
            Path to the main calculation output file.
            """
            return self.turbomole_output_path
                
        @property
        def next_coords(self):
            """
            Path to the output coordinate file that should be used for any subsequent calculations.
            """
            return Path(self.destination.calc_dir.output_directory, "coord")
        
        @property
        def scratch_base(self):
            """
            The basename of the directories to use as scratch.
            When running in multi-process mode (MPI), the scratch directory given to turbomole is only used as a suggestion; each node will create its own folder by APPENDING "-001" etc to the given file name.
            
            This is annoying, so when in MPI, we make sure each of these folders is created inside the specified scratch dir.
            """
            if self.calculation.scratch_directory is None:
                return None
            elif self.calculation.performance['parallel_mode'] == "MPI":
                return Path(self.calculation.scratch_directory, "node")
            else:
                return self.calculation.scratch_directory
                    
        @property
        def define_input_path(self):
            """
            Path to the input file used to power the define input generator.
            """
            return Path(self.destination.calc_dir.input_directory, "define.in")
        
        def define(self):
            """
            Run setup for turbomole.
            
            Normally this involves running define, but some turbomole calcs have a different setup (eg, UFF).
            """
            # IMPORTANT: The auxiliary basis set menus of define are broken in the following manner:
            # When attempting to assign an aux basis set, if the primary basis set cannot be found in the aux 
            # basis set library, define will panic and enter an uncontrollable state.
            # This is frustrating because it prevents us from ever assigning the correct aux basis set.
            #
            # This behaviour can be overcome if options like:
            #    $atoms
            #    n  1                                                                           \
            #       basis =n 6-31G*                                                             \
            #       jkbas =n cc-pVTZ
            #    ...
            # are already set in the control file, but this is challenging to do by hand.
            #
            # TODO: We need a solution to this.
            
            # Write define input file.
            with open(self.define_input_path, "wt") as define_input:
                define_input.write(self.calculation.define_input)
            
            # Get our wrapper script.
            wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/define.mako").render_unicode(program = self)
                        
            # Run control to generate input.
            try:
                subprocess.run(
                    ("bash",),
                    input = wrapper_body,
                    stdout = subprocess.PIPE,
                    stderr = subprocess.STDOUT,
                    universal_newlines = True,
                    timeout = self.calculation.performance['define_timeout'],
                    cwd = self.destination.calc_dir.prep_directory,
                    check = True
                )
                
            except TimeoutExpired as e:
                # Ran out of time, probably got stuck.
                raise Submission_error(self, "Program 'define' failed to finish executing in {} s, check output file '{}' for errors".format(self.calculation.define_timeout, self.define_output_path)) from e    
            
            except CalledProcessError as e:
                # Something went wrong.
                e.__context__ = None
                raise Submission_error(self, "define did not exit successfully:\n{}".format(e.stdout)) from e
            
            # Check the output from define to see if it ended happily or not.
            # Define uses (we think) the iso-8859-1 encoding.
            with open(self.define_output_path, "rt", encoding = "iso-8859-1") as define_output:
                if "define ended abnormally" in tail(define_output, 1)[0]:
                    raise Submission_error(self, "Program 'define' does not appear to have executed correctly, check output file '{}' for errors".format(self.define_output_path))
            
            # Sadly, some options are not supported by define and have to be appended manually.
            if self.calculation.properties['opt']['calc'] and self.calculation.properties['es']['calc']:
                self.add_control_option("$exopt {}".format(self.calculation.properties['es']['state_of_interest']))
                
            if self.calculation.analysis['anadens']['calculate']:
                self.add_control_option("$anadens\n calc {} from\n 1d0 {}\n {}1d0 {}".format(
                    self.calculation.analysis['anadens']['output'],
                    self.calculation.analysis['anadens']['first_density'],
                    self.calculation.analysis['anadens']['operator'],
                    self.calculation.analysis['anadens']['second_density'],
            ))
        
        
        def wrap_unsafe_module(self, run_func):
            """
            Wrap a function so that the current Output dir is first copied to a 'safe' location (one without whitespace in the path etc)
            and the function is called with the the new safe location as an argument.
            
            :param run_func: The function to call. The path to the temp dir is passed as an argument.
            """
            # First, get ourselves a temp directory.
            # NOTE: This is not actually guaranteed to be safe, because the temp directory
            # could be changed by the user. We don't handle this eventuality for now.
            with tempfile.TemporaryDirectory() as temp_dir:
                # Copy the current output dir to a 'safe' location.
                copytree(self.destination.calc_dir.output_directory, temp_dir)
                
                # Delete the main output file from this copied dir
                # (because it's not actually going to be written to, and otherise
                # it will overwrite our actual log file when we copy back(.
                Path(temp_dir, self.turbomole_output_path.name).unlink(True)
                
                # Now run the program.
                try:
                    run_func(temp_dir)
                
                finally:
                    # Copy back the temp dir so we can see what happened.
                    copytree(temp_dir, self.destination.calc_dir.output_directory)
        
            
        def add_control_option(self, option):
            """
            Manually write to the control file.
            
            The control file should already exist before this method is called.
            """
            # TODO: Turbomole includes the program adg (add data group, see manual) which could be used for this purpose?
            control_name = Path(self.destination.calc_dir.prep_directory, "control")
            tmp_control_name = Path(self.destination.calc_dir.prep_directory, "control.tmp")
            done = False
            
            with open(control_name, "r") as control_file_in, open(tmp_control_name, "w") as control_file_out:
                for line in control_file_in:
                    if "$end" in line.lower() and not done:
                        # This is the end, append option now.
                        control_file_out.write(option + "\n")
                        done = True
                    
                    # Write existing line.
                    control_file_out.write(line)
            
            tmp_control_name.rename(control_name)
        
            
        def UFF_define(self):
            """
            Run setup for turbomole UFF.
            
            Unlike real turbomole calcs, UFF doesn't use define for setup (but we still take some options from the control file weirdly...)
            """
            # The turbomole program mdprep could help here?
            # Get our wrapper script.
            wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/module.mako").render_unicode(program = self, module = self.calculation.modules[0])
            
            # Copy atoms to prep dir.
            shutil.copy2(self.coord_file_path, self.destination.calc_dir.prep_directory)
            
            # Run control to generate input.
            # We check for errors, but sadly uff doesn't seem to raise any?
            try:
                subprocess.run(
                    ("bash",),
                    input = wrapper_body,
                    stdout = subprocess.PIPE,
                    stderr = subprocess.STDOUT,
                    universal_newlines = True,
                    cwd = self.destination.calc_dir.prep_directory,
                    check = True
                )
                
            except CalledProcessError as e:
                # Something went wrong.
                e.__context__ = None
                raise Submission_error(self, "UFF (setup) did not exit successfully:\n{}".format(e.stdout)) from e
            
            except Exception as e:
                # Something else.
                raise e from None
            
            # Get our custom uff section.
            uff_input = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/uff.mako").render_unicode(program = self)
            
            # We now need to manually alter some options in control (sigh...).
            # Open control.
            with open(Path(self.destination.calc_dir.prep_directory, "control"), "r+") as control_file:
                # Read in existing control.
                control = control_file.read()
                
                # Replace existing uff section with our custom one.
                control = re.sub(self.UFF_SECTION, uff_input, control)
                
                # Seek back to start of file.
                control_file.seek(0)
                
                # Write modified control.
                control_file.write(control)
                
                # Remove any leftover data.
                control_file.truncate()
            
        
        def pre(self):
            """
            Pre-calculation setup for Turbomole.
            """
            # Call parent for setup first.
            super().pre()
            
            # Sanity check the $PATH variable.
            if " " in os.environ['PATH']:
                warnings.warn("The $PATH environmental variable contains whitespace. This is not supported by the dscf Turbomole module and will likely lead to calculation failure.")
            
            # If we have a directory calc, copy the old directory to our new Prep directory.
            if isinstance(self.calculation.input_coords, Calculation_directory_input):
                copytree(self.calculation.input_coords.calculation_directory, self.destination.calc_dir.prep_directory)
                
                # Clean up some old files which may be left by the previous calc if run with silico.
                for delete_file in chain(
                    (Path(self.destination.calc_dir.prep_directory, self.define_output_path.name), Path(self.destination.calc_dir.prep_directory, self.destination.calc_dir.log_file.name)),
                    self.destination.calc_dir.prep_directory.glob("*.log"),
                    self.destination.calc_dir.prep_directory.glob("job.*"),
                ):
                    
                    try:
                        delete_file.unlink()
                    except FileNotFoundError:
                        # This is ok.
                        pass
                    
                    
            else:
                # Our calc is not a directory calc, write our input file to our calculation Input directory.
                with open(self.coord_file_path, "wt") as coord_file:
                    coord_file.write(self.calculation.input_coords.to_format("tmol"))
                    
                # Create our prep dir.
                self.destination.calc_dir.prep_directory.mkdir()
                            
            # Run define.
            if "Turbomole-UFF" in self.calculation.CLASS_HANDLE:
                # Use alternative UFF setup.
                self.UFF_define()
            else:
                # Normal setup.
                self.define()
            
            # If we're using a scratch output dir, create it now.
            if self.scratch_output is not None:
                try:        
                    # Copy our input dir to the scratch version.
                    copytree(self.destination.calc_dir.prep_directory, self.scratch_output)
                except Exception as e:
                    raise Submission_error(self, "Failed to copy input to scratch subdirectory") from e
                
            else:
                # Copy input to normal output folder.
                copytree(self.destination.calc_dir.prep_directory, self.destination.calc_dir.output_directory)
                
            # Copy prep dir to input dir and delete.
            copytree(self.destination.calc_dir.prep_directory, self.destination.calc_dir.input_directory)
            shutil.rmtree(self.destination.calc_dir.prep_directory)
            
        def calculate(self):
            """
            Main submission method; the calculation will be run here (for which this method will block, possibly for hours+).
            """
            # Turbomole is not, in fact, a single program. Instead, it is a whole suite of programs and scripts that together
            # perform the calculation(s) requested. This is much more complicated than most other CC programs, because it means
            # we have to not only worry about setting the correct options in the control file (the input file), but also
            # about calling the correct programs in the correct order. Each of these sub-programs is called a 'module'.
            # Most modules can be setup before any of the calculation starts, but this is not true of mpgrad. This module
            # requires setup via the mp2prep script, but this script can only be called after the dscf module has already run.
            #
            # To run a Turbomole calculation then, we iteratively call each specified module, wrapping each in the necessary
            # setup scripts.
            modules = self.calculation.modules
            for index, module in enumerate(modules):
                # Start by announcing which module we're running.
                silico.log.get_logger().info("Running module {}/{}: {}".format(index+1, len(modules), module.name))
                
                # Get our wrapper script.
                wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/module.mako").render_unicode(program = self, module = module)
                
                # Write it to our input dir.
                with open(Path(self.destination.calc_dir.input_directory, "{}.{}.in".format(module.name, index +1)), "wt") as wrapper_file:
                    wrapper_file.write(wrapper_body)
                
                # Run Turbomole!
                def run_func(working_dir):
                    subprocess.run(
                        ("bash",),
                        input = wrapper_body,
                        stdout = subprocess.PIPE,
                        stderr = subprocess.STDOUT,
                        universal_newlines = True,
                        cwd = working_dir,
                        check = True
                    )
                    
                # Decide if we need a safe working dir or not.
                if module.unsafe:
                    self.wrap_unsafe_module(run_func)
                    
                else:
                    run_func(self.working_directory)
        
        def parse_results(self):
            """
            Parse the finished calculation result file(s).
            
            Certain flags will be set depending on the status of the calculation. Additionally, a PDF_report object will be stored to self.result
            """
            if "Turbomole-UFF" in self.calculation.CLASS_HANDLE:
                # We can't parse UFF results yet.
                return
            else:
                return super().parse_results()
            
        def get_result(self):
            """
            Get a result set from this calculation.
            """
            # For turbomole, we provide both the main log file and the calculation dir, because these might be in different locations.
            return parse_calculation(self.calc_output_file_path, self.working_directory)    
            