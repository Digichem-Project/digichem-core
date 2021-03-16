# General imports.
from pathlib import Path
import shutil
from subprocess import TimeoutExpired, CalledProcessError
import os
import re
from mako.lookup import TemplateLookup
import subprocess

# Silico imports.
from silico.submit.program.base import Program_target
from silico.config.configurable.option import Option
from silico.exception.base import Submission_error
from silico.misc.directory import copytree
import silico.report


class Turbomole(Program_target):
    """
    Top level class for submitting calculations to Turbomole.
    """
    
    CLASS_HANDLE = ("Turbomole",) 

    define_executable = Option(help = "Name/path of the define executable", default = "define", type = str)
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
            return Path(self.method.calc_dir.input_directory, "coord")
        
        @property
        def control_file_path(self):
            """
            Path to the input control file.
            """
            return Path(self.method.calc_dir.input_directory, "control")
        
        @property
        def define_output_path(self):
            """
            Path to the file in which define output is written.
            """
            return Path(self.method.calc_dir.output_directory, "define.out")
        
        @property
        def turbomole_output_path(self):
            """
            Path to the file in which turbomole output is written.
            """
            return Path(self.method.calc_dir.output_directory, self.calculation.molecule_name).with_suffix(".log")
        
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
            return Path(self.method.calc_dir.output_directory, "coord")
        
        @property
        def scratch_base(self):
            """
            The basename of the directories to use as scratch.
            When running in multi-process mode (MPI), the scratch directory given to turbomole is only used as a suggestion; each node will create its own folder by APPENDING "-001" etc to the given file name.
            
            This is annoying, so when in MPI, we make sure each of these folders is created inside the specified scratch dir.
            """
            if self.calculation.scratch_directory is None:
                return None
            elif self.calculation.parallel_mode == "MPI":
                return Path(self.calculation.scratch_directory, "node")
            else:
                return self.calculation.scratch_directory
                    
        @property
        def define_input_path(self):
            """
            Path to the input file used to power the define input generator.
            """
            return Path(self.method.calc_dir.input_directory, "define.in")
        
        def define(self):
            """
            Run setup for turbomole.
            
            Normally this involves running define, but some turbomole calcs have a different setup (eg, UFF).
            """
            # Write define input file.
            with open(self.define_input_path, "wt") as define_input:
                define_input.write(self.calculation.define_input)
            
            # Get our wrapper script.
            wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/define_wrapper.mako").render_unicode(program = self)
                        
            # Run control to generate input.
            try:
                subprocess.run(
                    ("bash",),
                    input = wrapper_body,
                    stdout = subprocess.PIPE,
                    stderr = subprocess.STDOUT,
                    universal_newlines = True,
                    timeout = self.calculation.define_timeout,
                    cwd = self.method.calc_dir.prep_directory,
                    check = True
                )
                
            except TimeoutExpired as e:
                # Ran out of time, probably got stuck.
                raise Submission_error(self, "Program 'define' failed to finish executing in {} s, check output file '{}' for errors".format(self.calculation.define_timeout, self.define_output_path)) from e    
            
            except CalledProcessError as e:
                # Something went wrong.
                e.__context__ = None
                raise Submission_error(self, "define did not exit successfully:\n{}".format(e.stdout)) from e
            
            except Exception as e:
                # Something else.
                raise e from None
            
            # Sadly, some options are not supported by define and have to be appended manually.
            if self.calculation.optimisation_state is not None:
                self.add_control_option("$exopt {}".format(self.calculation.optimisation_state))
            
        def add_control_option(self, option):
            """
            Manually write to the control file.
            
            The control file should already exist before this method is called.
            """
            control_name = Path(self.method.calc_dir.prep_directory, "control")
            tmp_control_name = Path(self.method.calc_dir.prep_directory, "control.tmp")
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
            # Get our wrapper script.
            wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/turbomole_wrapper.mako").render_unicode(program = self)
            
            # Copy atoms to prep dir.
            shutil.copy2(self.coord_file_path, self.method.calc_dir.prep_directory)
            
            # Run control to generate input.
            # We check for errors, but sadly uff doesn't seem to raise any?
            try:
                subprocess.run(
                    ("bash",),
                    input = wrapper_body,
                    stdout = subprocess.PIPE,
                    stderr = subprocess.STDOUT,
                    universal_newlines = True,
                    cwd = self.method.calc_dir.prep_directory,
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
            with open(Path(self.method.calc_dir.prep_directory, "control"), "r+") as control_file:
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
            
            # If we have a directory calc, copy the old directory to our new Prep directory.
            if self.calculation.DIRECTORY_CALCULATION:
                copytree(self.calculation.input, self.method.calc_dir.prep_directory)
                
                # Clean up some old files which may be left by the previous calc if run with silico.
                for delete_file in (self.define_output_path.name, self.method.calc_dir.log_file.name):
                    delete_path = Path(self.method.calc_dir.prep_directory, delete_file)
                    
                    try:
                        delete_path.unlink()
                    except FileNotFoundError:
                        # This is ok.
                        pass
                    
            else:
                # Our calc is not a directory calc, write our input file to our calculation Input directory.
                with open(self.coord_file_path, "wt") as coord_file:
                    coord_file.write(self.calculation.input_coords.to_format("tmol"))
                    
                # Create our prep dir.
                self.method.calc_dir.prep_directory.mkdir()
                            
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
                    copytree(self.method.calc_dir.prep_directory, self.scratch_output)
                except Exception as e:
                    raise Submission_error(self, "Failed to copy input to scratch subdirectory") from e
                
            else:
                # Copy input to normal output folder.
                copytree(self.method.calc_dir.prep_directory, self.method.calc_dir.output_directory)
                
            # Copy prep dir to input dir and delete.
            copytree(self.method.calc_dir.prep_directory, self.method.calc_dir.input_directory)
            shutil.rmtree(self.method.calc_dir.prep_directory)
            
        def calculate(self):
            """
            Main submission method; the calculation will be run here (for which this method will block, possibly for hours+).
            """                                        
            # Get our wrapper script.
            wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/turbomole_wrapper.mako").render_unicode(program = self)            
            
            # Run Turbomole!
            subprocess.run(
                ("bash",),
                input = wrapper_body,
                stdout = subprocess.PIPE,
                stderr = subprocess.STDOUT,
                universal_newlines = True,
                cwd = self.working_directory,
                check = True
            )
        
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
                
        def get_report(self):
            """
            Get a report suitable for parsing this type of calculation.
            """            
            return silico.report.from_result(
                self.result,
                options = self.calculation.silico_options,
                turbomole_calculation = self.calculation
            )
            