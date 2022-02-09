# General imports.
from pathlib import Path
import subprocess
from mako.lookup import TemplateLookup

# Silico imports.
import silico.logging
from silico.submit.program import Program_target
from silico.file.fchk import Chk_to_fchk
from silico.config.configurable.option import Option
from silico.parser.base import parse_calculation

class Gaussian(Program_target):
    """
    Top level class for submitting calculations to Gaussian.
    """
    
    CLASS_HANDLE = ("Gaussian",)
    
    # Configurable options.
    executable = Option(help = "Name/path of the main Gaussian executable", required = True, type = str)
    root_environ_name = Option(help = "The name of the environmental variable Gaussian looks for to find 'root'", required = True, type = str)
    root = Option(help = "Path to the directory where gaussian is installed", required = True, type = Path)
    init_file = Option(help = "Path to the gaussian .profile script which is run to set-up gaussian", required = True, type = str)
        
        
    ############################
    # Class creation mechanism #
    ############################
    
    class _actual(Program_target._actual):
        """
        Inner class for programs.
        """
        
        def __init__(self, *args, **kwargs):
            """
            Constructor for Gaussian programs.
            """
            super().__init__(*args, **kwargs)
            
        @property
        def com_file_path(self):
            """
            Path to the (ready-to-go) input file. Note that although this is known as the com file, it may infact have any extension (.com and .gjf are most common).
            """
            return Path(self.destination.calc_dir.input_directory, self.calculation.com_file_name)
    
        @property
        def log_file_path(self):
            """
            Default path to the .log output file written to by Gaussian, see log_file_path for where the log file is currently.
            """
            return Path(self.destination.calc_dir.output_directory, self.calculation.com_file_name).with_suffix(".log")
    
        @property
        def calc_output_file_path(self):
            """
            Path to the main calculation output file.
            
            For Gaussian calcs this is the .log file.
            """
            return self.log_file_path
                
        @property
        def next_coords(self):
            """
            Path to the output coordinate file that should be used for any subsequent calculations.
            """
            return self.log_file_path
            
        @property
        def fchk_file_path(self):
            """
            Path to the formatted checkpoint .fchk file.
            """
            return Path(self.destination.calc_dir.output_directory, self.calculation.chk_file_name).with_suffix(".fchk")
        
        @property
        def formchk_executable(self):
            """
            Path to the formchk executable of this Gaussian installation.
            """
            return Path(self.root, 'formchk')
    
        def pre(self):
            """
            Pre-calculation setup for Gaussian.
            """
            # Call parent for setup first.
            super().pre()
            
            # Write our input file to our calculation Input directory.
            with open(self.com_file_path, "wt") as com_file:
                com_file.write(self.calculation.com_file_body)
                
            # Set some file locations based on whether we're using scratch output or not.
            if self.scratch_output is not None:
                # We're using scratch output, the .chk and .rwf files will be in the scratch dir.
                self.chk_file_path = Path(self.scratch_output, self.calculation.chk_file_name)
                self.rwf_file_path = Path(self.scratch_output, self.calculation.rwf_file_name)
            else:
                # Not using scratch output, .rwf and .chk will be in normal output dir.
                self.chk_file_path = Path(self.destination.calc_dir.output_directory, self.calculation.chk_file_name)
                self.rwf_file_path = Path(self.destination.calc_dir.output_directory, self.calculation.rwf_file_name)
    
        def calculate(self):
            """
            Main submission method; the calculation will be run here (for which this method will block, possibly for hours+).
            """                            
            # First, get our wrapper script.
            gaussian_wrapper_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/gaussian/wrapper.mako").render_unicode(program = self)
            
            # Run Gaussian!
            subprocess.run(
                ['bash'],
                input = gaussian_wrapper_body,
                universal_newlines = True,
                cwd = self.working_directory,
                check = True,
                # Capture output.
                stdout = subprocess.PIPE,
                stderr = subprocess.STDOUT
                )
        
        def post(self):
            """
            Post submission method.
            """
            # Chk/fchk management. Do this before making report (in super()) to avoid making fchk twice.
            try:
                # Create an fchk file if asked.
                if self.calculation.convert_chk:
                    fchk_file = Chk_to_fchk(self.fchk_file_path, chk_file = self.chk_file_path, memory = self.calculation.memory, formchk_executable = self.formchk_executable)
                    fchk_file.get_file()
            except Exception:
                silico.logging.get_logger().error("Failed to create fchk file", exc_info = True)
            else:
                try:
                    # Now delete the chk file if we were asked to.
                    if not self.calculation.keep_chk:
                        self.chk_file_path.unlink()
                except FileNotFoundError:
                    # We can ignore not finding the file; we were trying to get rid of it anyway.
                    pass
                except Exception:
                    silico.logging.get_logger().error("Failed to delete chk file", exc_info = True)
            
            # Use our parent to create result and report files.
            super().post()
            
            # Remove rwf if we've been asked to.
            try:
                if not self.calculation.keep_rwf:
                    self.rwf_file_path.unlink()
            except FileNotFoundError:
                # We can ignore not finding the file; we were trying to get rid of it anyway.
                pass
            except Exception:
                silico.logging.get_logger().error("Failed to delete rwf file", exc_info = True)
                
        def cleanup(self, success):
            """
            """
            super().cleanup(success)
            
            # Update file locations.
            # Regardless of whether we were using scratch output or not, all files should now be in the normal output dir (or have been deleted).
            self.chk_file_path = Path(self.destination.calc_dir.output_directory, self.calculation.chk_file_name)
            self.rwf_file_path = Path(self.destination.calc_dir.output_directory, self.calculation.rwf_file_name)
            
        def get_result(self):
            """
            Get a result set from this calculation.
            """
            # We provide paths to aux files as well as the main output.
            return parse_calculation(self.calc_output_file_path, chk_file = self.chk_file_path, fchk_file = self.fchk_file_path, rwf_file = self.rwf_file_path)
        