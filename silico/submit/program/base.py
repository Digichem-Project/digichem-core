from pathlib import Path
from logging import getLogger
from timeit import default_timer as timer
import datetime
import shutil

from silico import misc
from silico.submit.structure.flag import Flag
from silico.config.configurable.option import Option
from silico.file.convert.babel import Openbabel_converter
from silico.exception.base import Submission_error
from silico.exception.uncatchable import Signal_caught
from silico.file.convert.main import Silico_input
from silico.extract.text import Text_summary_group_extractor
from silico.extract.csv import Long_CSV_group_extractor
from silico.extract.long import Atoms_long_extractor, Orbitals_long_extractor,\
    Beta_long_extractor, SCF_long_extractor, MP_long_extractor, CC_long_extractor,\
    Excited_state_long_extractor, Excited_state_transitions_long_extractor,\
    TDM_long_extractor, Vibrations_long_extractor,\
    Absorption_spectrum_long_extractor, Absorption_energy_spectrum_long_extractor,\
    IR_spectrum_long_extractor, SOC_long_extractor
from silico.submit import Configurable_target
from silico.misc.directory import copytree
import silico.report
from silico.parser import parse_calculation
from silico.report.main.pdf import PDF_report

class Program_target(Configurable_target):
    """
    Top-level class for classes that implement submission to a calculation program.
    """
    
    CLASS_HANDLE = ("program",)
    
    # Configurable options.
    parents = Option("methods", help = "A list of methods that this program is compatible with", required = True, type = list)
    
    ############################
    # Class creation mechanism #
    ############################
    
    class _actual(Configurable_target._actual):
        """
        Inner class for programs.
        """
        
        def __init__(self, method):
            """
            Constructor for program objects.
            
            :param method: A Method_target_actual object that is going to be submitted to.
            """
            # Set our method.
            self.method = method
            self.validate_parent(method)
            # Let our method know who we are.
            self.method.program = self
            # We don't have a calculation yet.
            self.calculation = None
        
            # A Result_set object that will be populated once this calculation has completed.
            self.result = None
            self.start_time = None
            self.end_time = None
            self.duration = None
            
            # An exception caught during the calculation. This will be re-raised once cleanup and analysis has been finished (attempted).
            self.error = None
            
        @property
        def scratch_output(self):
            """
            Path to the scratch folder in which the calculation will be run.
            
            When using scratch and all_output is set to True, this is the folder in which the calculation will be run (the CWD).
            If all_output is False (or no scratch is to be used at all), this property will == None.
            """
            if self.calculation.scratch_directory is None or not self.calculation.scratch_options['all_output']:
                return None
            else:
                return Path(self.calculation.scratch_directory, "Output")
            
        def submit(self):
            """
            Submit this program; running the specified calculation.
            
            Inheriting classes can normally avoid overwriting this method by defining a calculate() method instead.
            
            :param calculation: A Calculation_target_actual object that is going to be submitted.
            """
            # Call method submit first to setup the environment.
            self.method.submit()
            
            # Pre-calc (write input files etc).
            self.pre()
            
            # Run Program (script).
            self.start()            
                
            try:
                
                # Go.
                self.calculate()
                
            except (Signal_caught, KeyboardInterrupt):
                # We've been told to stop (probably by SLURM because we went over time etc).
                self.end(False)
                
                # Continue stopping.
                raise
                
            except Exception as e:
                # Something went wrong.
                self.end(False)
                
                # Raise.
                # Store for later so we can try and generate PDFs and results.
                self.error = Submission_error(self, "Error executing calculation program")
                self.error.__cause__ = e
                
            else:
                # Finished normally.
                self.end(True)
                
            # If we got an error during the calc, re-raise it now.
            if self.error is not None:
                raise self.error
    
        @property
        def success(self):
            """
            Get whether the calculation has completed successfully.
            
            This property has one of 3 values:
                - True, if the calculation finished successfully.
                - False, if the calculation finished normally but the quantum chemistry program signalled an error (eg, failed convergence, unknown keyword etc).
                - None, if the calculation has not yet finished or if the calculation failed catastrophically (and we could not load enough results to determine success or not). 
            """
            if self.result is None:
                return None
            else:
                return self.result.safe_get('metadata', 'success')
                        
            
        def start(self):
            """
            Signals the start of a calculation, this method is called automatically by submit().
            """
            # Set the start and running flags.
            self.method.calc_dir.set_flag(Flag.STARTED)
            self.method.calc_dir.set_flag(Flag.RUNNING)
            
            # Save our start time.
            self.start_timer = timer()
            self.start_date = datetime.datetime.now()
            
            # Log.
            getLogger(silico.logger_name).info("Calculation start on {} ".format(misc.date_to_string(self.start_date)))
        
        def end(self, success = True):
            """
            Signals the end of a calculation, this method is called automatically by submit().
            
            :param success: Should be False if the program finished with error.
            """
            try:        
                # Set our error flag if something went wrong
                if not success:
                    self.method.calc_dir.set_flag(Flag.ERROR)
                
                self.end_timer = timer()
                self.end_date = datetime.datetime.now()
                
                # Assemble our string.
                message = "Calculation end on {}" if success else "Abnormal calculation end on {}"
                message = message.format(misc.date_to_string(self.end_date)) 
                
                if success:
                    getLogger(silico.logger_name).info(message)
                else:
                    getLogger(silico.logger_name).error(message)
                    
                # Work out how much time has passed.
                self.duration = datetime.timedelta(seconds = self.end_timer - self.start_timer)
                getLogger(silico.logger_name).info("Calculation duration: {} ({} total seconds)".format(misc.timedelta_to_string(self.duration), self.duration.total_seconds()))
                
                # Unset our running flag.
                self.method.calc_dir.del_flag(Flag.RUNNING)
                
                ########
                # Post #
                ########
                # We perform post prior to cleanup to save potentially copying large files (.rwf, .chk) from scratch to output.
                # Set Flag.
                self.method.calc_dir.set_flag(Flag.POST)
                
                try:
                    self.post()
                except Exception as e:
                    # Save so we can do cleanup.
                    self.error = e
                
                # Delete Flag.
                self.method.calc_dir.del_flag(Flag.POST)
            
                ###########
                # Cleanup #
                ###########
                # Set our cleanup flag.
                self.method.calc_dir.set_flag(Flag.CLEANUP)
            
                # Call user specified cleanup.
                self.cleanup(success)
                
                # If we were using scratch output, copy back now.
                if self.scratch_output is not None:
                    
                    # Copy.
                    copytree(self.scratch_output, self.method.calc_dir.output_directory)
                    
                    # Delete the scratch output.
                    shutil.rmtree(self.scratch_output)
                
                # If we've been asked to, we'll try and save what remains of the scratch directory (might contain something useful).
                if self.calculation.scratch_options['use_scratch'] and ((self.calculation.scratch_options['rescue'] and not success) or self.calculation.scratch_options['keep']):
                    # Check to see if there's anything in scratch.
                    scratch_content = -1
                    try:
                        scratch_content = len(list(self.calculation.scratch_directory.iterdir()))
                    finally:
                        # Move our scratch dir (or at least try to) if it's not empty (or we couldn't determine if it's empty).
                        if scratch_content != 0:
                            shutil.move(str(self.calculation.scratch_directory), str(self.method.calc_dir.scratch_directory))
                    
            finally:
                # Remove the scratch directory, forcibly if need be.
                if self.calculation.scratch_options['use_scratch'] and (success or self.calculation.scratch_options['force_delete']):
                    try:
                        shutil.rmtree(self.calculation.scratch_directory)
                    except FileNotFoundError:
                        # This is ok, means the scratch has already been deleted (or moved).
                        pass
                    except Exception:
                        # Some other error occurred that prevented us from deleting the scratch.
                        # This is odd, annoying and possibly evidence of a bug, but shouldn't stop us from continuing.
                        getLogger(silico.logger_name).warning("Failed to delete scratch directory '{}'".format(self.calculation.scratch_directory), exc_info = True)
                    
                # Done cleanup.
                self.method.calc_dir.del_flag(Flag.CLEANUP)
        
        @property
        def calc_output_file_path(self):
            """
            Path to the main calculation output file.
            """
            raise NotImplementedError()
        
        @property
        def working_directory(self):
            """
            The working directory in which the calculation will be performed. 
            """
            # Decide on where we are running.
            if self.scratch_output is not None:
                return self.scratch_output
            else:
                return self.method.calc_dir.output_directory
        
        def pre(self):
            """
            Perform pre-setup for a calculation.
            """
            # Make our scratch directory if we're using scratch.
            if self.calculation.scratch_directory is not None:
                try:
                    # Make the folder.
                    self.calculation.scratch_directory.mkdir(parents = True)
                except FileExistsError:
                    # The scratch folder already existing is actually pretty serious; it's supposed to be unique to us, so if it already exists something's probably gone wrong.
                    # However, some SLURM implementations automatically create our scratch dir for us; because this scenario is common (?), we currently print a warning and continue.
                    # This may change in the future.
                    getLogger(silico.logger_name).warning("Could not create scratch directory '{}' because it already exists; continuing anyway".format(self.calculation.scratch_directory)) 
                except Exception:
                    raise Submission_error(self, "unable to create scratch directory")
            
            # Also make our scratch output directory if we're using one.
            if self.scratch_output is not None:
                try:        
                    self.scratch_output.mkdir()
                except Exception as e:
                    raise Submission_error(self, "Failed to make scratch output subdirectory") from e
            
            # Write our input file in .si format for easy reuse.
            if hasattr(self.calculation, "input_coords"):
                with open(Path(self.method.calc_dir.input_directory, self.calculation.molecule_name).with_suffix(".si"), "wt") as input_file:
                    self.calculation.input_coords.to_file(input_file)
    
        def calculate(self):
            """
            This method should be implemented in child classes to perform the designated calculation.
            """
            raise NotImplementedError()
        
        def cleanup(self, success):
            """
            This method should be implemented in child classes to perform cleanup after the calculation.
            
            :param success: True if the calculation finished normally.
            """
            pass
        
        def get_result(self):
            """
            Get a result set from this calculation.
            """
            return parse_calculation(self.calc_output_file_path)
        
        def get_report(self):
            """
            Get a report suitable for parsing this type of calculation.
            """
            return PDF_report(self.results, options = self.calculation.silico_options)
            
            
        def parse_results(self):
            """
            Parse the finished calculation result file(s).
            
            Certain flags will be set depending on the status of the calculation. Additionally, a PDF_report object will be stored to self.result
            """
            # Try and load our results.
            # We need to know whether the calculation was successful or not, so we make no effort to catch exceptions here.
            try:
                try:
                    self.result = self.get_result()
                except Exception as e:
                    raise Submission_error(self, "failed to process completed calculation results") from e
                
                # See if our calculation was successful or not.
                if self.result.metadata.success and self.error is None:
                    self.method.calc_dir.set_flag(Flag.SUCCESS)
                else:
                    # No good.
                    #self.method.calc_dir.set_flag(Flag.ERROR)
                    raise Submission_error(self, "an error occurred during the calculation; check calculation output for what went wrong")
                    
                # Also check optimisation convergence.
                if self.result.metadata.optimisation_converged is not None and "Optimisation" in self.result.metadata.calculations:
                    if self.result.metadata.optimisation_converged:
                        # Converged.
                        self.method.calc_dir.set_flag(Flag.CONVERGED)
                    else:
                        self.method.calc_dir.set_flag(Flag.NOT_CONVERGED)
                        raise Submission_error(self, "the optimisation did not converge")
                        
            except Exception:
                # No good.
                self.method.calc_dir.set_flag(Flag.ERROR)
                raise
        
        def post(self):
            """
            Perform post analysis and cleanup, this method is called after a calculation has finished.
            """                
            
            # First, make our result directory.
            try:
                self.method.calc_dir.result_directory.mkdir()
            except FileExistsError:
                pass
            
            # Next, load calc results.
            try:
                self.parse_results()
            except Exception as e:
                self.error = e
            
            # If we've been asked to write result files, do so.
            try:
                if self.calculation.write_summary:
                    self.write_summary_files()
            except Exception:
                getLogger(silico.logger_name).warning("Failed to write calculation result summary files", exc_info = True)
                
            # Write XYZ file.
            try:
                self.write_xyz_file()
            except Exception:
                getLogger(silico.logger_name).warning("Failed to write XYZ result file", exc_info = True)
                
            # Write .si file.
            try:
                self.write_si_file()
            except Exception:
                getLogger(silico.logger_name).warning("Failed to write silico (.si) result file", exc_info = True)
                
            # Similarly, if we've been asked to write a report, do that.
            # TEMP: Don't write reports if we fail, see issue #29.
            if self.error is None:
                try:
                    if self.calculation.write_report:
                        self.write_report_files()
                except Exception:
                    getLogger(silico.logger_name).warning("Failed to write calculation report", exc_info = True)
            else:
                getLogger(silico.logger_name).info("Skipping report generation because calculation did not finish successfully")
                
                

        
        def write_summary_files(self):
            """
            Write text result files (like with cresult) from this calculation.
            """
            # First, write a text summary.
            Text_summary_group_extractor(ignore = True, config = self.calculation.silico_options).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name +".summary"))
            
            # Atoms.
            Long_CSV_group_extractor(Atoms_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".atoms.csv"))
            
            # Alpha. We'll use a different file name depending on whether we are restricted or unrestricted.
            Long_CSV_group_extractor(Orbitals_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".{}.csv".format("orbitals" if self.result.metadata.orbital_spin_type == "restricted" else "alpha")))
            # And beta (which can only be for unrestricted.
            Long_CSV_group_extractor(Beta_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".beta.csv"))
            
            # Energies.
            Long_CSV_group_extractor(SCF_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".SCF.csv"))
            Long_CSV_group_extractor(MP_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".MP.csv"))
            Long_CSV_group_extractor(CC_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".CC.csv"))
            
            # Excited states.
            Long_CSV_group_extractor(Excited_state_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".ES.csv"))
            Long_CSV_group_extractor(Excited_state_transitions_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".transitions.csv"))
            Long_CSV_group_extractor(TDM_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".TDM.csv"))
            Long_CSV_group_extractor(Absorption_spectrum_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".UV-Vis.csv"))
            Long_CSV_group_extractor(Absorption_energy_spectrum_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".absorptions.csv"))
            Long_CSV_group_extractor(SOC_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".SOC.csv"))
            
            # And vibrations.
            Long_CSV_group_extractor(Vibrations_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".vibrations.csv"))
            Long_CSV_group_extractor(IR_spectrum_long_extractor(ignore = True, config = self.calculation.silico_options)).write_single(self.result, Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".IR.csv"))
            
        def write_report_files(self):
            """
            Write report files for this calculation.
            """
            # Load report.
            report = self.get_report()
            # The full report.
            report.write(Path(self.method.calc_dir.report_directory, self.calculation.molecule_name + ".pdf"))
            # And atoms.
            report.write(Path(self.method.calc_dir.report_directory, self.calculation.molecule_name + ".atoms.pdf"), report_type = "atoms")
            
        def write_xyz_file(self):
            """
            Write an XYZ file from the finished calculation results.
            """
            # Get our converter.
            conv = Openbabel_converter.from_file(self.next_coords, input_file_type = self.calculation.OUTPUT_COORD_TYPE)
            
            # Open output file.
            with open(Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".xyz"), "wt") as xyz_file:
                # Write.
                xyz_file.write(conv.convert("xyz"))
                
        def write_si_file(self):
            """
            Write a .si file from the finished calculation results.
            """
            # Get our convertor.
            conv = Silico_input.from_file(self.next_coords, file_type = self.calculation.OUTPUT_COORD_TYPE, name = self.calculation.input_coords.name, charge = self.calculation.input_coords.charge, multiplicity = self.calculation.input_coords.multiplicity)
            
            # Write new file.
            with open(Path(self.method.calc_dir.result_directory, self.calculation.molecule_name + ".si"), "wt") as si_file:
                # Write.
                conv.to_file(si_file)
            
            
        