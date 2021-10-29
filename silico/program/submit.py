# General imports.
from pathlib import Path

# Silico imports.
from silico.program.base import Program
from silico.misc.base import to_bool
from silico.interface.urwid.submit.base import Calculation_submitter
from silico.exception.base import Silico_exception
from silico.submit.calculation.base import Calculation_target
from silico.exception.uncatchable import Submission_paused


class Submit_program(Program):
    """
    The Silico submission program.
    """
    
    name = "Silico Calculation Submitter"
    command = "submit"
    description = "submit calculation files"
    help = "Submit calculations"
    usage = """%(prog)s submit ...
       or: %(prog)s result ...
       or: %(prog)s report ..."""
    
    @classmethod
    def arguments(self, sub_parsers_object):
        """
        Add this program's arguments to an argparser subparsers object.
        
        :param sub_parsers_object: An argparser subparsers object to add to (as returned by add_subparsers().
        :returns: The argparser object (which supports add_argument()).
        """
        sub_parser = super().arguments(sub_parsers_object)
        
        sub_parser.add_argument("calculation_files", help = "Calculation input files to submit", nargs = "*", type = Path)
        sub_parser.add_argument("-o", "--output", help = "Base directory to perform calculations in. Defaults to the current directory", default = Path("./"))
        sub_parser.add_argument("-m", "--methods", help = "Methods to perform, identified either by name or by ID", nargs = "*", default = [])
        sub_parser.add_argument("-C", "--charge", help = "Set the molecular charge of all input files. Note that certain calculations will override this value", default = None, type = int)
        sub_parser.add_argument("-M", "--multiplicity", help = "Set the multiplicity of all input files. Note that certain calculations will override this value", default = None, type = int)
        sub_parser.add_argument("--gen3D", help = "Whether to generate 3D coordinates (this will scramble existing atom coordinates). The default is yes, but only if it can be safely determined that the loaded coordinates are not already in 3D)", type = to_bool, default = None)
        
        return sub_parser
    
    def get_methods(self):
        """
        Get a list of methods that we've been asked to submit to.
        """
        return [self.config.methods.resolve_method_string(method_string) for method_string in getattr(self.args, 'methods', [])]
    
    def main(self):
        """
        """
        # Get upset if we have no files.
        if len(self.args.calculation_files) == 0:
            raise Silico_exception("No files to submit")
                        
        # Get upset if we have no methods.
        if len(self.args.methods) == 0:
            raise Silico_exception("No methods to submit")
        
        # Resolve our methods.
        methods = self.get_methods()
        
        try:
            # Arrange our calcs into a linked list.
            first = Calculation_target.link(methods, global_silico_options = self.config)
            
        except Exception:
            raise Silico_exception("Error processing methods to submit")
            
        # Print any warning messages (but only once each).
        warnings = [None]
        for configurables in methods:
            for configurable in configurables:
                    if configurable.warning not in warnings:
                        self.logger.warning(configurable.warning)
        
        # The number of calcs we successfully submitted.
        done = 0
        
        for input_file_path in self.args.calculation_files:
            try:
                # Prepare.
                first.prepare_from_file(self.args.output, input_file_path, molecule_charge = self.args.charge, molecule_multiplicity = self.args.multiplicity, gen3D = self.args.gen3D)
                
                # Go.
                first.submit()
                
                done += 1
                
            except Submission_paused:
                # This is fine.
                done += 1
            
            except Exception:
                # Something went wrong.
                # We don't stop here though, we might have more calcs we can submit.
                self.logger.error("Failed to submit file '{}'".format(input_file_path), exc_info = True)
        
        return done
        
    def load_interface(self, window):
        """
        Function called to get an urwid widget to display for when this subprogram is called interactively.
        """
        # First, resolve any methods we've been given.
        initial_methods = self.get_methods()
        
        return Calculation_submitter(
            window.top,
            self.config.methods,
            initial_files = getattr(self.args, 'calculation_files', []),
            initial_charge = getattr(self.args, 'charge', None),
            initial_mult = getattr(self.args, 'multiplicity', None),
            initial_methods = initial_methods
        )
        
        