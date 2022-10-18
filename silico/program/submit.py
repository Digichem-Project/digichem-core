# General imports.
from pathlib import Path
import argparse

# Silico imports.
from silico.program.base import Program
from silico.misc.base import to_bool
from silico.interface.urwid.submit import Submit_interface
from silico.exception.base import Silico_exception
from silico.submit.calculation.base import Calculation_target
from silico.exception.uncatchable import Submission_paused
from silico.input import si_from_file
from silico.submit.base import parse_method_from_file


class Method_action(argparse.Action):
    """
    """
    
    def __init__(self, *args, method_type, **kwargs):
        self.method_type = method_type
        super().__init__(*args, **kwargs)
    
    def __call__(self, parser, namespace, values, option_string = None):
        """
        """
        try:
            methods = getattr(namespace, self.dest)
            
        except AttributeError:
            methods = []
            setattr(namespace, self.dest, methods)
            
        for method in values:
            methods.append((self.method_type, method))
            
        

class Submit_program(Program):
    """
    The Silico submission program.
    """
    
    name = "Silico Calculation Submitter"
    command = "submit"
    aliases = ["s", "sub"]
    description = "Submit calculation files"
    help = "Submit calculations"
    usage = """%(prog)s [options] coord1.file [coord2.file] -c d1/p1/c1 [d2/p2/c2]
       %(prog)s [options] coord1.file [coord2.file] -m method1.file [method2.file]
       """
    
    @classmethod
    def arguments(self, sub_parsers_object):
        """
        Add this program's arguments to an argparser subparsers object.
        
        :param sub_parsers_object: An argparser subparsers object to add to (as returned by add_subparsers().
        :returns: The argparser object (which supports add_argument()).
        """
        sub_parser = super().arguments(sub_parsers_object)
        
        sub_parser.add_argument("coordinate_files", help = "Coordinate input files to submit", nargs = "*", type = Path)
        sub_parser.add_argument("-O", "--output", help = "Base directory to perform calculations in. Defaults to the current directory", default = Path("./"))
        sub_parser.add_argument("-m", "--method-files", help = "Methods to perform, identified by file name", nargs = "*", default = [], dest = "methods", action = Method_action, method_type = "file")
        sub_parser.add_argument("-c", "--method-codes", help = "Methods to perform, identified either by name or by ID (such, as 1/1/1)", nargs = "*", default = [], dest = "methods", action = Method_action, method_type = "code")
        sub_parser.add_argument("-C", "--charge", help = "Set the molecular charge of all input files. Note that certain calculations will override this value", default = None, type = int)
        sub_parser.add_argument("-M", "--multiplicity", "--mult", help = "Set the multiplicity of all input files. Note that certain calculations will override this value", default = None, type = int)
        sub_parser.add_argument("--gen3D", help = "Whether to generate 3D coordinates from any given 2D input coordinate files (which will scramble existing atom coordinates). The default is yes, but only if it can be safely determined that the loaded coordinates are not already in 3D)", type = to_bool, default = True)
        sub_parser.add_argument("-p", "--prepare-only", help = "Whether to only perform setup for the calculation without actually performing it", action = "store_true")
        sub_parser.add_argument("--exit-status", help = "If given and the submission is successful, the exit status set by silico will be a positive integer equal to the number of files successfully submitted. If not given (the default), the exit status will be zero if all given files were submitted successfully, and a negative integer otherwise.", action = "store_true")
        
        return sub_parser
    
    def __init__(self, coords, methods, args, config, logger):
        """
        Constructor for submit programs.
        
        :param output: A pathlib Path to a directory to write output to. Each submitted molecule will have a subdirectory under the output directory.
        :param coords: A list of Silico_coords objects representing coordinates to submit.
        :param methods: A list of method tuples to submit to.
        """
        super().__init__(args = args, config = config, logger = logger)
        self.coords = coords
        self.methods = methods
        
    @classmethod
    def load_from_argparse(self, args, config, logger):
        """
        Create a Program object from the data provided by argparse.
        
        :param args: The command-line arguments the program was started with.
        :param config: A loaded Silico config object.
        :param logger: The logger to use for output.
        """
        # First get our list of methods.
        methods = []
        for method_type, method_id in args.methods:
            if method_type == "code":
                methods.append(config.methods.resolve_method_string(method_id))
                
            elif method_type == "file":
                methods.append(parse_method_from_file(method_id, config))
        
        # Load coordinates.
        coords = [si_from_file(file, gen3D = args.gen3D, charge = args.charge, multiplicity = args.multiplicity) for file in args.coordinate_files]
                
        return self(coords, methods, args = args, config = config, logger = logger)
    
    def main(self):
        """
        """
        # Get upset if we have no files.
        if len(self.coords) == 0:
            raise Silico_exception("No files specified")
                        
        # Get upset if we have no methods.
        if len(self.methods) == 0:
            raise Silico_exception("No methods specified")
        
        try:
            # Arrange our calcs into a linked list.
            first = Calculation_target.link(self.methods, global_silico_options = self.config, prepare_only = self.args.prepare_only)
            
        except Exception:
            raise Silico_exception("Error processing methods to submit")
            
        # Print any warning messages (but only once each).
        warnings = [None]
        for configurables in self.methods:
            for configurable in configurables:
                    if configurable.meta['warning'] not in warnings:
                        warnings.append(configurable.meta['warning'])
                        self.logger.warning(configurable.meta['warning'])
        
        # The number of calcs we successfully submitted.
        done = 0
        
        for coord in self.coords:
            try:
                self.logger.info("Submitting file '{}'...".format(coord.file_name))
                
                # Prepare.
                first.prepare(self.args.output, coord)
                
                # Go.
                first.submit()
                
                done += 1
                
            except Submission_paused:
                # This is fine.
                done += 1
            
            except Exception:
                # Something went wrong.
                # We don't stop here though, we might have more calcs we can submit.
                self.logger.error("Failed to submit file '{}'".format(coord.file_name), exc_info = True)
        
        self.logger.info("Successfully submitted {} file(s)".format(done))
        
        # What we return depends on the --exit-status argument.
        if not self.args.exit_status:
            # This the orthodox return type.
            # Return zero if all our calcs submitted successfully, -1 otherwise.
            return 0 if done == len(self.coords) else -1
        
        else:
            # This is the unorthodox return type.
            # Return the number of calcs submitted successfully.
            return done
        
    def load_interface(self, window):
        """
        Function called to get an urwid widget to display for when this subprogram is called interactively.
        """
        return Submit_interface(window, self)
        
        