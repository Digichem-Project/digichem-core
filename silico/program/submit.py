# General imports.
from pathlib import Path

# Silico imports.
from silico.program.base import Program
from silico.misc.base import to_bool
from silico.interface.urwid.submit.base import Calculation_submitter
from silico.exception.base import Silico_exception
from silico.submit.calculation.base import Calculation_target
from silico.exception.uncatchable import Submission_paused
from silico.file.convert.main import Silico_input


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
        sub_parser.add_argument("--gen3D", help = "Whether to generate 3D coordinates (this will scramble existing atom coordinates). The default is yes, but only if it can be safely determined that the loaded coordinates are not already in 3D)", type = to_bool, default = True)
        
        return sub_parser
    
    def __init__(self, coords, methods, args, config, logger):
        """
        Constructor for submit programs.
        
        :param output: A pathlib Path to a directory to write output to. Each submitted molecule will have a subdirectory under the output directory.
        :param coords: A list of Silico_input objects representing coordinates to submit.
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
        methods = [config.methods.resolve_method_string(method_string) for method_string in getattr(args, 'methods', [])]
        
        # Load coordinates.
        coords = [Silico_input.from_file(file, gen3D = args.gen3D, charge = args.charge, multiplicity = args.multiplicity) for file in args.calculation_files]
                
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
            first = Calculation_target.link(self.methods, global_silico_options = self.config)
            
        except Exception:
            raise Silico_exception("Error processing methods to submit")
            
        # Print any warning messages (but only once each).
        warnings = [None]
        for configurables in self.methods:
            for configurable in configurables:
                    if configurable.warning not in warnings:
                        self.logger.warning(configurable.warning)
        
        # The number of calcs we successfully submitted.
        done = 0
        
        for coord in self.coords:
            try:
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
        
        self.logger.info("Successfully submitted {} calculations".format(done))
        
        return done
        
    def load_interface(self, window):
        """
        Function called to get an urwid widget to display for when this subprogram is called interactively.
        """
        return Calculation_submitter(window, self)
        
        