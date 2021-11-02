# General imports.
from pathlib import Path
import dill
import logging

# Silico imports.
from silico.program.base import Program
import silico
from silico.exception.uncatchable import Submission_paused
from silico.exception.base import Silico_exception


class Resume_program(Program):
    """
    The Silico method status program.
    """
    
    name = "Silico submission resume"
    command = "resume"
    description = "resume an in-progress silico submission"
    help = "Resume submission (used automatically by part of the submission mechanism)"
    
    
    def __init__(self, destination, *, args, config, logger):
        """
        Constructor for resume programs.
        
        :param destination: A destination_target object from which submission will resume.
        :param args: The command-line arguments the program was started with.
        :param config: A loaded Silico config object.
        :param logger: The logger to use for output.
        """
        super().__init__(args = args, config = config, logger = logger)
        self.destination = destination
        
        # Delete the pickled file to clean up (also prevents us running twice on the same file, which would be bad. Maybe we should do some file locking anyway?)
        try:
            args.resume_file.unlink()
        except Exception:
            logger.error("Failed to delete pickle file", exc_info = True)
        
    @classmethod
    def init_from_argparse(self, args):
        """
        Perform program setup from the data provided by argparse.
        
        :param args: Argparser object.
        :return: A tuple of (args, config, logger).
        """
        # First, sort out our logger.
        logger = logging.getLogger(silico.logger_name)
        
        # Load our pickled class.
        with open(args.resume_file, "rb") as pickle_file:
            destination = dill.load(pickle_file)
            
        # Setup program wide stuff.
        self.program_init(args, destination.program.calculation.silico_options, logger)
        
        return (destination, args, destination.program.calculation.silico_options, logger)
    
    @classmethod
    def load_from_argparse(self, destination, args, config, logger):
        """
        Create a program instance.
        
        :param args: The command-line arguments the program was started with.
        :param config: A loaded Silico config object.
        :param logger: The logger to use for output.
        """
        return self(destination, args = args, config = config, logger = logger)
    
    @classmethod
    def from_argparse(self, args):
        """
        Create a Program object from the data provided by argparse.
        """
        destination, args, config, logger = self.init_from_argparse(args)
        return self._from_argparse(destination, args, config, logger)
    
    @classmethod
    def arguments(self, sub_parsers_object):
        """
        Add this program's arguments to an argparser subparsers object.
        
        :param sub_parsers_object: An argparser subparsers object to add to (as returned by add_subparsers().
        :returns: The argparser object (which supports add_argument()).
        """
        sub_parser = sub_parsers_object.add_parser("resume",
            description = self.description,
            epilog = self.epilog,
            help = self.help
        )
        # Set main function.
        sub_parser.set_defaults(func = self.main)
        
        sub_parser.add_argument("resume_file", help = "Path to the file to resume from (this should be a pickled calculation class)", type = Path)
    
        return sub_parser
    
    def main(self):
        """
        Logic for our program.
        """
        try:
            self.destination.resume()
        except Submission_paused:
            # This is fine.
            pass
        except Exception:
            raise Silico_exception("Error during submission")
        
        
        