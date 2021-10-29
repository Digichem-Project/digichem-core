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
    
    
    def __init__(self, args, logger):
        """
        Unlike normal programs, we do not load our local config option.
        
        Instead, we use the config loaded by our original calculation.
        """
        # Load the pickled class.
        with open(args.resume_file, "rb") as pickle_file:
            self.destination = dill.load(pickle_file)
    
        # Delete the pickled file to clean up (also prevents us running twice on the same file, which would be bad. Maybe we should do some file locking anyway?)
        try:
            args.resume_file.unlink()
        except Exception:
            logger.error("Failed to delete pickle file", exc_info = True)
            
        super().__init__(self, args, self.destination.program.calculation.silico_options, logger)
        
    @classmethod
    def from_argparse(self, args, logger_name = None):
        """
        Create an instance of this program from the data provided by argparse.
        
        :param args: Argparser object.
        :param logger_name: The name of a logger to use for this program's output.
        :returns: An instance of this program.
        """
        logger_name = logger_name if logger_name is not None else silico.logger_name
        # First, sort out our logger.
        logger = logging.getLogger(logger_name)
        
        return self(args, logger)
    
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
        
        
        