#!/usr/bin/env python3

"""Main entry point to Silico."""

# This should suppress a matplotlib warning when we compile with pyinstaller.
import warnings
warnings.filterwarnings("ignore", "(?s).*MATPLOTLIBDATA.*", category = UserWarning)
# This one is caused by some minor bug when plotting graphs.
warnings.filterwarnings("ignore", "(?s).*Source ID .* was not found when attempting to remove it.*")
# This one is caused by the way tight_layout works. It's fine to ignore here, but not so good if we want to plot interactively.
warnings.filterwarnings("ignore", "(?s).*tight_layout: falling back to Agg renderer*", category = UserWarning)

# This warning is when we try to create a graph with no height, normally caused by an energy convergence graph with no change in energy.
warnings.filterwarnings("ignore", "(?s).*Attempting to set identical bottom == top*", category = UserWarning)
# This warning appear to arise from matplotlib internally, may become a problem in the future...
warnings.filterwarnings("ignore", '(?s).*got unexpected keyword argument "dpi"*')


import sys
# These two suppress weasyprint warnings about incompatible libraries (which we ignore when freezing).
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    warnings.filterwarnings("ignore", "(?s).*@font-face support needs Pango >= 1.38*", category = UserWarning)
    warnings.filterwarnings("ignore", "(?s).*There are known rendering problems and missing features with cairo < 1.15.4*", category = UserWarning)

# Increase verbosity of logging to catch errors during setup.
# This will be altered once config data has been loaded.
import silico.log
silico.log.set_logging_level("DEBUG")

# General imports.
import pydevd;pydevd.settrace()

# Silico imports.
from silico.program.config.main import Config_program
from silico.program.convert import Convert_program
from silico.program.interactive import Interactive_program
from silico.program.report import Report_program
from silico.program.result import Result_program
from silico.program.resume import Resume_program
from silico.program.status import Status_program
from silico.program.submit import Submit_program
from silico.program.base import Program


# A list of each of the sub programs we know about.
PROGRAMS = [
    Config_program,
    Convert_program,
    Interactive_program,
    Report_program,
    Result_program,
    Resume_program,
    Status_program,
    Submit_program]

def get_argparser():
    """
    A helper function which returns an argaprser object for all our subprograms.
    
    This function is separate from main to help with autodocing.
    """
    # Configure our argument parser.
    # This is the top-level parser; each sub program will specify its own sub-parser.
    parser = Program.top_level_arguments()
    subparser = parser.add_subparsers(dest="prog")
    subparser.required = True
        
    # Create sub parsers for each subprogram. Each will define its own parser.
    for sub_program in PROGRAMS:
        sub_program.arguments(subparser)
        
    return parser

def main():
    """
    Main entry point for the program.
    """
    # ----- Program init -----
    parser = get_argparser()
        
#     # A hack to pick a default program if none has been chosen.
#     if len(sys.argv) < 2 or ( sys.argv[1] not in ["--help", "-h", "-v"] and sys.argv[1] not in subparser.choices ):
#         sys.argv.insert(1, "submit")
#         sys.argv.insert(2, "-I")
    
    # Process command line arguments.
    args = parser.parse_args()
    
    # Get our logger.
    logger = silico.log.get_logger()
    
    try:
        # First, get our chosen subprogram.
        inner_program = args.prog_cls.from_argparse(args)
        outer_program = inner_program
        
        #inner_program.program_init()
        
        inner_program.logger.debug("Startup completed")
    
    except Exception:
        # Couldn't start.
        logger.error("failed to start program", exc_info = True)
        return -2
    
    # What we do next depends on whether we're running interactively or not.
    if getattr(args, 'interactive', False) and type(outer_program) != Interactive_program:
        # Because we're running interactively, we'll actually use the interactive program rather than the chosen program.
        outer_program = Interactive_program(inner_program.args, inner_program.config, inner_program.logger, initial = inner_program)

    return outer_program.main_wrapper()


    
# If we've been invoked as a program, call main().    
if __name__ == '__main__':
    sys.exit(main())
    