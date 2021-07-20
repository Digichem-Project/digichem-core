#!/usr/bin/env python3.8

# Main entry point to Silico.

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
# These    two suppress weasyprint warnings about incompatible libraries (which we ignore when freezing).
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    warnings.filterwarnings("ignore", "(?s).*@font-face support needs Pango >= 1.38*", category = UserWarning)
    warnings.filterwarnings("ignore", "(?s).*There are known rendering problems and missing features with cairo < 1.15.4*", category = UserWarning)

# Init openbabel.
import silico
silico.init_obabel()

# General imports.
import argparse

# Silico imports.
import silico.program.submit
import silico.program.result
import silico.program.report
import silico.program.resume
import silico.program.status
import silico.program.convert
import silico.program.config.main

# Make configurables available
from silico.submit.method import *
from silico.submit.program import *
from silico.submit.calculation import *
from silico.submit.basis import *

# Printable name of this program.
NAME = "Silico"
DESCRIPTION = "computational chemistry management"
EPILOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))
USAGE = """%(prog)s submit ...
   or: %(prog)s result ...
   or: %(prog)s report ..."""

def main():
    """
    Main entry point for the program.
    """
    # ----- Program init -----
    
    # Configure our argument parser.
    # This is the top-level parser; each sub program will specify its own sub-parser.
    parser = argparse.ArgumentParser(
        description = DESCRIPTION,
        epilog = EPILOG)
    parser.add_argument("-v", "--version", action = "version", version = str(silico.version))
    subparser = parser.add_subparsers(dest="prog")
    
    # Create sub parsers for each sub-program. Each will define its own parser.
    silico.program.submit.arguments(subparser)
    silico.program.convert.arguments(subparser)
    silico.program.result.arguments(subparser)
    silico.program.report.arguments(subparser)
    #silico.program.spreadsheet.arguments(subparser)
    silico.program.status.arguments(subparser)
    silico.program.config.main.arguments(subparser)
    silico.program.resume.arguments(subparser)
    
#     # Before we call parse_args(), we check to see if a sub command has been chosen.
#     # If not, we'll add 'submit' as the default.
#     # This is a bit hacky, but works ok.
#     if len(sys.argv) <= 1 or (
#         sys.argv[1] not in ["submit", "result", "report", "resume", "status", "convert"]
#         and "-h" not in sys.argv
#         and "--help" not in sys.argv
#         and "-v" not in sys.argv
#         and "--version" not in sys.argv
#     ):
#         sys.argv.insert(1, "submit")
    
    # Process command line arguments.
    args = parser.parse_args()
        
    # Go.
    # args.func is defined by each subprogram.
    return args.func(args)

    
# If we've been invoked as a program, call main().    
if __name__ == '__main__':
    sys.exit(main())
    