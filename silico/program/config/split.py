# Helper program for setting configurables.

# General imports.

# Silico imports.
import silico.program

# Printable name of this program.
NAME = "Silico config splitter"
DESCRIPTION = "split a number of config files into multiple children"
EPILOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))

def arguments(subparsers):
    """
    Add this program's arguments to an argparser object.
    """
    parser = subparsers.add_parser("split",
        description = DESCRIPTION,
        parents = [silico.program.standard_args],
        epilog = EPILOG,
        help = "Split config"
    )
    # Set main function.
    parser.set_defaults(func = main)
    
    parser.add_argument("configurables", help = "a (number of) config files to split into multiple children", nargs = "+")
    parser.add_argument("-o,", "--option", help = "the name of the option to split on", required = True)
    parser.add_argument("-v,", "--values", help = "a list of values to split", nargs = "*", default = [])
    
def main(args):
    """
    Main entry point for the resume program.
    """
    silico.program.main_wrapper(_main, args = args)
    
def _main(args, config, logger):
    """
    Inner portion of main (wrapped by a try-catch-log hacky boi).
    """