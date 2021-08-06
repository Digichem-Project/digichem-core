# Helper program for setting configurables.

# General imports.

# Silico imports.
import silico.program.config.validate

# Printable name of this program.
NAME = "Silico config helper"
DESCRIPTION = "setup options for Silico"
EPILOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))

def arguments(subparser):
    """
    Add this program's arguments to an argparser object.
    """
    parser = subparser.add_parser("config",
        description = DESCRIPTION,
        epilog = EPILOG,
        help = "Setup options"
    )
    subparsers = parser.add_subparsers(dest="prog")
    
    # Create sub parsers for each sub-program. Each will define its own parser.
    silico.program.config.validate.arguments(subparsers)