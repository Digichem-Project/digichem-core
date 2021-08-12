# Check the values of loaded config options are valid.

# General imports.

# Silico imports.
import silico.program

# Printable name of this program.
NAME = "Silico config validator"
DESCRIPTION = "Check the values of loaded config options."
EPILOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))

def arguments(subparsers):
    """
    Add this program's arguments to an argparser object.
    """
    parser = subparsers.add_parser("validate",
        description = DESCRIPTION,
        parents = [silico.program.standard_args],
        epilog = EPILOG,
        help = "Validate configs"
    )
    # Set main function.
    parser.set_defaults(func = main)
    
    parser.add_argument("type", help = "the type of configurables to validate. If none are given, all will be validated.", nargs = "*", choices = ("basis_sets", "methods", "programs", "calculations", "all"), default = "all")
    
def main(args):
    """
    Main entry point for the resume program.
    """
    silico.program.main_wrapper(_main, args = args)
        
def _main(args, config, logger):
    """
    Inner portion of main (wrapped by a try-catch-log hacky boi).
    """
    if args.type == "all":
        args.type = ("basis_sets", "methods", "programs", "calculations")

    for configurable_type in args.type:
        for index, configurable in enumerate(getattr(config, configurable_type)):
            configurable.validate()
            logger.info("Validated: {}) {}".format(index+1, configurable.name))
            
    logger.info("All configs validated")