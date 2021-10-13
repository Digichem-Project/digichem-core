# The csubmit program (this is the successor to gsubmit).

# General imports.
from pathlib import Path
import readline
import shlex
import copy

# Silico imports.
import silico.program
from silico.submit.calculation import Calculation_target
from silico.exception import Silico_exception, Configurable_exception
from silico.exception.uncatchable import Submission_paused
#from silico.interface.urwid.tree import Calculation_browser
from silico.misc.base import to_bool
from silico.interface.urwid.base import Window
from silico.interface.urwid.submit.base import Calculation_submitter



# Printable name of this program.
NAME = "Silico Calculation Submitter"
DESCRIPTION = "submit calculation files"
EPILOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))
USAGE = """%(prog)s [options]
   or: %(prog)s [options] file.com [file2.com ...]
   or: %(prog)s [options] file.com [file2.com ...] -c destination/program/calc [destination2/program2/calc2 ...]"""

def _get_input(prompt):
    """
    Helper function that gets input from the user.
    """
    return input(prompt)
    
def _parse_calc_string(method_string, destinations, programs, calculations):
    """
    Parse a calculation request string, of the format M/P/C, P/C or C.
    
    :raises IndexError: If the calculation string is invalid.
    :return: The relevant destination, program and calculation.
    """
    split_string = method_string.split("/")
    
    # We parse backwards.
    # The final part is the calculation.
    if len(split_string) > 0:
        calculation = calculations.get_config(split_string[-1])
    else:
        raise Silico_exception("'{}' is not a valid calculation identifier".format(method_string))
    
    if len(split_string) > 1:
        program = programs.get_config(split_string[-2])
    else:
        try:
            program = programs.get_config(calculation.programs[0])
        except IndexError:
            raise Configurable_exception(calculation, "calculation has no programs set")
        
    if len(split_string) > 2:
        destination = destinations.get_config(split_string[-3])
    else:
        try:
            destination = destinations.get_config(program.destinations[0])
        except IndexError:
            raise Configurable_exception(program, "program has no destinations set")
    
    # We need to deepcopy our configurables in case we re-use one of them.
    return (copy.deepcopy(destination), copy.deepcopy(program), copy.deepcopy(calculation))


def arguments(subparsers):
    """
    Add this program's arguments to an argparser object.
    """
    parser = subparsers.add_parser(
        'submit',
        description = DESCRIPTION,
        parents = [silico.program.standard_args],
        usage = USAGE,
        epilog = EPILOG,
        help = "Submit calculations")
    # Set main function.
    parser.set_defaults(func = main)
    
    parser.add_argument("calculation_files", help = "Gaussian input files to submit", nargs = "*", type = Path)
    parser.add_argument("-o", "--output", help = "Base directory to perform calculations in. Defaults to the current directory", default = Path("./"))
    parser.add_argument("-m", "--methods", help = "Methods to perform, identified either by name or by ID. To use a destination and/or program other than the default, use the format M/P/C (eg, 2/1/1)", nargs = "*", default = [])
    parser.add_argument("-C", "--charge", help = "Set the molecular charge of all input files. Note that certain calculations will override this value", default = None, type = float)
    parser.add_argument("-M", "--multiplicity", help = "Set the multiplicity of all input files. Note that certain calculations will override this value", default = None, type = int)
    parser.add_argument("--gen3D", help = "Whether to generate 3D coordinates (this will scramble existing atom coordinates). The default is yes, but only if it can be safely determined that the loaded coordinates are not already in 3D)", type = to_bool, default = None)

def main(args):
    """
    Main entry point for the resume program.
    """
    silico.program.main_wrapper(_main, args = args)

def _main(args, config, logger):
    """
    Inner portion of main (wrapped by a try-catch-log hacky boi).
    """
    # Set tab completion to on.
    readline.parse_and_bind("tab: complete")
    
    # Check to see if we've got any files to submit (and prompt if we can).
    while len(args.calculation_files) == 0:
        # Prompt for some files.
        print("Specify a list of files to submit (separate each with a space)")
        args.calculation_files = [Path(calc_file) for calc_file in shlex.split(_get_input("Files: "))]
        
    # Get upset if we have no files.
    if len(args.calculation_files) == 0:
        raise Silico_exception("No files to submit")
    
    
    window = Window(Calculation_submitter(config.methods), title = "Silico {}".format(silico.version))
    window.run_loop(config.urwid_palette)
    exit()
    
    methods = []
    while True:
        try:
            # Parse our calculation strings, getting the actual calculation objects.
            methods = [_parse_calc_string(method_string, known_destinations, known_programs, known_calculations) for method_string in args.methods]
        except:
            raise Silico_exception("Failed to parse calculations to submit")
                            
        # Quit the loop if we have some methods or if we can't ask the user for some.
        if len(methods) != 0:
            break
        
        # Ask the user for calcs.
        Calculation_browser(config.methods, config.urwid_palette).main()
        # args.methods = shlex.split(Calculation_browser.run(known_destinations, known_programs, known_calculations, Calculation_browser.yaml_to_palette(config['palette'])))
        if len(args.methods) == 0:
            break
        
    # Print any warning messages (but only once each).
    warnings = [None]
    for configurables in methods:
        for configurable in configurables:
                if configurable.warning not in warnings:
                    logger.warning(configurable.warning)
                    
    
    # Get upset again if we have no calculations.
    if len(args.calculations) == 0:
        raise Silico_exception("No calculations to submit")
    
    try:
        # Arrange our calcs into a linked list.
        first = Calculation_target.link(methods, global_silico_options = config)
    except Exception:
        raise Silico_exception("Error processing calculations to submit")
    
    # The number of calcs we successfully submitted.
    done = 0
    
    for input_file_path in args.calculation_files:
        try:
            # Prepare.
            first.prepare_from_file(args.output, input_file_path, molecule_charge = args.charge, molecule_multiplicity = args.multiplicity, gen3D = args.gen3D)
            
            # Go.
            first.submit()
            
            done += 1
            
        except Submission_paused:
            # This is fine.
            done += 1
        
        except Exception:
            # Something went wrong.
            # We don't stop here though, we might have more calcs we can submit.
            logger.error("Failed to submit file '{}'".format(input_file_path), exc_info = True)
            
        
    
    return done
    