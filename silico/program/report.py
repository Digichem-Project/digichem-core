# The creport program.

# General imports.
from pathlib import Path

# Silico imports.
import silico.program
from silico.exception.base import Silico_exception
import silico.report
from silico.parser import parse_calculations
import itertools

# Printable name of this program.
NAME = "Calculation Report Generator"
DESCRIPTION = "generate PDF reports from finished Gaussian calculations"
EPILOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))
USAGE = "%(prog)s [options] file.log [file.chk] [file.fchk] [-o report.pdf]"

def arguments(subparsers):
    """
    Add this program's arguments to an argparser object.
    """
    parser = subparsers.add_parser("report",
        description = DESCRIPTION,
        parents = [silico.program.standard_args],
        usage = USAGE,
        epilog = EPILOG,
        help = "Write reports"
    )
    # Set main function.
    parser.set_defaults(func = main)
    
    parser.add_argument("log_files", help = "a (number of) calculation result file(s) (.log) to extract results from", nargs = "+")

    
    emission_group = parser.add_argument_group("emission energy", "options for specifying additional input files that allow for the calculation of relaxed emission energy")
    emission_group.add_argument("--adiabatic_ground", help = "path to a secondary result file that contains results for the ground state to be used to calculate the adiabatic emission energy", default = None)
    emission_group.add_argument("--vertical_ground", help = "same as --adiabatic_ground, but for vertical emission", default = None)
    emission_group.add_argument("--emission", help = "path to a secondary result file that contains results for the excited state to be used to calculate adiabatic and/or vertical emission energy", default = None)
    emission_group.add_argument("--emission_state", help = "the excited state in 'emission' to use, either a number indicating the state (eg, '1' for lowest state) or the state label (eg, 'S(1)', 'T(1)').", default = None)
                            
    output_group = parser.add_mutually_exclusive_group()
    output_group.add_argument("--pdf_file", help = "a filename/path to a pdf file to write to (this is an alternative to the 'output' option). Other output files will placed in the same directory as the 'pdf_file'", nargs = "?", default = None)
    output_group.add_argument("-o", "--output", help = "a filename/path to write to. If this filename ends in a .pdf extension, this is taken as the filename of the pdf file to write to. Otherwise, output is assumed to be the name of a directory and the pdf file is written inside this directory", nargs = "?", type = Path, default = None)
    
    parser.add_argument("--name", help = "name of the molecule/system to use in the report", default = None)
    parser.add_argument("--type", help = "the type of report to make", choices = ["full", "atoms"], default = "full")
    
    aux_input_group = parser.add_argument_group("auxiliary input", "options for specifying additional input files. These options are all optional, but if given the ordering should match that of the given log files")
    aux_input_group.add_argument("--chk", help = "a (number of) Gaussian chk file(s) that will be used to generate all image files required. Note that this option requires Gaussian to be installed and formchk & cubegen to be in your path", nargs = "*", default = [])
    aux_input_group.add_argument("--fchk", help = "a (number of) Gaussian fchk file(s) that will be used to generate all image files required. Note that this option requires Gaussian to be installed and cubegen to be in your path", nargs = "*", default = [])
    
    image_group = parser.add_argument_group("image options", "advanced options for altering the appearance of images in the report")
    image_group.add_argument("--render_style", help = "change the style used to render molecules and orbitals", choices=['pastel', 'light-pastel', 'dark-pastel','sharp', 'gaussian', 'vesta'], default = None)
    image_group.add_argument("-D", "--dont_create_new_images", help = "don't attempt to create any new image files. If existing image files can be found, they will still be included in the report", action = "store_true", default = None)
    image_group.add_argument("-O", "--overwrite_existing_images", help = "force overwriting and re-rendering of all image files", action = "store_true", default = None)
    image_group.add_argument("--dont_auto_crop", help = "disable automatic cropping of excess whitespace around the border of images. If given, rendered images will have the same resolution, but molecules may only occupy a small portion of the true image", action = "store_false", default = None)
    
    # 'Temporary' function that maps our custom command-line arguments to our config object.
    def arg_to_config(args, config):
        # Add our report specific command line arguments to our config dictionary.
        config.add_config({
            'molecule_image': {
                'rendering_style': args.render_style,
                'auto_crop': args.dont_auto_crop,
                'use_existing': not args.overwrite_existing_images if args.overwrite_existing_images is not None else None
            },
            'image': {
                'dont_modify': args.dont_create_new_images
            }
        })
    

def main(args):
    """
    Main entry point for the resume program.
    """
    silico.program.main_wrapper(_main, args = args)

def _main(args, config, logger):
    """
    Inner portion of main (wrapped by a try-catch-log hacky boi).
    """
    config.resolve()
    
    if args.alignment is not None and not args.overwrite_existing_images:
        logger.warning("Alignment method has been changed but not overwriting existing images; use '-OK method' to ensure molecule images are re-rendered to reflect this change")
    
    named_input_files = {
        "chk_file_path": args.chk,
        "fchk_file_path": args.fchk,
        "adiabatic_emission_ground_result": args.adiabatic_ground,
        "vertical_emission_ground_result": args.vertical_ground,
        "emission_excited_result": args.emission
    }
    
    try:
        # Transpose our lists of aux inputs.
        aux_files = [{"chk_file":chk_file, "fchk_file":fchk_file} for chk_file, fchk_file in list(itertools.zip_longest(args.chk, args.fchk))]
        
        # First, load results.
        result = parse_calculations(*args.log_files, aux_files = aux_files)
        
        # Then get an appropriate report.
        report = silico.report.from_result(result, options = config)
        
#         report = silico.report.from_files(
#             args.log_file,
#             options = config,
#             named_input_files = named_input_files,
#             emission_excited_state = args.emission_state
#         )
    except Exception as e:
        raise Silico_exception("Failed to load results")
    
    log_file = args.log_files[0]
    if log_file is not None:
        input_name = Path(log_file)
    elif len(args.calculation_files) > 0:
        input_name = Path(args.calculation_files[0])
    else:
        input_name = "Report"
    
    # The filename (excluding parent directories) we'll write to if the user doesn't give us one.
    if args.type == "full":
        default_base_name = input_name.with_suffix(".pdf").name
    else:
        default_base_name = input_name.with_suffix(".atoms.pdf").name
    
    # Decide on our output file.
    if args.pdf_file is not None:
        # We've been given a specific file to write to. No problems.
        args.pdf_file = Path(args.pdf_file)
    elif args.output is not None:
        # We've been given the more general 'output' option, we need to decide whether this aught to be a file or directory.
        if args.output.suffix == ".pdf":
            # It's a file! Set pdf_file appropriately.
            args.pdf_file = args.output
        else:
            # It's a directory. Decide on a fitting name from our input file.
            args.pdf_file = Path(args.output, default_base_name)
    else:
        # No output was given, we'll use a default location (in the same place as the main input file).
        args.pdf_file = Path(input_name.parent, "report", default_base_name)
        
    # Now make (compile?) and save our report.
    logger.info("Generating report '{}'".format(args.pdf_file))
    try:
        report.write(args.pdf_file, report_type = args.type)
    except Exception as e:
        raise Silico_exception("Failed to generate report '{}'".format(args.pdf_file)) from e
    
    # Done.
    logger.info("Done generating report '{}'".format(args.pdf_file))

