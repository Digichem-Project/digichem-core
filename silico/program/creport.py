#!/usr/bin/env python3

# The creport program.

# This should suppress a matplotlib warning when we compile with pyinstaller.
import warnings
warnings.filterwarnings("ignore", "(?s).*MATPLOTLIBDATA.*", category = UserWarning)
# This one is caused by some minor bug when plotting graphs.
warnings.filterwarnings("ignore", "(?s).*Source ID .* was not found when attempting to remove it.*")
# This one is caused by the way tight_layout works. It's fine to ignore here, but not so good if we want to plot interactively.
warnings.filterwarnings("ignore", "(?s).*tight_layout: falling back to Agg renderer*", category = UserWarning)

import sys
# These	two suppress weasyprint warnings about incompatible libraries (which we ignore when freezing).
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
	warnings.filterwarnings("ignore", "(?s).*@font-face support needs Pango >= 1.38*", category = UserWarning)
	warnings.filterwarnings("ignore", "(?s).*There are known rendering problems and missing features with cairo < 1.15.4*", category = UserWarning)

# Init openbabel.
import silico
silico.init_obabel()

# General imports.
import argparse
from pathlib import Path

# Silico imports.
import silico.program
from silico.report.pdf import PDF_report
from silico.exception.base import Silico_exception


# Printable name of this program.
NAME = "Calculation Report Generator"
DESCRIPTION = "generate PDF reports from finished Gaussian calculations"
EPILOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))
USAGE = "%(prog)s [options] file.log [file.chk] [file.fchk] [-o report.pdf]"

def main():
	"""
	Main entry point for the creport program.
	"""
	# ----- Program init -----
	
	# First configure our argument reader.
	parser = argparse.ArgumentParser(
		description = DESCRIPTION,
		usage = USAGE,
		epilog = EPILOG)
	parser.add_argument("calculation_files", help = "calculation result files (.log, .chk, .fchk etc) to extract results from", nargs = "*", default = [])
	
	emission_group = parser.add_argument_group("emission energy", "options for specifying additional input files that allow for the calculation of relaxed emission energy")
	emission_group.add_argument("--adiabatic_ground", help = "path to a secondary result file that contains results for the ground state to be used to calculate the adiabatic emission energy", default = None)
	emission_group.add_argument("--adiabatic_excited", help = "path to a secondary result file that contains results for the excited state to be used to calculate the adiabatic emission energy", default = None)
	emission_group.add_argument("--adiabatic_state", help = "the excited state in 'adiabatic_excited' to use, either a number indicating the state (eg, '1' for lowest state) or the state label (eg, 'S(1)', 'T(1)'). If not given and 'adiabatic_excited' contains excited states, the lowest will be used, otherwise the total energy of 'adiabatic_excited' will be used", default = None)
	emission_group.add_argument("--vertical_ground", help = "same as --adiabatic_ground, but for vertical emission", default = None)
	emission_group.add_argument("--vertical_excited", help = "same as --adiabatic_excited, but for vertical emission", default = None)
	emission_group.add_argument("--vertical_state", help = "same as --adiabatic_state, but for vertical emission", default = None)
							
	output_group = parser.add_mutually_exclusive_group()
	output_group.add_argument("--pdf_file", help = "a filename/path to a pdf file to write to (this is an alternative to the 'output' option). Other output files will placed in the same directory as the 'pdf_file'", nargs = "?", default = None)
	output_group.add_argument("-o", "--output", help = "a filename/path to write to. If this filename ends in a .pdf extension, this is taken as the filename of the pdf file to write to. Otherwise, output is assumed to be the name of a directory and the pdf file is written inside this directory", nargs = "?", type = Path, default = None)
	
	parser.add_argument("--type", help = "the type of report to make", choices = ["full", "atoms"], default = "full")
	parser.add_argument("-v", "--version", action = "version", version = str(silico.version))
	
	aux_input_group = parser.add_argument_group("additional input", "options for specifying additional input files used to generate the report. These options are all optional. Note that these input files can be given as positional arguments, in which case their file type is determined automatically. Use these named arguments to explicitly set the file type")
	aux_input_group.add_argument("--log_file", help = "a Gaussian log file that results will be read from", nargs = "?", default = None)
	aux_input_group.add_argument("--chk_file", help = "a Gaussian chk file that will be used to generate all image files required. Note that this option requires Gaussian to be installed and formchk & cubegen to be in your path", nargs = "?", default = None)
	aux_input_group.add_argument("--fchk_file", help = "a Gaussian fchk file that will be used to generate all image files required. Note that this option requires Gaussian to be installed and cubegen to be in your path", nargs = "?", default = None)
	
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
				'auto_crop': args.dont_auto_crop
			},
			'image': {
				'dont_create_new': args.dont_create_new_images,
				'use_existing': not args.overwrite_existing_images if args.overwrite_existing_images is not None else None
			}
		})
	
	# Use our generic init function.
	args, config, logger = silico.program.init_program(
		arg_parser = parser,
		arg_to_config = arg_to_config)
	
	if args.alignment is not None and not args.overwrite_existing_images:
		logger.warning("Alignment method has been changed but not overwriting existing images; use '-OK method' to ensure molecule images are re-rendered to reflect this change")
	
	
	# ----- Program begin -----
	return silico.program.main_wrapper(
		args = args,
		config = config,
		logger = logger,
		inner_func = _main
	)

def _main(args, config, logger):
	"""
	Inner portion of main (wrapped by a try-catch-log hacky boi).
	"""
	try:
		report = PDF_report.from_calculation_files(
			*args.calculation_files,
			gaussian_log_file = args.log_file,
			chk_file_path = args.chk_file,
			fchk_file_path = args.fchk_file,
			prog_version = silico.version,
			adiabatic_emission_ground_result = args.adiabatic_ground,
			adiabatic_emission_excited_result = args.adiabatic_excited,
			adiabatic_emission_excited_state = args.adiabatic_state,
			vertical_emission_ground_result = args.vertical_ground,
			vertical_emission_excited_result = args.vertical_excited,
			vertical_emission_excited_state = args.vertical_state,
			options = config
		)
	except Exception as e:
		raise Silico_exception("Failed to load results")
		
	
	# The filename (excluding parent directories) we'll write to if the user doesn't give us one.
	if args.type == "full":
		default_base_name = Path(report.gaussian_log_file.name).with_suffix(".pdf").name
	else:
		default_base_name = Path(report.gaussian_log_file.name).with_suffix(".atoms.pdf").name
	
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
		args.pdf_file = Path(report.gaussian_log_file.parent, "report", default_base_name)
		
	# Now make (compile?) and save our report.
	logger.info("Generating report '{}'".format(args.pdf_file))
	try:
		report.write(args.pdf_file, report_type = args.type)
	except Exception as e:
		raise Silico_exception("Failed to generate report '{}'; {}".format(args.pdf_file, str(e)))
		#logger.error("Failed to generate report '{}'; {}".format(args.pdf_file, str(e)), exc_info = True)
		#return -1
	
	# Done.
	logger.info("Done generating report '{}'".format(args.pdf_file))

	
# If we've been invoked as a program, call main().	
if __name__ == '__main__':
	sys.exit(main())
