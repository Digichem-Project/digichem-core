# The silico file conversion program.
# This program is largely a wrapper around openbabel, with a few additional formats supported.

# General imports.


# Silico imports.
import silico.program
from silico.file.convert import Silico_input
from silico.misc.file_wrapper import Multi_file_wrapper
from silico.file.convert.babel import Openbabel_converter
from silico.exception.base import Silico_exception
from logging import getLogger
from silico.misc.base import to_bool

# Printable name of this program.
NAME = "Calculation File Converter"
DESCRIPTION = "convert computational files between formats"
EPILOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))


def arguments(subparsers):
	"""
	Add this program's arguments to an argparser object.
	"""
	parser = subparsers.add_parser(
		'convert',
		description = DESCRIPTION,
		parents = [silico.program.standard_args],
		epilog = EPILOG,
		help = "Convert file formats")
	# Set main function.
	parser.set_defaults(func = main)
	
	parser.add_argument("input_file", help = "Input file to read and convert.")
	parser.add_argument("-i", "--input_format", help = "Input format (.com, .xyz, .tmol etc)")
	parser.add_argument("-o", "--output_format", help = "Output format (.com, .xyz, .tmol etc)")
	parser.add_argument("-O", "--output_file", help = "Output file", default = "-")
	parser.add_argument("-C", "--charge", help = "The molecular charge to set in the output format. Note that not all formats support a charge", default = None, type = float)
	parser.add_argument("-M", "--multiplicity", help = "The multiplicity to set in the output format. Note that not all formats support a multiplicity", default = None, type = int)
	parser.add_argument("--gen3D", help = "Whether to generate 3D coordinates (this will scramble existing atom coordinates). The default is yes, but only if it can be safely determined that the loaded coordinates are not already in 3D)", type = to_bool , default = None)
	
def main(args):
	"""
	Main entry point for the resume program.
	"""
	silico.program.main_wrapper(_main, args = args)

def _main(args, config, logger):
	"""
	Inner portion of main (wrapped by a try-catch-log hacky boi).
	"""
	# Load the file we were given.
	parser = Silico_input.from_file(args.input_file, args.input_format, charge = args.charge, multiplicity = args.multiplicity, gen3D = args.gen3D)
		
	# If we weren't given an output format, try and guess one.
	if args.output_format is None:
		try:
			args.output_format = Openbabel_converter.type_from_file_name(args.output_file)
		except Silico_exception:
			# Couldn't guess a format, we'll just assume .si
			getLogger(silico.logger_name).info("No output format given and could not guess from file name; using .si as default ", exc_info = True)
			args.output_format = "si"
	
	# Convert and write.
	with Multi_file_wrapper(args.output_file, "wt") as outfile:
		outfile.write(parser.to_format(args.output_format))
		
		