# Program to extract results and import into a partially populated excel-like spreadsheet.

# General imports.

# Silico imports.
import silico.program

# Printable name of this program.
NAME = "Calculation to Spreadsheet Converter"
DESCRIPTION = "extract results from calculation output files and insert them into partially populated spreadsheets"
EPILOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))

def arguments(subparsers):
	"""
	Add this program's arguments to an argparser object.
	"""
	parser = subparsers.add_parser(
		"spreadsheet",
		aliases = ["sheet", "excel"],
		description = DESCRIPTION,
		parents = [silico.program.standard_args],
		epilog = EPILOG,
		help = "Extract results into spreadsheets"
	)
	# Set main function.
	parser.set_defaults(func = main)
	
	parser.add_argument("calculation_files", help = "calculation result file (.log etc) to extract results from")
		
	parser.add_argument("-s", "--spreadsheet", help = "path to a spreadsheet file that will be populated with data", default = None)
	parser.add_argument("-w", "--worksheet", help = "name of the worksheet (in spreadsheet) to populate with data", default = "Data")
	

def main(args):
	"""
	Main entry point for the resume program.
	"""
	silico.program.main_wrapper(_main, args = args)

def _main(args, config, logger):
	"""
	Inner portion of main (wrapped by a try-catch-log hacky boi).
	"""
	

	