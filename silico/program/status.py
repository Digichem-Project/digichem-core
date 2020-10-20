# The creport program.

# General imports.
from pathlib import Path

# Silico imports.
import silico.program
from silico.report.pdf import PDF_report
from silico.exception.base import Silico_exception
from silico.result.excited_states import Excited_state_list


# Printable name of this program.
NAME = "Method status"
DESCRIPTION = "check status of known submission methods"
EPILOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))

def arguments(subparsers):
	"""
	Add this program's arguments to an argparser object.
	"""
	parser = subparsers.add_parser(
		'status',
		description = DESCRIPTION,
		parents = [silico.program.standard_args],
		epilog = EPILOG,
		help = "Check status")
	# Set main function.
	parser.set_defaults(func = main)
	
	parser.add_argument("methods", help = "Selected methods to show status for", nargs = "*", default = ())
	

def main(args):
	"""
	Main entry point for the resume program.
	"""
	silico.program.main_wrapper(_main, args = args)

def _main(args, config, logger):
	"""
	Inner portion of main (wrapped by a try-catch-log hacky boi).
	"""
	# Load our calculation definitions.
	try:
		config.resolve()
		known_methods = config.methods
	except Exception:
		raise Silico_exception("Failed to load calculation methods")
	
	# Now get the one's we've been asked to list.
	if len(args.methods) == 0:
		# None specified, use all.
		methods = known_methods
		
	else:
		methods = type(known_methods)([known_methods.get_config(method_id) for method_id in args.methods])
		
	for method in methods:
		try:
			print("{}: {}".format(method.NAME, method.status))
		except NotImplementedError:
			# No status for this method.
			print("{}: Status not available for this method".format(method.NAME))
		except Exception:
			logger.error("Failed to fetch status information for method '{}'".format(method.NAME), exc_info = True)
	
	