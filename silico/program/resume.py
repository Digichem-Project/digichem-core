# General imports.
import numpy
import argparse
from pathlib import Path
import pickle
from logging import getLogger

# Silico imports.
import silico.program
from silico.exception.uncatchable import Submission_paused
from silico.exception.base import Silico_exception

# This script is part of the silico resume mechanism. It is called as part of silico submit automatically, you would not normally run this script yourself.

# Printable name of this program.
NAME = "Silico submission resume"
DESCRIPTION = "resume an in-progress silico submission"
EPILOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))

def arguments(subparsers):
	"""
	Add this program's arguments to an argparser object.
	"""
	parser = subparsers.add_parser("resume",
		description = DESCRIPTION,
		epilog = EPILOG,
		help = "Resume submission (used automatically by part of the submission mechanism)"
	)
	parser.add_argument("resume_file", help = "Path to the file to resume from (this should be a pickled calculation class)", type = Path)

def main(args):
	"""
	Main entry point for the resume program.
	"""
	# ----- Program init -----
		
	# Load the pickled class.
	with open(args.resume_file, "rb") as pickle_file:
		method = pickle.load(pickle_file)

	# Delete the pickled file to clean up (also prevents us running twice on the same file, which would be bad. Maybe we should do some file locking anyway?)
	try:
		args.resume_file.unlink()
	except Exception:
		getLogger(silico.logger_name).error("Failed to delete pickle file", exc_info = True)
		
	# Set program stuff.
	silico.program.init_from_config(getLogger(silico.logger_name), method.program.calculation.silico_options)
	silico.program.init_signals(getLogger(silico.logger_name))
	
	# Set numpy errors (not sure why this isn't the default...)
	numpy.seterr(invalid = 'raise', divide = 'raise')
	
	return silico.program.run(_main, method = method)


def _main(method):
	"""
	Main entry point for the resume program.
	"""
	# And resume.
	try:
		method.resume()
	except Submission_paused:
		# This is fine.
		pass
	except Exception:
		raise Silico_exception("Error during submission")

