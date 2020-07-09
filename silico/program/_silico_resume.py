#!/usr/bin/env python3

# This should suppress a matplotlib warning when we compile with pyinstaller.
import warnings
from silico.submit.structure.flag import Flag
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
import numpy
import argparse
from pathlib import Path
import pickle
from logging import getLogger

# Silico imports.
import silico.program
from silico.submit.method import *
from silico.submit.program import *
from silico.submit.calculation import *
from silico.exception.uncatchable import Submission_paused
from silico.exception.base import Silico_exception

# This script is part of the csubmit resume mechanism. It is called as part of csubmit automatically, you would not normally run this script yourself.

# Printable name of this program.
NAME = "Silico csubmit resume"
DESCRIPTION = "resume an in-progress csubmit submission"
EPLIOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))

def main():
	"""
	Main entry point for the resume program.
	"""
	# ----- Program init -----
	
	# First configure our argument reader.
	parser = argparse.ArgumentParser(
		description = DESCRIPTION,
		epilog = EPLIOG)
	parser.add_argument("resume_file", help = "Path to the file to resume from (this should be a pickled calculation class)", type = Path)
	
	args = parser.parse_args()
		
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
	
	# ----- Program begin -----
# 	return silico.program.main_wrapper(
# 		_main,
# 		method = method,
# 		arg_parser = parser,
# 	)
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
	finally:
		method.calc_dir.set_flag(Flag.DONE)


# If we've been invoked as a program, call main().	
if __name__ == '__main__':
	sys.exit(main())