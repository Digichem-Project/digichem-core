#!/usr/bin/env python3.6

# The csubmit program (this is the successor to gsubmit).

# This should suppress a matplotlib warning when we compile with pyinstaller.
import warnings
import copy
from silico.interface.urwid.tree import Calculation_browser
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
import readline
import shlex
import textwrap

# Silico imports.
import silico.program
from silico.submit.method import *
from silico.submit.program import *
from silico.submit.calculation import *
from silico.submit.basis import *
from silico.exception import Silico_exception, Configurable_exception
from silico.exception.uncatchable import Submission_paused
from silico.misc.node_printer import Node_printer, Node
import silico.misc.tree


# Printable name of this program.
NAME = "Silico Calculation Submitter"
DESCRIPTION = "submit calculation files"
EPILOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))
USAGE = """%(prog)s [options] -i
   or: %(prog)s [options] file.com [file2.com ...] -i
   or: %(prog)s [options] file.com [file2.com ...] -c method/program/calc [method2/program2/calc2 ...]"""


def _list(methods, programs, calculations, all = False):
	"""
	Helper function that list all known calculations.
	"""
	# Get our nodes.
	top_node = Node.from_list((methods, programs, calculations))
	
	print('\033[1m' + "-" * 80 + '\033[0m')
	print('\033[1m' + " " * 27 +"Available Calculations:" + '\033[0m')
	print('\033[1m' + "-" * 80 + '\033[0m')
	print(Node_printer((top_node,)).get(all = all))
	print("-" * 80)
	print("Calculations are selected by specifying the 3 relevant numbers separated by a slash (eg, 1/1/1)")
	print("The first two numbers default to '1' if not given, so '3' and '1/3' are equivalent to '1/1/3'")
	print("Multiple calculations are selected by separating each with a space (eg, 1/1/1 1/1/2 1/1/3...)")
	print("-" * 80)
				
def _get_input(prompt):
	"""
	Helper function that gets input from the user.
	"""
	return input(prompt)
	
def _parse_calc_string(calc_string, methods, programs, calculations):
	"""
	Parse a calculation request string, of the format M/P/C, P/C or C.
	
	:raises IndexError: If the calculation string is invalid.
	:return: The relevant method, program and calculation.
	"""
	split_string = calc_string.split("/")
	
	# We parse backwards.
	# The final part is the calculation.
	if len(split_string) > 0:
		#calculation = Calculation_target.from_name_in_list(split_string[-1], calculations)
		calculation = calculations.get_config(split_string[-1])
	else:
		raise Silico_exception("'{}' is not a valid calculation identifier".format(calc_string))
	
	if len(split_string) > 1:
		#program = Program_target.from_name_in_list(split_string[-2], programs)
		program = programs.get_config(split_string[-2])
	else:
		try:
			#program = Program_target.from_name_in_list(calculation.programs[0], programs)
			program = programs.get_config(calculation.programs[0])
		except IndexError:
			raise Configurable_exception(calculation, "calculation has no programs set")
		
	if len(split_string) > 2:
		#method = Method_target.from_name_in_list(split_string[-3], methods)
		method = methods.get_config(split_string[-3])
	else:
		try:
			#method = Method_target.from_name_in_list(program.methods[0], methods)
			method = methods.get_config(program.methods[0])
		except IndexError:
			raise Configurable_exception(program, "program has no methods set")
	
	# We need to deepcopy our configurables incase we re-use one of them.
	return (copy.deepcopy(method), copy.deepcopy(program), copy.deepcopy(calculation))

def _get_warning_confirmation(method):
	"""
	Helper function that prompts the user for confirmation when they choose a method with a usage warning.
	"""
	print('\033[1m' +'\033[91m' + "-" * 80 + '\033[0m')
	print(" " * 33 + '\033[1m' + '\033[91m' + "!!! WARNING !!!" + '\033[0m')
	print('\033[1m' +'\033[91m' + "-" * 80 + '\033[0m')
	print("\n".join(textwrap.wrap(method.warning, 80)))
	return True if input('\033[1m' + "Are you sure?" + '\033[0m' + " (y/N): ").lower() == "y" else False

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
	parser.add_argument("calculation_files", help = "Gaussian input files to submit", nargs = "*", type = Path)
	parser.add_argument("-o", "--output", help = "Base directory to perform calculations in. Defaults to the current directory", default = Path("./"))
	parser.add_argument("-c", "--calculations", help = "Calculations to perform, identified either by name or by ID. To use a method and/or program other than the default, use the format M/P/C (eg, 2/1/1)", nargs = "*", default = [])
	parser.add_argument("-l", "--list", help = "List all known calculations; give twice for more output", action = "count", default = 0)
	parser.add_argument("-i", "--interactive", help = "Run in interactive mode, prompting for missing input", action = "store_true")
		
	# ----- Program begin -----
	return silico.program.main_wrapper(
		_main,
		arg_parser = parser
	)

def _main(args, config, logger):
	"""
	Inner portion of main (wrapped by a try-catch-log hacky boi).
	"""
	# Set tab completion to on.
	readline.parse_and_bind("tab: complete")
	
	# Load our calculation definitions.
	try:
		config.resolve()
		known_methods = config.methods
		known_programs = config.programs
		known_calculations = config.calculations
	except Exception:
		raise Silico_exception("Failed to load calculations")
	
	logger.debug("Loaded {} methods, {} programs and {} calculations".format(len(known_methods), len(known_programs), len(known_calculations)))
	
	# Print our list if we've been asked to.
	if args.list > 0:
		_list(known_methods, known_programs, known_calculations, all = args.list > 1)
		return 0
	
	# Check to see if we've got any files to submit (and prompt if we can).
	while len(args.calculation_files) == 0 and args.interactive:
		# Prompt for some files.
		print("Specify a list of files to submit (separate each with a space)")
		args.calculation_files = [Path(calc_file) for calc_file in shlex.split(_get_input("Files: "))]
		
	# Get upset if we have no files.
	if len(args.calculation_files) == 0:
		raise Silico_exception("No files to submit (use the -i option for interactive mode)")
		#logger.error("No files to submit (use the -i option for interactive mode)")
		#return 1
	
	
	calculations = []
	while True:
		try:
			# Parse our calculation strings, getting the actual calculation objects.
			calculations = [_parse_calc_string(calc_string, known_methods, known_programs, known_calculations) for calc_string in args.calculations]
		except:
			raise Silico_exception("Failed to parse calculations to submit")
			#logger.error("Failed to parse calculations to submit", exc_info = True)
			#exit(1)
						
		# If a 'dangerous' (one with a warning set) method has been chosen (and we're interactive), get confirmation.
		for method, program, calculation in calculations:
			if args.interactive and method.warning is not None:
				# Get confirm.
				if not _get_warning_confirmation(method):
					# User said no, reset out list so we'll ask them again.
					args.calculations = []
					calculations = []
		
		# Quit the loop if we have some calculations or if we can't ask the user for some.
		if len(calculations) != 0 or not args.interactive:
			break
		
		# Ask the user for calcs.
		#args.calculations = shlex.split(silico.misc.tree.run(Node.from_list((known_methods, known_programs, known_calculations))))
		args.calculations = shlex.split(Calculation_browser.run(Node.from_list((known_methods, known_programs, known_calculations))))
	
	# Get upset again if we have no calculations.
	if len(args.calculations) == 0:
		raise Silico_exception("No calculations to submit (use the -i option for interactive mode)")
		#logger.error("No calculations to submit (use the -i option for interactive mode)")
		#return 1
		
	
	try:
		# Arrange our calcs into a linked list.
		calculation = Calculation_target.prepare_list(calculations)
	except Exception:
		raise Silico_exception("Error processing calculations to submit")
		#logger.error("Error processing calculations to submit", exc_info = True)
		#return 1
	
	# The number of calcs we successfully submitted.
	done = 0
	
	for input_file_path in args.calculation_files:
		try:
			calculation.submit_from_file(args.output, input_file_path)
		except Submission_paused:
			# This is fine.
			pass
		except Exception:
			# Something went wrong.
			# We don't stop here though, we might have more calcs we can submit.
			logger.error("Failed to submit file '{}'".format(input_file_path), exc_info = True)
		done += 1
	
	return done
	
	
# If we've been invoked as a program, call main().	
if __name__ == '__main__':
	sys.exit(main())
	