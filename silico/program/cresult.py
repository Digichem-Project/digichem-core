#!/usr/bin/env python3

# The cresult (Result Extraction) program.
# This program is designed to extract results from gaussian calculations and display/save them as more convenient formats.

# This should suppress a matplotlib warning when we compile with pyinstaller.
import warnings
warnings.filterwarnings("ignore", "(?s).*MATPLOTLIBDATA.*", category = UserWarning)
# This one is caused by some minor bug when plotting graphs.
warnings.filterwarnings("ignore", "(?s).*Source ID .* was not found when attempting to remove it.*")
# This one is caused by the way tight_layout works. It's fine to ignore here, but not so good if we want to plot interactively.
warnings.filterwarnings("ignore", "(?s).*tight_layout: falling back to Agg renderer*", category = UserWarning)

# Init openbabel.
import silico
silico.init_obabel()

# General imports.
import argparse
import cclib
from multiprocessing import Pool
from itertools import filterfalse
from functools import partial
import logging

# Silico imports.
import silico.program
from silico.misc.base import List_grouper
from silico.result.result import Result_set
from silico.misc.file_wrapper import Multi_file_wrapper
from silico.extract.text import Text_summary_group_extractor
from silico.extract.csv import CSV_summary_group_extractor, Long_CSV_group_extractor
from silico.extract.table import Table_summary_group_extractor, Long_table_group_extractor
from silico.exception.base import Silico_exception

# Printable name of this program.
NAME = "Calculation Result Extractor"
DESCRIPTION = "extract results from calculation output files and convert them to more convenient intermediate formats"
EPILOG = "{} V{}. Written by {}. Last updated {}.".format(NAME, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))
USAGE = """%(prog)s [options] file1.log [file2.log...] [-a %%FORMAT file.tbl] [-b %%FORMAT file.tbl] [-c %%FORMAT file.csv] [-d %%FORMAT file.csv] [-t %%FORMAT file.txt]
   eg: %(prog)s file1.log -t %%META %%GEOM %%ATOMS
   eg: %(prog)s *.log -c results.csv
   eg: %(prog)s *.log -a %%META %%ORB=+2,+1,+0,-1 | less -S"""

def list_known_sections(extractor_groups):
	"""
	Helper function, called when -l or --list is specified.
	"""
	for extractor_group in extractor_groups:
		# Get a list of sections.
		for extractor_class in extractor_group.recursive_subclasses():
			if hasattr(extractor_class, 'CLASS_HANDLE'):
				print(extractor_class.CLASS_HANDLE[0])
	
		
def _get_result_set(filename, **kwargs):
	"""
	Helper function, designed to be called in parallel to read input files.
	"""
	logger = logging.getLogger(silico.logger_name)
	# First open our file.
	try:
		with Multi_file_wrapper(filename, "rt") as input_file:
			logger.info("Parsing calculation result file '{}'".format(filename))
			return Result_set.from_cclib(
				cclib.io.ccread(input_file.file),
				gaussian_log_file = input_file.name if input_file.name != "-" else None,
				name = input_file.name,
				**kwargs)
	except Exception:
			logger.warning("Unable to parse calculation result file '{}'; skipping".format(filename), exc_info = True)

def _get_extractor_group_class(argument_group):
	"""
	Helper function, retrieves an extractor class.
	"""
	# TODO: We can ditch this function, just have the argparse action thingy set the class there and then.
	if argument_group is None or "-t" in argument_group.option_strings:
		return Text_summary_group_extractor
	elif "-c" in argument_group.option_strings:
		return CSV_summary_group_extractor
	elif "-d" in argument_group.option_strings:
		return Long_CSV_group_extractor
	elif "-a" in argument_group.option_strings:
		return Table_summary_group_extractor
	elif "-b" in argument_group.option_strings:
		return Long_table_group_extractor
	else:
		raise NotImplementedError("TODO")
	
def _split_file_names(names):
	"""
	Take a list of file names (as would be given on the command line) and determine which are files and which are sections.
	
	Sections are defined as starting with a '%' character and instruct the extractor which parts of the result set to extract.
	Sections can include the '=' character, with string following it being interpreted as filters for that section. Multiple filters are separated by commas (',').
	File names are everything else and are paths to files to write to. The name '-' (a single dash) is also recognised as stdout.
	
	:param names: The list of names to split.
	:return: A tuple of (file_names, sections). Each section is a dict of the form {'name', 'sub_criteria'}. Sub_criteria is a list of strings.
	"""
	# We're going to split our 'names' into those that are file names (paths) and those that represent extractors (which extract/write certain parts of our result object).
	file_names = []
	sections = []
	
	# Iterate through our list.
	for name in names:
		# The percent % character is what to look for.
		if name[:1] != "%":
			# This is a file name.
			file_names.append(name)
		elif name[:2] == "%%":
			# This is also a file name (double % is the escape mechanism, it means a file name that happens to start with a %).
			file_names.append(name[1:])
		else:
			# This is a section name.
			# Split on '=', the part before is the name of a section (which should match a CLASS_HANDLE), parts after are criteria to pass to that section class.
			equals_split = name.split("=", maxsplit = 1)
			
			# The name is the first part (remembering to remove the starting %.
			section_name = equals_split[0].strip()[1:]
			
			# The remainder are criteria, separated by commas (',').
			criteria = [[sub_criterion.strip() for sub_criterion in criterion.split(';')] for criterion in equals_split[1].split(',')] if len(equals_split) > 1 else [[]]
			
			# Add to our list of sections.
			sections.extend([{'name': section_name, 'sub_criteria': sub_criteria} for sub_criteria in criteria])
			
	return (file_names, sections)

def main():
	"""
	Main entry point for the cresult program.
	"""
	# ----- Program init -----
	
	# First configure out argument reader.
	parser = argparse.ArgumentParser(
		description = DESCRIPTION,
		usage = USAGE,
		epilog = EPILOG
	)
	parser.add_argument("calculation_files", help = "calculation result files (.log etc) to extract results from", nargs = "*", default = [])
	
	emission_group = parser.add_argument_group("emission energy", "options for specifying additional input files that allow for the calculation of relaxed emission energy")
	emission_group.add_argument("--adiabatic_ground", help = "path to a secondary result file that contains results for the ground state to be used to calculate the adiabatic emission energy", default = None)
	emission_group.add_argument("--adiabatic_excited", help = "path to a secondary result file that contains results for the excited state to be used to calculate the adiabatic emission energy", default = None)
	emission_group.add_argument("--adiabatic_state", help = "the excited state in 'adiabatic_excited' to use, either a number indicating the state (eg, '1' for lowest state) or the state label (eg, 'S(1)', 'T(1)'). If not given and 'adiabatic_excited' contains excited states, the lowest will be used, otherwise the total energy of 'adiabatic_excited' will be used", default = None)
	emission_group.add_argument("--vertical_ground", help = "same as --adiabatic_ground, but for vertical emission", default = None)
	emission_group.add_argument("--vertical_excited", help = "same as --adiabatic_excited, but for vertical emission", default = None)
	emission_group.add_argument("--vertical_state", help = "same as --adiabatic_state, but for vertical emission", default = None)
	
	parser.add_argument("-i", "--ignore", help = "ignore missing sections rather than stopping with an error", action = "store_true")
	parser.add_argument("-v", "--version", action = "version", version = str(silico.version))
	#parser.add_argument("-l", "--list", help = "list the available subsections for each requested output format and exit", action = "store_true")
	
	# Now process our output formats. The list of understood formats is generated automatically from the classes that inherit from silico.result.output.base.Output
	output_formats = parser.add_argument_group("output formats", "Specify output formats to write results to. Optionally give a number of filenames after the format (eg, -t output1.results output2.results) to write to those filenames, otherwise output is written to stdout. stdout can also be explicitly requested using the dash ('-') character (eg, -t -). Additionally, subsections can be requested by specifying the name of the section preceded by a percent ('%') character (eg, -t atoms.text %GEOM)")
	
	# Add our known formats
	output_formats.add_argument("-t", dest = "outputs", metavar = ("%FORMAT", "file.txt"), help = "human readable text format; shows various summaries of results for each result file", nargs = "*", default = [], action = List_grouper)
	output_formats.add_argument("-c", dest = "outputs", metavar = ("%FORMAT", "file.csv"), help = "CSV tabular format; shows one row per result; useful for comparing many results at once", nargs = "*", action = List_grouper)
	output_formats.add_argument("-d", dest = "outputs", metavar = ("%FORMAT", "file.csv"), help = "CSV summary format; shows a separate table for each property (atoms, MOs etc); one row for each item (atom, orbital etc)", nargs = "*", action = List_grouper)
	output_formats.add_argument("-a", dest = "outputs", metavar = ("%FORMAT", "file.tbl"), help = "tabulated text format; the same as -c but formatted with an ASCII table, recommended that output piped to 'less -S'", nargs = "*", action = List_grouper)
	output_formats.add_argument("-b", dest = "outputs", metavar = ("%FORMAT", "file.tbl"), help = "tabulated summary text format; the same as -d but formatted with an ASCII table", nargs = "*", action = List_grouper)
	
	# Use our generic init function.
	args, config, logger = silico.program.init_program(
		arg_parser = parser,
		arg_to_config = None)
	
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
	# First, process our args.
	extractor_groups = []
		
	# If we've not been given any outputs, use a default.
	if len(args.outputs) == 0:
		args.outputs.append({
			'group': None,
			'values': []
			})
	
	# We do this first because reading all the calc files can be very slow, and if there is an error in the extractors section then it will be all wasted.
	try:
		for argument in args.outputs:
			# This is the 'name' of the argument as given to argparse, eg 'text', 'table' etc,
			argument_group = argument['group']
			# Files is the list of options given to argparse, it is a mixture of file paths and section names ('%GEOM' for example)
			argument_files = argument['values']
			
			# Decide on which top-level class to use.
			# We'll store this in a dict for later.
			extractor_group = {}
			
			# We'll call from_class_handle() on this class to get a list of classes which will get our output for us.
			extractor_group_cls =_get_extractor_group_class(argument_group)
		
			# Split our 'files' into real files and sections.
			output_file_paths, output_sections = _split_file_names(argument_files)
			
			# Save the files for later.
			extractor_group['output_file_paths'] = output_file_paths
		
			# Now get a list of extractor classes that the user has asked for.
			extractors = [extractor_group_cls.from_class_handle(output_section['name'])(*output_section['sub_criteria'], ignore = args.ignore if args.ignore is True else None) for output_section in output_sections]
			
			# Now get our main extractor, which takes out list of sub extractors as an argument to its constructor, and save it for later.
			extractor_group['extractor'] = extractor_group_cls(*extractors, ignore = args.ignore if args.ignore is True else None)
			
			# Add to our list.
			extractor_groups.append(extractor_group)
	except Exception:
		raise Silico_exception("Failed to process extraction commands")
		#logger.error("failed to process extraction commands", exc_info = True)
		#return -1
	
	# If we were asked to list, do do.
	#if args.list:
	#	list_known_sections([group_dict['extractor'] for group_dict in extractor_groups])
	#	exit(0)
	
	# Get upset if we've got nothing to read.
	if len(args.calculation_files) == 0:
		raise Silico_exception("No calculation files to read (input files should appear before -a -b -c -d or -t options)")
		#logger.error("no calculation files to read (input files should appear before -a -b -c -d or -t options)")
		#return -1
		
	# Get our list of results. Do this in parallel.
	try:
		with Pool() as pool:
			results = list(
				filterfalse(lambda x: x is None,
					pool.map(
						partial(
							_get_result_set,
							alignment_class_name = config['alignment'],
							adiabatic_emission_ground_result = args.adiabatic_ground,
							adiabatic_emission_excited_result = args.adiabatic_excited,
							adiabatic_emission_excited_state = args.adiabatic_state,
							vertical_emission_ground_result = args.vertical_ground,
							vertical_emission_excited_result = args.vertical_excited,
							vertical_emission_excited_state = args.vertical_state,
						),
						args.calculation_files
					)
				)
			)
	except Exception:
		raise Silico_exception("Failed to read calculation files")
		#logger.error("failed to read calculation files", exc_info = True)
		#return -1
	
	# Now process.
	try:
		for extractor_group in extractor_groups:
			# If we don't have any files, add stdout as a default
			if len(extractor_group['output_file_paths']) == 0:
				extractor_group['output_file_paths'].append("-")	
			
			# Finally, write to each of the given paths.
			extractor_group['extractor'].write(results, extractor_group['output_file_paths'])	
	except Exception:
		raise Silico_exception("Failed to read calculation files")
		#logger.error("failed to extract results", exc_info = True)
		#return -1

# If we've been invoked as a program, call main().	
if __name__ == '__main__':
	main()
	