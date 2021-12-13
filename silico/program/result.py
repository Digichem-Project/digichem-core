# General imports.
import os

# Silico imports.
from silico.program.base import Program
from silico.exception.base import Silico_exception
from silico.result.alignment.base import Alignment
from silico.format.text import Text_summary_group_format
from silico.format.csv import CSV_summary_group_format,\
    CSV_property_group_format
from silico.format.table import Table_summary_group_format,\
    Property_table_group_format
from silico.parser.base import parse_multiple_calculations
from silico.interface.urwid.result.base import Result_interface


class Result_program(Program):
    """
    The Silico result parsing program.
    """
    
    name = "Calculation Result Parser"
    command = "result"
    description = "extract results from calculation output files and convert them to more convenient intermediate formats"
    help = "Parse results"


    @classmethod
    def arguments(self, sub_parsers_object):
        """
        Add this program's arguments to an argparser subparsers object.
        
        :param sub_parsers_object: An argparser subparsers object to add to (as returned by add_subparsers().
        :returns: The argparser object (which supports add_argument()).
        """
        sub_parser = super().arguments(sub_parsers_object)
        
        sub_parser.add_argument("log_files", help = "a (number of) calculation result file(s) (.log) to extract results from", nargs = "*", default = [])
         
        sub_parser.add_argument("-i", "--ignore", help = "ignore missing sections rather than stopping with an error", action = "store_true", default = False)
        sub_parser.add_argument("-C", "--num_CPUs", help = "the number of CPUs to use in parallel to parse given log_files, defaults to the number of CPUs on the system", type = int, nargs = "?", default = os.cpu_count())
        
        output_group = sub_parser.add_argument_group("output format", "the format to write results to. Only one option from the following may be chosen")
        output_format = output_group.add_mutually_exclusive_group()
        output_format.add_argument("-t", "--text", help = "human readable text format; shows various summaries of results for each result file", dest = "format", action = "store_const", const = Text_summary_group_format, default = Text_summary_group_format)
        output_format.add_argument("-c", "--csv", help = "CSV tabular format; shows one row per result; useful for comparing many results at once", dest = "format", action = "store_const", const = CSV_summary_group_format)
        output_format.add_argument("-d", "--csv-property", help = "CSV property format; shows a separate table for each property (atoms, MOs etc); one row for each item (atom, orbital etc)", dest = "format", action = "store_const", const = CSV_property_group_format)
        output_format.add_argument("-a", "--table", help = "tabulated text format; the same as -c but formatted with an ASCII table, recommended that output be piped to 'less -S'", dest = "format", action = "store_const", const = Table_summary_group_format)
        output_format.add_argument("-b", "--table-property", help = "tabulated property text format; the same as -d but formatted with an ASCII table", dest = "format", action = "store_const", const = Property_table_group_format)
        
        sub_parser.add_argument("-p", "--properties", "--prop", help = "a list of properties (SCF, MOS, atoms etc) to parse", nargs = "*", default = [])
        
        #sub_parser.add_argument("-O", "--output", help = "a (number of) filename/path to write results to. If none is give, results will be written to stdout, which can also be explicitly requested with '-'", nargs = "*", default = ("-",))
        sub_parser.add_argument("-O", "--output", help = "a filename/path to write results to. If none is give, results will be written to stdout, which can also be explicitly requested with '-'", default = "-")
    
        return sub_parser
    
    def __init__(self, results, args, config, logger):
        """
        Constructor for submit programs.
        
        :param results: A list of result sets to format.
        """
        super().__init__(args = args, config = config, logger = logger)
        self.results = results

    @classmethod
    def load_from_argparse(self, args, config, logger):
        """
        Create a Program object from the data provided by argparse.
        
        :param args: The command-line arguments the program was started with.
        :param config: A loaded Silico config object.
        :param logger: The logger to use for output.
        """
        # Get our alignment class.
        alignment = Alignment.from_class_handle(config['alignment'])
        
        # Parse the log files we've been given (in parallel).
        # Get our list of results. Do this in parallel.
        try:
            results = parse_multiple_calculations(*args.log_files, alignment_class = alignment, init_func = self.subprocess_init, processes = args.num_CPUs)
            
        except Exception:
            raise Silico_exception("Failed to read calculation files")
                
        return self(results, args = args, config = config, logger = logger)
    
    def process_properties(self, properties_strings):
        """
        Process a list of properties strings, each of which identifies a formatter.
        
        :param properties_strings: A list of strings to process.
        :returns: A list of dictionaries of {'name', 'sub_criteria'}
        """
        # A list of dictionaries of {'name', 'sub_criteria'}
        requested_properties = []
        
        # Iterate through our list.
        for properties_string in properties_strings:
            # Split on '=', the part before is the name of a section (which should match a CLASS_HANDLE), parts after are criteria to pass to that section class.
            equals_split = properties_string.split("=", maxsplit = 1)
            
            # The name is the first part.
            property_name = equals_split[0].strip()
            
            # The remainder are criteria, separated by commas (',').
            criteria = [[sub_criterion.strip() for sub_criterion in criterion.split(';')] for criterion in equals_split[1].split(',')] if len(equals_split) > 1 else [[]]
            
            # Add to our list of sections.
            requested_properties.extend([{'name': property_name, 'sub_criteria': sub_criteria} for sub_criteria in criteria])
            
        return requested_properties
    
    def load_interface(self, window):
        """
        Get the interface widget we'll use for display.
        
        :param window The parent window object.
        """
        return Result_interface(window, self)

    def main(self):
        """
        Logic for our program.
        """
        # Get upset if we have no log files.
        if len(self.results) == 0:
            raise Silico_exception("No result files specified/available")
        
        # Also if we have no-where to write to.
        if self.args.output is None:
            raise Silico_exception("No output location specified")
        
        # Process our given properties, splitting each into the property name and associated filters.
        requested_properties = self.process_properties(self.args.properties)
        
        # First, get the individual parsers that the user asked for (these are based on the values given to the --properties option).
        formats = [self.args.format.from_class_handle(requested_property['name'])(*requested_property['sub_criteria'], ignore = self.args.ignore, config = self.config) for requested_property in requested_properties]
        
        # Now get our main formating object.
        output_format = self.args.format(*formats, ignore = self.args.ignore, config = self.config)
                
        # Now process.
        try:
            output_format.write(self.results, (self.args.output,))
            
        except Exception as e:
            raise Silico_exception("Failed to parse/write results") from e
