# The silico file conversion program.
# This program is largely a wrapper around openbabel, with a few additional formats supported.

# General imports.

# Silico imports.
from silico.program.base import Program
from silico.exception.base import Silico_exception
from silico.misc.base import to_bool
from silico.input import si_from_file
from silico.file.babel import Openbabel_converter
from silico.misc.file_wrapper import Multi_file_wrapper
from silico.interface.urwid.convert import Convert_interface


class Convert_program(Program):
    """
    The Silico format conversion program.
    """
    
    name = "Calculation File Converter"
    command = "convert"
    aliases = ["c", "con"]
    description = "convert computational files between formats"
    help = "Convert files"
    
    @classmethod
    def arguments(self, sub_parsers_object):
        """
        Add this program's arguments to an argparser subparsers object.
        
        :param sub_parsers_object: An argparser subparsers object to add to (as returned by add_subparsers().
        :returns: The argparser object (which supports add_argument()).
        """
        sub_parser = super().arguments(sub_parsers_object)
        
        sub_parser.add_argument("input_file", help = "Input file to read and convert.", nargs = "?")
        sub_parser.add_argument("-i", "--input_format", help = "Input format (.com, .xyz, .tmol etc)")
        sub_parser.add_argument("-o", "--output_format", help = "Output format (.com, .xyz, .tmol etc)")
        sub_parser.add_argument("-O", "--output_file", help = "Output file", default = "-")
        sub_parser.add_argument("-C", "--charge", help = "The molecular charge to set in the output format. Note that not all formats support a charge", default = None, type = float)
        sub_parser.add_argument("-M", "--multiplicity", help = "The multiplicity to set in the output format. Note that not all formats support a multiplicity", default = None, type = int)
        sub_parser.add_argument("--gen3D", help = "Whether to generate 3D coordinates (this will scramble existing atom coordinates). The default is yes, but only if it can be safely determined that the loaded coordinates are not already in 3D)", type = to_bool , default = True)
    
        return sub_parser
    
    def main(self):
        """
        Logic for our program.
        """
        # Panic if we don't have an input file.
        if self.args.input_file is None:
            raise Silico_exception("No input file specified")
        
        # Load the file we were given.
        parser = si_from_file(self.args.input_file, self.args.input_format, charge = self.args.charge, multiplicity = self.args.multiplicity, gen3D = self.args.gen3D)
        
        # If we weren't given an output format, try and guess one.
        if self.args.output_format is None:
            try:
                self.args.output_format = Openbabel_converter.type_from_file_name(self.args.output_file)
            except Silico_exception:
                # Couldn't guess a format, we'll just assume .si
                self.logger.info("No output format given and could not guess from file name; using .si as default ", exc_info = True)
                self.args.output_format = "si"
        
        # Convert and write.
        with Multi_file_wrapper(self.args.output_file, "wt") as outfile:
            outfile.write(parser.to_format(self.args.output_format))
            
    def load_interface(self, window):
        """
        Get the interface widget we'll use for display.
        
        :param window The parent window object.
        """
        return Convert_interface(window, self)
        