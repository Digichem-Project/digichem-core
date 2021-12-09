# General imports.
from pathlib import Path
import itertools

# Silico imports.
from silico.program.base import Program
from silico.exception.base import Silico_exception
from silico.parser import parse_calculations
from silico.report.main.pdf import PDF_report
from silico.interface.urwid.report.base import Report_interface


class Report_program(Program):
    """
    The Silico report generation program.
    """
    
    name = "Calculation Report Generator"
    command = "report"
    description = "generate PDF reports from finished calculations"
    usage = "%(prog)s [options] file.log [file2.log] [-O report.pdf]"
    help = "Write reports"
    
    @classmethod
    def arguments(self, sub_parsers_object):
        """
        Add this program's arguments to an argparser subparsers object.
        
        :param sub_parsers_object: An argparser subparsers object to add to (as returned by add_subparsers().
        :returns: The argparser object (which supports add_argument()).
        """
        sub_parser = super().arguments(sub_parsers_object)
        
        sub_parser.add_argument("log_files", help = "a (number of) calculation result file(s) (.log) to extract results from", nargs = "*")
                            
        output_group = sub_parser.add_mutually_exclusive_group()
        output_group.add_argument("--pdf_file", help = "a filename/path to a pdf file to write to (this is an alternative to the 'output' option). Other output files will placed in the same directory as the 'pdf_file'", nargs = "?", default = None)
        output_group.add_argument("-O", "--output", help = "a filename/path to write to. If this filename ends in a .pdf extension, this is taken as the filename of the pdf file to write to. Otherwise, output is assumed to be the name of a directory and the pdf file is written inside this directory", nargs = "?", type = Path, default = None)
        
        sub_parser.add_argument("--name", help = "name of the molecule/system to use in the report", default = None)
        sub_parser.add_argument("--type", help = "the type of report to make", choices = ["full", "atoms"], default = "full")
        
        aux_input_group = sub_parser.add_argument_group("auxiliary input", "options for specifying additional input files. These options are all optional, but if given the ordering should match that of the given log files")
        aux_input_group.add_argument("--chk", help = "a (number of) Gaussian chk file(s) that will be used to generate all image files required. Note that this option requires Gaussian to be installed and formchk & cubegen to be in your path", nargs = "*", default = [])
        aux_input_group.add_argument("--fchk", help = "a (number of) Gaussian fchk file(s) that will be used to generate all image files required. Note that this option requires Gaussian to be installed and cubegen to be in your path", nargs = "*", default = [])
        
        image_group = sub_parser.add_argument_group("image options", "advanced options for altering the appearance of images in the report")
        image_group.add_argument("--render_style", help = "change the style used to render molecules and orbitals", choices=['pastel', 'light-pastel', 'dark-pastel','sharp', 'gaussian', 'vesta'], default = None)
        image_group.add_argument("-d", "--dont_create_new_images", help = "don't attempt to create any new image files. If existing image files can be found, they will still be included in the report", action = "store_true", default = None)
        image_group.add_argument("-o", "--overwrite_existing_images", help = "force overwriting and re-rendering of all image files", action = "store_true", default = None)
        image_group.add_argument("--dont_auto_crop", help = "disable automatic cropping of excess whitespace around the border of images. If given, rendered images will have the same resolution, but molecules may only occupy a small portion of the true image", action = "store_true", default = None)
    
        return sub_parser
    
    def process_arguments(self, args):
        """
        Process command line arguments and set defaults.
        """
        
    @classmethod
    def arg_to_config(self, args, config):
        """
        A class method that will be called to add command line arguments from 'args' to the configuration object 'config'.
        
        :param args: An argparser object.
        :param config: A Silico config object to add to.
        """
        # Add our report specific command line arguments to our config dictionary.
        if args.render_style is not None:
            config['molecule_image']['rendering_style'] = args.render_style
        if args.overwrite_existing_images is not None:
                config['molecule_image']['use_existing'] = not args.overwrite_existing_images
        if args.dont_auto_crop is not None:
                config['molecule_image']['auto_crop'] = not args.dont_auto_crop
        
        if args.dont_create_new_images is not None:
            config['image']['dont_modify'] = args.dont_create_new_images
        
        
    def load_interface(self, window):
        """
        Get the interface widget we'll use for display.
        
        :param window The parent window object.
        """
        return Report_interface(window, self)
        
    
    def main(self):
        """
        Logic for our program.
        """
        # Get upset if we have no log files.
        if len(self.args.log_files) == 0:
            raise Silico_exception("No log files specified")
        
        try:
            # Transpose our lists of aux inputs.
            aux_files = [{"chk_file":chk_file, "fchk_file":fchk_file} for chk_file, fchk_file in list(itertools.zip_longest(self.args.chk, self.args.fchk))]
            
            # First, load results.
            result = parse_calculations(*self.args.log_files, aux_files = aux_files)
            
            # Then get a report.
            report = PDF_report(result, options = self.config)
            
        except Exception as e:
            raise Silico_exception("Failed to load results")
        
        # Forcibly set a name if we've been given one.
        if self.args.name is not None:
            result.metadata.name = self.args.name
        
        log_file = self.args.log_files[0]
        if log_file is not None:
            input_name = Path(log_file)
        elif len(self.args.calculation_files) > 0:
            input_name = Path(self.args.calculation_files[0])
        else:
            input_name = "Report"
        
        # The filename (excluding parent directories) we'll write to if the user doesn't give us one.
        if self.args.type == "full":
            #default_base_name = input_name.with_suffix(".pdf").name
            default_base_name = result.metadata.name + ".pdf"
        else:
            #default_base_name = input_name.with_suffix(".atoms.pdf").name
            default_base_name = result.metadata.name + ".atoms.pdf"
        
        # Decide on our output file.
        if self.args.pdf_file is not None:
            # We've been given a specific file to write to. No problems.
            self.args.pdf_file = Path(self.args.pdf_file)
        elif self.args.output is not None:
            # We've been given the more general 'output' option, we need to decide whether this aught to be a file or directory.
            if self.args.output.suffix == ".pdf":
                # It's a file! Set pdf_file appropriately.
                self.args.pdf_file = self.args.output
            else:
                # It's a directory. Decide on a fitting name from our input file.
                self.args.pdf_file = Path(self.args.output, default_base_name)
        else:
            # No output was given, we'll use a default location (in the same place as the main input file).
            self.args.pdf_file = Path(input_name.parent, "report", default_base_name)
            
        # Now make (compile?) and save our report.
        self.logger.info("Generating report '{}'...".format(self.args.pdf_file))
        try:
            report.write(self.args.pdf_file, report_type = self.args.type)
        except Exception as e:
            raise Silico_exception("Failed to generate report '{}'".format(self.args.pdf_file)) from e
        
        # Done.
        self.logger.info("Done generating report '{}'".format(self.args.pdf_file))
        