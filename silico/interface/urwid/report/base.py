# Main urwid interface for generating reports.

# General imports.
from pathlib import Path

# Silico imports.
from silico.config.configurable.option import Option
from silico.interface.urwid.file.list import File_list
from silico.interface.urwid.section import Section
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.main.base import Program_view


class Report_generator(Program_view):
    """
    Class for controlling interface to reports generator.
    """
    
    output = Option(help = "Base directory to which output will be written.", type = Path,)
    name = Option(help = "Name of the molecule/system to use in the report.", type = str)
    report_type = Option(help = "The type of report to render.", type = str, choices = ["full", "atoms"])
    render_style = Option(help = "The visual style with which to render orbital and molecular images.", type = str, choices=['pastel', 'light-pastel', 'dark-pastel','sharp', 'gaussian', 'vesta'])
    render_images = Option("render", help = "Whether to render new images.", type = bool)
    overwrite = Option(help = "Whether to overwrite existing image files. Useful for forcing re-rendering, for example if the render_style has been changed.", type = bool)
    auto_crop = Option(help = "Whether to automatically crop rendered images to remove excess whitespace. If false, rendered images will have the same resolution, but molecules may only occupy a small portion of the true image", type = bool)
    
    
    def __init__(self, window, program):
        """
        Constructor for calculation submitter objects.
        
        :param window: The program window we are being rendered in.
        :param program: A program object.
        """
        # Keep track of our individual widgets.
        self.file_list = File_list(window.top, initial_files = program.args.log_files)
        
        super.__init__(window, program)
    
    def get_body(self):
        """
        Get the widget used to display the body of this program.
        
        :returns: An urwid widget to display.
        """
        return Tab_pile([
            Section(self.file_list, "Calculation Files"),
        ])