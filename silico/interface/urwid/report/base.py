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
    
    output = Option(help = "Base directory to which output will be written.", type = Path, default = "./")
    name = Option(help = "Name of the molecule/system to use in the report.", type = str, default = None)
    report_type = Option(help = "The type of report to render.", type = str, default = "full", choices = ["full", "atoms"])
    render_style = Option(help = "The visual style with which to render orbital and molecular images.", type = str, choices=['pastel', 'light-pastel', 'dark-pastel','sharp', 'gaussian', 'vesta'], default = "pastel")
    render_images = Option("render", help = "Whether to render new images.", type = bool, default = True)
    overwrite = Option(help = "Whether to overwrite existing image files. Useful for forcing re-rendering, for example if the render_style has been changed.", type = bool, default = False)
    auto_crop = Option(help = "Whether to automatically crop rendered images to remove excess whitespace. If false, rendered images will have the same resolution, but molecules may only occupy a small portion of the true image", type = bool, default = True)
    
    
    def __init__(self, top, initial_files = None, ):
        """
        Constructor for calculation submitter objects.
        
        :param top: The topmost widget used for display.
        :param methods: A configurable list of known methods (starting with the topmost destination).
        :param initial_files: A list of file paths to initially populate with.
        :param initial_methods: A list of methods (tuples of destination, program, calculation) to initially populate with.
        """
        self.top = top
        
        # Keep track of our individual widgets.
        self.file_list = File_list(top, initial_files = initial_files)
        
        self.pile = Tab_pile([
            Section(self.file_list, "Calculation Files"),
        ])
        
        super().__init__(self.pile, "Calculation Submission", border = False)
        