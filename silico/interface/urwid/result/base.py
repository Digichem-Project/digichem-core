# General imports.
from pathlib import Path
import urwid

# Silico imports.
from silico.interface.urwid.main.program import Program_view
from silico.config.configurable.option import Option
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.section import Section
from silico.interface.urwid.result.list import Result_list
from silico.interface.urwid.setedit.base import Option_setedit


class Result_interface(Program_view):
    """
    Class for controlling interface to result program.
    """
    
    output = Option(help = "Path to the directory or file to write output to.", type = Path)
    ignore_missing = Option(help = "Ignore missing sections rather than stopping with an error.", type = bool)
    
    def __init__(self, window, program):
        """
        Constructor for calculation submitter objects.
        
        :param window: The program window we are being rendered in.
        :param program: A program object.
        """
        # Keep track of our individual widgets.
        self.file_list = Result_list(window.top, initial_results = program.results, subprocess_init = program.subprocess_init)
        
        # Set our options from our program object.
        self.output = program.args.output if program.args.output != "-" else None
        self.ignore_missing = program.args.ignore
        self.file_list.selector.alignment = program.config['alignment']
        self.file_list.selector.num_CPUs = program.args.num_CPUs
        
        super().__init__(window, program)
        
    def get_body(self):
        """
        Get the widget used to display the body of this program.
        
        :returns: An urwid widget to display.
        """
        return Tab_pile([
            Section(self.file_list, "Calculation Log Files"),
            ('pack', urwid.LineBox(Option_setedit.from_configurable_option(self.window.top, self, type(self).output).get_widget()))
        ])
        
    def setup(self):
        """
        Setup our program to ready it to be run.
        """
        