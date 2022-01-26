# Main urwid interface for generating reports.

# General imports.
import urwid

# Silico imports.
from silico.config.configurable.option import Option
from silico.interface.urwid.file.list import File_list
from silico.interface.urwid.section import Section
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.main.program import Program_view
from silico.interface.urwid.edit.popup import Output_edit


class Report_interface(Program_view):
    """
    Class for controlling interface to reports generator.
    """
    
    name = Option(help = "Name of the molecule/system to use in the report.", type = str)
    report_type = Option(help = "The type of report to render.", type = str, choices = ["full", "atoms"])
    
    
    def __init__(self, window, program):
        """
        Constructor for calculation submitter objects.
        
        :param window: The program window we are being rendered in.
        :param program: A program object.
        """
        # Keep track of our individual widgets.
        self.file_list = File_list(window.top, initial_files = program.args.log_files, can_choose_folders = True)
        # The location to write the formatted results to.
        self.output_widget = Output_edit(window.top, program.args.output, folder = False)
        
        # Set our options from our program object.
        self.name = program.args.name
        self.report_type = program.args.type
        
        super().__init__(window, program)
        
    def setup(self):
        """
        Setup our program to ready it to be run.
        """
        self.program.args.output = self.output_widget.value
        self.program.args.name = self.name
        self.program.args.type = self.report_type
        self.program.args.log_files = self.file_list.get_values()
    
    def get_body(self):
        """
        Get the widget used to display the body of this program.
        
        :returns: An urwid widget to display.
        """
        return Tab_pile([
            Section(self.file_list, "Calculation Log Files"),
                ('pack', Section(urwid.Pile([
                    urwid.Columns([
                        ('pack', urwid.Text("Output location:")),
                        urwid.AttrMap(self.output_widget, "editable"),
                    ], dividechars = 1)
                
                ]), "Output Options")
            )
        ])
        
        