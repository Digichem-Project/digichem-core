# General imports.
import urwid
import shlex

# Silico imports.
from silico.interface.urwid.main.program import Program_view
from silico.config.configurable.option import Option
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.section import Section
from silico.interface.urwid.result.list import Result_list
from silico.interface.urwid.edit.popup import Choices_edit, Output_edit
from silico.format.text import Text_summary_group_format
from silico.format.csv import CSV_property_group_format,\
    CSV_summary_group_format
from silico.format.table import Table_summary_group_format,\
    Property_table_group_format


class Result_interface(Program_view):
    """
    Class for controlling interface to result program.
    """
    
    stop_on_missing = Option(help = "Stop on missing properties rather than ignoring them.", type = bool)
    
    formats = {
        "text": Text_summary_group_format,
        "csv": CSV_summary_group_format,
        "csv-property": CSV_property_group_format,
        "table": Table_summary_group_format,
        "table-property": Property_table_group_format
    }
    
    def __init__(self, window, program):
        """
        Constructor for calculation submitter objects.
        
        :param window: The program window we are being rendered in.
        :param program: A program object.
        """
        # Keep track of our individual widgets.
        # A list of results to format.
        self.file_list = Result_list(window.top, initial_results = program.results, subprocess_init = program.subprocess_init)
        # The location to write the formatted results to.
        #self.output_widget = File_edit(window.top, program.args.output if program.args.output != "-" else None, "Output location")
        self.output_widget = Output_edit(window.top, program.args.output if program.args.output != "-" else None)
        
        # The output format.
        inital = list(self.formats.keys())[list(self.formats.values()).index(program.args.format)]
        self.format_widget = Choices_edit(window.top, [desc for desc, cls in self.formats.items()], inital, "Output format")
        # Specific filters selected by the user.
        self.filters_widget = urwid.Edit("", " ".join(program.args.filters))
        
        # Set our options from our program object.
        self.stop_on_missing = program.args.stop
        self.file_list.selector.alignment = program.config['alignment']
        self.file_list.selector.num_CPUs = program.args.num_CPUs
        
        super().__init__(window, program)
        
    def get_body(self):
        """
        Get the widget used to display the body of this program.
        
        :returns: An urwid widget to display.
        """
        return Tab_pile([
            Section(self.file_list, "Calculation Log Files"),('pack', 
                Section(urwid.Pile([
                    urwid.Columns([
                        ('pack', urwid.Text("Format:")),
                        urwid.AttrMap(self.format_widget, "editable"),
                    ], dividechars = 1),
                    
                    urwid.Columns([
                        ('pack', urwid.Text("Filters:")),
                        urwid.AttrMap(self.filters_widget, "editable"),
                    ], dividechars = 1),
                    
                    urwid.Columns([
                        ('pack', urwid.Text("Save location:")),
                        urwid.AttrMap(self.output_widget, "editable"),
                    ], dividechars = 1)
                
                ]), "Output Options")
            )
        ])
        
    def submit(self):
        self.program.logger.info("Writing results...")
        return super().submit()
    
    def post(self):
        self.program.logger.info("Done writing results")
        
    def setup(self):
        """
        Setup our program to ready it to be run.
        """
        self.program.args.output = self.output_widget.value
        self.program.args.stop = self.stop_on_missing
        self.program.results = self.file_list.get_values()
        self.program.args.format = self.formats[self.format_widget  .value]
        self.program.args.filters = shlex.split(self.filters_widget.get_edit_text()) if self.filters_widget.get_edit_text() != "" else []
        
        
        
        