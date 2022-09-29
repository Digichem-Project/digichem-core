# Main urwid interface for generating reports.

# General imports.
import urwid

# Silico imports.
from silico.config.configurable.option import Option
from silico.interface.urwid.file.list import File_list
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.program import Program_view
from silico.interface.urwid.popup import Output_edit
from silico.interface.urwid.layout import Pane
from silico.interface.urwid.setedit.configurable import make_settings_page_from_configurable_option


class Report_interface(Program_view):
    """
    Class for controlling interface to reports generator.
    """
    
    name = Option(help = "Name of the molecule/system to use in the report.", type = str)
    report_type = Option(help = "The type of report to render.", type = str, choices = ["full", "atoms"], default = "full")
    report_style = Option(help = "The style of the report to render.", type = str, choices = ["journal", "traditional"], default = "journal")
    
    @property
    def additional_option_pages(self):
        """
        A dict of additional 'pages' of options to edit.
        """
        options = self.program.config.get_options()
        return dict([
            make_settings_page_from_configurable_option(self.window.top, self.program.config, self.program.config._configurable_options, options['absorption_spectrum']),
            make_settings_page_from_configurable_option(self.window.top, self.program.config, self.program.config._configurable_options, options['emission_spectrum']),
            make_settings_page_from_configurable_option(self.window.top, self.program.config, self.program.config._configurable_options, options['excited_states_diagram']),
            make_settings_page_from_configurable_option(self.window.top, self.program.config, self.program.config._configurable_options, options['IR_spectrum']),
            make_settings_page_from_configurable_option(self.window.top, self.program.config, self.program.config._configurable_options, options['rendered_image']),
            make_settings_page_from_configurable_option(self.window.top, self.program.config, self.program.config._configurable_options, options['orbital_diagram']),
            make_settings_page_from_configurable_option(self.window.top, self.program.config, self.program.config._configurable_options, options['report']),
            make_settings_page_from_configurable_option(self.window.top, self.program.config, self.program.config._configurable_options, options['skeletal_image'])
        ])
    
    
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
        self.report_style = program.args.style
        
        super().__init__(window, program)
        
    def setup(self):
        """
        Setup our program to ready it to be run.
        """
        self.program.args.output = self.output_widget.value
        self.program.args.name = self.name
        self.program.args.type = self.report_type
        self.program.args.log_files = self.file_list.get_values()
        self.program.args.style = self.report_style
    
    def get_body(self):
        """
        Get the widget used to display the body of this program.
        
        :returns: An urwid widget to display.
        """
        return Tab_pile([
            Pane(self.file_list, "Calculation Log Files"),
                ('pack', Pane(urwid.Pile([
                    urwid.Columns([
                        ('pack', urwid.Text("Output location:")),
                        urwid.AttrMap(self.output_widget, "editable"),
                    ], dividechars = 1)
                
                ]), "Output Options")
            )
        ])
        
        