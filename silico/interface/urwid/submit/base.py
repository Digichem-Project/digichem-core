# Main urwid interface for submitting calculations.

# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.program import Program_view
from silico.interface.urwid.coord.list import Coordinate_list
from silico.interface.urwid.method.list import Method_list
from silico.interface.urwid.edit.popup import Output_edit
from silico.interface.urwid.layout import Pane
from silico.config.configurable.option import Option


class Submit_interface(Program_view):
    """
    Class for controlling interface to calculation submission.
    """
    
    prepare_only = Option(help = "Whether to only perform setup for the calculation without actually performing it", type = bool, default = False)
    
    def __init__(self, window, program):
        """
        Constructor for calculation submitter objects.
        
        :param window: The program window we are being rendered in.
        :param program: A program object.
        """
        # Keep track of our individual widgets.
        self.coordinate_list = Coordinate_list(window.top, initial_coords = program.coords, initial_charge = program.args.charge, initial_mult = program.args.multiplicity, gen3D = program.args.gen3D)
        self.method_list = Method_list(window.top, method_library = program.config, initial_methods = program.methods)
        # The location to write the formatted results to.
        self.output_widget = Output_edit(window.top, program.args.output, folder = True)
        self.prepare_only = program.args.prepare_only
        
        super().__init__(window, program)
        
    def get_body(self):
        """
        Get the widget used to display the body of this program.
        
        :returns: An urwid widget to display.
        """
        return Tab_pile([
            Pane(self.coordinate_list, "Input Coordinates"),
            Pane(self.method_list, "Calculation Methods"),
                ('pack', Pane(urwid.Pile([
                    urwid.Columns([
                        ('pack', urwid.Text("Output location:")),
                        urwid.AttrMap(self.output_widget, "editable"),
                    ], dividechars = 1)
                
                ]), "Output Options")
            )
        ])
        
    def setup(self):
        """
        Setup our program to ready it to be run.
        """
        self.program.coords = self.coordinate_list.get_values()
        self.program.methods = self.method_list.get_values()
        self.program.args.output = self.output_widget.value
        self.program.args.prepare_only = self.prepare_only
        