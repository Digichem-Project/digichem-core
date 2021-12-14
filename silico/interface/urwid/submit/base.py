# Main urwid interface for submitting calculations.

# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.section import Section
from silico.interface.urwid.main.program import Program_view
from silico.interface.urwid.coord.list import Coordinate_list
from silico.interface.urwid.method.list import Method_list
from silico.interface.urwid.edit.popup import File_edit


class Submit_interface(Program_view):
    """
    Class for controlling interface to calculation submission.
    """
    
    def __init__(self, window, program):
        """
        Constructor for calculation submitter objects.
        
        :param window: The program window we are being rendered in.
        :param program: A program object.
        """
        # Keep track of our individual widgets.
        self.coordinate_list = Coordinate_list(window.top, initial_coords = program.coords, initial_charge = program.args.charge, initial_mult = program.args.multiplicity, gen3D = program.args.gen3D)
        self.method_list = Method_list(window.top, program.config.methods, initial_methods = program.methods)
        # The location to write the formatted results to.
        self.output_widget = File_edit(window.top, program.args.output, "Output location")
        
        super().__init__(window, program)
        
    def get_body(self):
        """
        Get the widget used to display the body of this program.
        
        :returns: An urwid widget to display.
        """
        return Tab_pile([
            Section(self.coordinate_list, "Input Coordinates"),
            Section(self.method_list, "Calculation Methods"),
                ('pack', Section(urwid.Pile([
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
        