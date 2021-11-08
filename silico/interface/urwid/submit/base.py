# Main urwid interface for submitting calculations.

# General imports.

# Silico imports.
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.submit.coord import Coordinate_list
from silico.interface.urwid.section import Section
from silico.interface.urwid.submit.method import Method_list
from silico.config.configurable.option import Option
from silico.interface.urwid.main.base import Program_view
from silico.interface.urwid.dialogue import Confirm_dialogue


class Calculation_submitter(Program_view):
    """
    Class for controlling interface to calculation submission.
    """
    
    output = Option(help = "Base directory to which output will be written.", type = str, default = "./")
    
    def __init__(self, window, program):
        """
        Constructor for calculation submitter objects.
        
        :param window: The program window we are being rendered in.
        :param program: A program object.
        """
        # Keep track of our individual widgets.
        self.coordinate_list = Coordinate_list(window.top, initial_coords = program.coords, initial_charge = program.args.charge, initial_mult = program.args.multiplicity, gen3D = program.args.gen3D)
        self.method_list = Method_list(window.top, program.config.methods, initial_methods = program.methods)
        
        super().__init__(window, program)
        
    def get_body(self):
        """
        Get the widget used to display the body of this program.
        
        :returns: An urwid widget to display.
        """
        return Tab_pile([
            Section(self.coordinate_list, "Input Coordinates"),
            Section(self.method_list, "Calculation Methods"),
        ])
        
    def setup(self):
        """
        Setup our program to ready it to be run.
        """
        self.program.coords = self.coordinate_list.get_values()
        self.program.methods = self.method_list.get_values()
        self.program.output = self.output
        
    def post(self, retval):
        """
        Method called once our main program has finished running.
        """
        self.window.top.popup(Confirm_dialogue("Submission Complete", "Successfully submitted {} calculations".format(retval), self.window.top, submit_callback = 1))
        