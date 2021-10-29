# Main urwid interface for submitting calculations.

# General imports.

# Silico imports.
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.submit.coord import Coordinate_list
from silico.interface.urwid.section import Section
from silico.interface.urwid.submit.method import Method_list
from silico.interface.urwid.top import View
from silico.config.configurable.option import Option


class Calculation_submitter(View):
    """
    Class for controlling interface to calculation submission.
    """
    
    output = Option(help = "Base directory to which output will be written.", type = str, default = "./")
    
    def __init__(self, top, methods, initial_files = None, initial_charge = None, initial_mult = None, initial_methods = None):
        """
        Constructor for calculation submitter objects.
        
        :param top: The topmost widget used for display.
        :param methods: A configurable list of known methods (starting with the topmost destination).
        :param initial_files: A list of file paths to initially populate with.
        :param initial_methods: A list of methods (tuples of destination, program, calculation) to initially populate with.
        """
        self.top = top
        self.methods = methods
        
        # Keep track of our individual widgets.
        self.coordinate_list = Coordinate_list(top, initial_charge = initial_charge, initial_mult = initial_mult, initial_files = initial_files)
        self.method_list = Method_list(top, self.methods, initial_methods = initial_methods)
        
        self.pile = Tab_pile([
            Section(self.coordinate_list, "Input Coordinates"),
            Section(self.method_list, "Calculation Methods"),
        ])
        
        super().__init__(self.pile, "Calculation Submission", border = False)
        