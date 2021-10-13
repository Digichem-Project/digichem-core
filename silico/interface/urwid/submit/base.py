# Main urwid interface for submitting calculations.

# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.submit.coord import Coordinate_list
from silico.interface.urwid.base import Swappable, Section
from silico.interface.urwid.browser.base import Method_selector
from silico.interface.urwid.file.base import File_selector
import pathlib



class Calculation_submitter(Tab_pile, Swappable):
    """
    Class for controlling interface to calculation submission.
    """
    
    def __init__(self, methods):
        """
        Constructor for calculation submitter objects.
        
        :param methods: A configurable list of known methods (starting with the topmost destination).
        """
        Swappable.__init__(self)
        self.methods = methods
        
        # Keep track of our individual widgets.
        self.coordinate_list = Coordinate_list()
        #self.method_list = File_list()
        self.confirm = urwid.Button("Submit")
        
        # Make some browser windows to use.
        # We'll reuse these to keep focus.
        self.file_selector = File_selector(pathlib.Path.cwd())
        self.method_selector = Method_selector(self.methods)
                
        super().__init__([
            Section(self.coordinate_list, "Input Coordinates"),
            #Section(self.method_list, "Selected Methods"),
            (1, urwid.AttrMap(urwid.Filler(urwid.Padding(self.confirm, 'center', 11)), "goodButton"))
        ])
        
        
    def keypress(self, size, key):
        """
        Handler for keypress events.
        """
        if key in [' ', 'enter']:
            # Get the item in focus.
            if self.focus.base_widget == self.file_list:
                self.top.swap(self.file_selector)
                
            elif self.focus.base_widget == self.method_list:
                self.top.swap(self.method_selector)
                
            else:
                super().keypress(size, key)
            
        else:
            return super().keypress(size, key)