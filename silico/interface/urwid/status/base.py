# General imports.
import urwid

# Silico imports
from silico.interface.urwid.swap.swappable import Swappable
from silico import logging


class Status_interface(Swappable):
    """
    Class for controlling interface to queue status information.
    """
    
    def __init__(self, window, program):
        """
        Constructor for calculation submitter objects.
        
        :param window: The program window we are being rendered in.
        :param program: A program object.
        """
        self.window = window
        self.program = program
        
        # Keep track of our individual widgets.
        #self.time_widget = urwid.AttrMap(urwid.IntEdit(("body", "Refresh rate /s: "), 5), "editable")
        self.status_widget = urwid.Text("")
        
        # Set our options from our program object.
        
        super().__init__(self.window.top, self.get_body(), title = program.name, border = False)
        print("Hello")
        logging.get_logger().info("Hello2")
        #self.program.main()
        
    
    def get_body(self):
        """
        Get the widget used to display the body of this program.
        
        :returns: An urwid widget to display.
        """
        return urwid.Pile([
            #('pack', self.time_widget),
            urwid.Filler(self.status_widget)
        ])