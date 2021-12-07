# General imports.
import logging

# Silico imports.

import silico
from silico.interface.urwid.swap.top import Swappable


class Program_view(Swappable):
    """
    A widget for running a sub-program.
    """
    
    def __init__(self, window, program):
        """
        Constructor for Program_view objects.
        
        :param window: A Silico_window widget used for display.
        :param program: A program object to run.
        """
        self.window = window
        self.program = program
        super().__init__(self.window.top, self.get_body(), title = program.name, border = False)
        
    def get_body(self):
        """
        Get the widget used to display the body of this program.
        
        :returns: An urwid widget to display.
        """
        raise NotImplementedError("Implement in subclass")
    
    def setup(self):
        """
        Setup our program to ready it to be run.
        """
        # This default implementation does nothing.
    
    def submit(self):
        """
        Submit this program widget, running the program it wraps.
        """
        try:
            self.setup()
            retval = self.program.main()
            self.post(retval)
            
        except Exception:
            logger = logging.getLogger(silico.logger_name)
            logger.exception("Sub-program {} stopped with error".format(self.program.command))
            return False
        
        except KeyboardInterrupt:
            # The user wanted to stop.
            logger = logging.getLogger(silico.logger_name)
            logger.error("Sub-program {} interrupted by user (ctrl-c)".format(self.program.command))
            return False
        
    def post(self, retval):
        """
        Method called once our main program has finished running.
        """
        # This default implementation does nothing.


