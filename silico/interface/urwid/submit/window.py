# General imports.
import urwid


class Submitter_window(urwid.ListBox):
    """
    A widget that holds data visible in the calculation submitter.
    """
    
    def selectable(self):
        """
        We are secretly a button that will open another window.
        """
        return True
    
    def keypress(self, size, key):
        """
        Dummy keypress function.
        """
        return key