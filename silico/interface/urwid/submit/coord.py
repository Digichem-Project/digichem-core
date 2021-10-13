
# Silico imports.
import urwid


class Coordinate_item(urwid.Columns):
    """
    An item that appears in a Coordinate_list, represents a loaded coordinates file.
    """
    
    def __init__(self, coords):
        """
        
        :param coords: Input coords (a Silico_input object).
        """
        self.coords = coords

class Coordinate_list(urwid.listbox):
    """
    A widget for displaying a list of loaded coordinates.
    """
    
    def __init__(self):
        """
        """
        super().__init__([])