
# Silico imports.
from silico.interface.urwid.row_list.base import Row_list, Row_item, Row_widget
from silico.interface.urwid.file.base import File_selector
from silico.file.convert.main import Silico_input
from silico.interface.urwid.dialogue import Confirm_dialogue
from silico.interface.urwid.misc import IntEditZero

# General imports.
import pathlib
import urwid

class Coordinate_widget(Row_widget):
    """
    Widget for displaying Coordinate_item objects.
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.charge_edit = None
        self.mult_edit = None

    def load_inner(self):
        """
        Load the widget we'll use to display our main body.
        """
        self.charge_edit = IntEditZero(("body", "charge:"), self.row_item.coords.charge)
        self.mult_edit = IntEditZero(("body", "mult:"), self.row_item.coords.multiplicity) 
        
        return urwid.Columns([
            ("weight", 2, urwid.Text(self.row_item.coords.name)),
            urwid.Text(self.row_item.coords.formula),
            urwid.AttrMap(self.charge_edit, "editable"),
            urwid.AttrMap(self.mult_edit, "editable"),
        ], dividechars = 1)


class Coordinate_item(Row_item):
    """
    An item that appears in a Coordinate_list, represents a loaded coordinates file.
    """
    
    def __init__(self, path, row_list, movable=True):
        """
        Constructor for Coordinate_item objects.
        """
        super().__init__(str(path), row_list, movable=movable)
        
        # Load the coordinates from file.
        self.coords = Silico_input.from_file(path)
        
    def load_widget(self):
        """
        Load the widget we'll use for display.
        """
        return Coordinate_widget(self)


class Coordinate_list(Row_list):
    """
    A widget for displaying a list of loaded coordinates.
    """
    
    def __init__(self, top, rearrangeable=True):
        """
        """
        self.top = top
        
        # Widget for choosing files.
        self.file_selector = File_selector(pathlib.Path.cwd())
        
        Row_list.__init__(self, self.swap_to_selector, rearrangeable = rearrangeable)
        
    def swap_to_selector(self):
        """
        Swap to our file selector so the user can choose files.
        """
        self.top.swap_into_window(self.file_selector, cancel_callback = self.cancel, submit_callback = self.submit)
        
    def submit(self):
        """
        Add selected files to our list.
        """
        errors = []
        # Go through the files the user selected.
        for selected_node in self.file_selector.browser.selected_nodes:
            try:
                self.add_row(selected_node.get_value())
            
            except Exception as error:
                # We couldn't load the file for some reason, show an error and ignore.
                errors.append((selected_node.get_value(), error))
                
        # Clear the selected files.
        self.file_selector.browser.reset()
        
        # If we got any errors, show them to the user.
        if len(errors) > 0:
            error_text = "The following {} file(s) could not be loaded and will be ignored:\n\n".format(len(errors)) + "\n\n".join(["{}".format(error) for path, error in errors])
            
            self.top.swap(Confirm_dialogue("Error Loading Coordinates", error_text, self.top, submit_callback = 2))
            # Returning False stops us swapping back before we show the dialogue.
            return False
        
    def cancel(self):
        """
        Reset our file selector.
        """
        self.file_selector.browser.reset()
    
    def load_row(self, value):
        return Coordinate_item(value, self, self.rearrangeable)
    
    
    
    