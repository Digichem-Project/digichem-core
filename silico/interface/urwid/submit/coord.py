
# General imports.
import pathlib
import urwid

# Silico imports.
from silico.interface.urwid.row_list.base import Row_item, Row_widget,\
    Row_browser
from silico.interface.urwid.file.browser import Coord_selector
from silico.file.convert.main import Silico_input
from silico.interface.urwid.dialogue import Confirm_dialogue
from silico.interface.urwid.misc import IntEditZero


class Coordinate_widget(Row_widget):
    """
    Widget for displaying Coordinate_item objects.
    """
    
    def __init__(self, *args, **kwargs):
        self.charge_edit = None
        self.mult_edit = None
        super().__init__(*args, **kwargs)

    def load_inner(self):
        """
        Load the widget we'll use to display our main body.
        """
        self.charge_edit = IntEditZero(("body", "charge:"), self.row_item._value.charge)
        self.mult_edit = IntEditZero(("body", "mult:"), self.row_item._value.multiplicity) 
        
        return urwid.Columns([
            ("weight", 2, urwid.Text(self.row_item._value.auto_name)),
            urwid.Text(self.row_item._value.formula),
            urwid.AttrMap(self.charge_edit, "editable"),
            urwid.AttrMap(self.mult_edit, "editable"),
        ], dividechars = 1)


class Coordinate_item(Row_item):
    """
    An item that appears in a Coordinate_list, represents a loaded coordinates file.
    """
    
    def __init__(self, coord, row_list, movable = True,):
        """
        Constructor for Coordinate_item objects.
        
        :param coord: A Silico_input object.
        :param row_list: Our parent row_list object.
        :param movable: Whether we should show rearranging buttons.
        """
        super().__init__(coord, row_list, movable = movable)
        
    def load_widget(self):
        """
        Load the widget we'll use for display.
        """
        return Coordinate_widget(self)
    
    def get_value(self):
        """
        Get the current value of this coordinate item object.
        """
        self._value.charge = self.get_widget().charge_edit.value()
        self._value.multiplicity = self.get_widget().mult_edit.value()
        return self._value


class Coordinate_list(Row_browser):
    """
    A widget for displaying a list of loaded coordinates.
    """
    
    def __init__(self, top, start_dir = None, rearrangeable = True, initial_coords = None, initial_charge = None, initial_mult = None, gen3D = True):
        """
        Constructor for Coordinate_list objects.
        
        :param top: A top-most widget to use for display.
        :param start_dir: The directory to show when first opening the file selector.
        :param rearrangeable: Whether the list of loaded files can be rearranged by the user.
        :param initial_coords: A list of coordinates to show by default.
        :param initial_charge: An optional integer charge to set for initial_files.
        :param initial_mult: An optional integer multiplicity to set for initial_files.
        """
        start_dir = start_dir if start_dir is not None else pathlib.Path.cwd()
        self.top = top
        
        # Widget for choosing files, stored at self.selector
        file_selector = Coord_selector(start_dir, title = "Select Input Coordinates")
        file_selector.charge = initial_charge
        file_selector.multiplicity = initial_mult
        file_selector.generate_3D = gen3D
        
        # Add our starting files.
        initial_coords = [] if initial_coords is None else initial_coords
        
        super().__init__(file_selector, top = top, rearrangeable = rearrangeable, initial = initial_coords)
        
    def submit(self):
        """
        Add selected files to our list.
        """
        errors = []
        # Go through the files the user selected.
        for selected_node in self.selector.browser.selected_nodes:
            try:
                self.add_row(self.value_from_node(selected_node))
            
            except Exception as error:
                # We couldn't load the file for some reason, show an error and ignore.
                errors.append((selected_node.get_value(), error))
                
        # Clear the selected files.
        self.selector.browser.reset()
        
        # If we got any errors, show them to the user.
        if len(errors) > 0:
            error_text = "The following {} file(s) could not be loaded and will be ignored:\n\n".format(len(errors)) + "\n\n".join(["{}".format(error) for path, error in errors])
            self.top.popup(Confirm_dialogue("Error Loading Coordinates", error_text, self.top, error = True))
        
    def value_from_node(self, node):
        """
        Get a value to add from a selected node.
        """
        return Silico_input.from_file(node.get_value(), charge = self.selector.charge, multiplicity = self.selector.multiplicity, gen3D = self.selector.generate_3D)
        
    def load_row(self, value):
        return Coordinate_item(value, self)
    
