
# General imports.
import pathlib
import urwid

# Silico imports.
from silico.interface.urwid.row_list.base import Row_list, Row_item, Row_widget
from silico.interface.urwid.file.base import Coord_selector
from silico.file.convert.main import Silico_input
from silico.interface.urwid.dialogue import Confirm_dialogue
from silico.interface.urwid.misc import IntEditZero, FloatEditZero


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
    
    def __init__(self, path, row_list, movable = True, initial_charge = None, initial_mult = None, gen3D = True):
        """
        Constructor for Coordinate_item objects.
        
        :param path: The path to the coordinate file we represent.
        :param row_list: Our parent row_list object.
        :param movable: Whether we should show rearranging buttons.
        :param initial_charge: The initial charge.
        :param initial_mult: The initial multiplicity.
        :param gen3D: Whether to try and convert 2D formats to 3D coords.
        """
        super().__init__(path, row_list, movable = movable)
        
        # Load the coordinates from file.
        self.coords = Silico_input.from_file(path, charge = initial_charge, multiplicity = initial_mult, gen3D = gen3D)
        
    def load_widget(self):
        """
        Load the widget we'll use for display.
        """
        return Coordinate_widget(self)


class Coordinate_list(Row_list):
    """
    A widget for displaying a list of loaded coordinates.
    """
    
    def __init__(self, top, start_dir = None, rearrangeable = True, initial_files = None, initial_charge = None, initial_mult = None):
        """
        Constructor for Coordinate_list objects.
        
        :param top: A top-most widget to use for display.
        :param start_dir: The directory to show when first opening the file selector.
        :param rearrangeable: Whether the list of loaded files can be rearranged by the user.
        :param initial_files: A list of files to show by default.
        :param initial_charge: An optional integer charge to set for initial_files.
        :param initial_mult: An optional integer multiplicity to set for initial_files.
        """
        start_dir = start_dir if start_dir is not None else pathlib.Path.cwd()
        self.top = top
        
        # Widget for choosing files.
        self.file_selector = Coord_selector(start_dir, title = "Select Input Coordinates")
        self.file_selector.charge = initial_charge
        self.file_selector.multiplicity = initial_mult
        
        # Add our starting files.
        initial_files = [] if initial_files is None else initial_files
        
        super().__init__(self.swap_to_selector, rearrangeable = rearrangeable, initial = initial_files)
        
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
        return Coordinate_item(value, self, self.rearrangeable, gen3D = self.file_selector.generate_3D, initial_charge = self.file_selector.charge, initial_mult = self.file_selector.multiplicity)
    
    
    
    