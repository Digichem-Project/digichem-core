
# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.row_list import Row_item, Row_widget,\
    Row_browser
from silico.input import si_from_file
from silico.interface.urwid.misc import IntEditZero
import silico.log
from silico.interface.urwid.coord.browser import Coord_selector


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
            ("weight", 2, urwid.Text(self.row_item._value.implicit_name)),
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
        
        :param coord: A Silico_coords object.
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
        top = top
        
        # Widget for choosing files, stored at self.selector
        file_selector = Coord_selector(top, start_dir)
        file_selector.charge = initial_charge
        file_selector.multiplicity = initial_mult
        file_selector.generate_3D = gen3D
        
        # Add our starting files.
        initial_coords = [] if initial_coords is None else initial_coords
        
        super().__init__(file_selector, top = top, rearrangeable = rearrangeable, initial = initial_coords)
        
    def add_to_list(self):
        """
        Add selected files to our list.
        """
        # Get the logger item we'll use for communication.
        logger = silico.log.get_logger()
        logger.info("Loading coordinates")
        
        # Go through the files the user selected.
        for path in self.selector.selected:
            try:
                self.add_row(self.value_from_selected(path))
             
            except Exception:
                # We couldn't load the file for some reason, show an error and ignore.
                logger.error("Failed to load coordinate file '{}':".format(path), exc_info = True)
                 
        # Clear the selected files.
        self.selector.reset()
         
        logger.info("Finished loading coordinates")
        
    def value_from_selected(self, selected_value):
        """
        Get a value to add from a value selected by the user.
        """
        return si_from_file(selected_value, charge = self.selector.charge, multiplicity = self.selector.multiplicity, gen3D = self.selector.generate_3D)
        
    def load_row(self, value):
        return Coordinate_item(value, self)
    
