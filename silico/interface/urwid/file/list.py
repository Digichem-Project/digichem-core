
# General imports.
import pathlib
import urwid

# Silico imports.
from silico.interface.urwid.row_list import Row_item, Row_widget,\
    Row_browser
from silico.interface.urwid.file.browser import File_selector


class File_widget(Row_widget):
    """
    Widget for displaying Coordinate_item objects.
    """
    
    def load_inner(self):
        """
        Load the widget we'll use to display our main body.
        """
        return urwid.Text(str(self.row_item.get_value()))


class File_item(Row_item):
    """
    An item that appears in a File_list, represents a path to a file.
    """
    
    def __init__(self, path, row_list, movable = True):
        """
        Constructor for File_item objects.
        
        :param path: The path to the file we represent (pathlib.Path).
        :param row_list: Our parent row_list object.
        :param movable: Whether we should show rearranging buttons.
        """
        super().__init__(path, row_list, movable = movable)
        
    def load_widget(self):
        """
        Load the widget we'll use for display.
        """
        return File_widget(self)


class File_list(Row_browser):
    """
    A widget for displaying a list of paths to files.
    """
    
    def __init__(self, top, start_dir = None, rearrangeable = True, initial_files = None, can_choose_folders = False):
        """
        Constructor for Coordinate_list objects.
        
        :param top: A top-most widget to use for display.
        :param start_dir: The directory to show when first opening the file selector.
        :param rearrangeable: Whether the list of loaded files can be rearranged by the user.
        :param initial_files: A list of files to show by default.
        """
        start_dir = start_dir if start_dir is not None else pathlib.Path.cwd()
        
        # Widget for choosing files.
        file_selector = File_selector(top, start_dir, can_choose_folders = can_choose_folders)
        
        # Add our starting files.
        initial_files = [] if initial_files is None else initial_files
        
        super().__init__(file_selector, top = top, rearrangeable = rearrangeable, initial = initial_files)
    
    def load_row(self, value):
        return File_item(value, self, self.rearrangeable)
    
    
    
    