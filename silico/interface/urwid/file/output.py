# Widgets for selecting output/save locations.

# General imports.
from pathlib import Path
import urwid

# Silico imports.
from silico.interface.urwid.file.browser import File_browser, Selector_mixin
from silico.interface.urwid.swap.swappable import Swappable
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.layout import Pane


class Output_tree_list_box(File_browser):
    """
    A modified tree list box that allows selecting output locations.
    """
     
    def __init__(self, selector, starting_dir, show_hidden = False, refresh = True, folder = False):
        """
        :param selector: The parent selector object we are a part of.
        :param starting_dir: Path to the directory to show open by default.
        :param show_hidden: Whether to show hidden files.
        :param refresh: Whether to refresh directories on close.
        :param folder: Whether to allow selection of folders or files.s
        """
        self.selector = selector
        self.folder = folder
        
        # Keep track of our file and directory separately.
        self.selected_file = None
        self.selected_directory = None
        
        super().__init__(starting_dir, show_hidden = show_hidden, refresh = refresh,  can_choose_folders = True, can_choose_files = not folder, can_choose_multiple = False)
    
    def select(self, focus_node, focus_widget):
        """
        Select (or deselect) the node currently in focus.
        """
        super().select(focus_node, focus_widget)
        selected = self.selected[0] if len(self.selected) > 0 else None
        
        if selected is not None:
            # Get the 'file' name part of what we selected
            if self.folder:
                # We are selecting a folder; our 'name' is the last part of the path (which should be a folder).
                self.selected_file = selected.name
                self.selected_directory = selected.parents[0]
                
            else:
                # We are selecting a file, whether we have a file name depends on whether a directory or file was chosen (we allow selecting both).
                if hasattr(self.selected_nodes[0], "has_children"):
                    # A directory, don't set file name.
                    self.selected_file = None
                    self.selected_directory = selected
                    
                else:
                    # A file, set file name.
                    self.selected_file = selected.name
                    self.selected_directory = selected.parents[0]
            
        else:
            self.selected_file = None
            self.selected_directory = None
            
        # If we have a file_name set, update.
        if self.selected_file is not None:
            self.selector.file_name_widget.set_edit_text(str(self.selected_file))
            
        self.selector.update_from_browser()
        

class Output_selector(Swappable, Selector_mixin):
    """
    A tree list box widget for selecting a location to save something.
    """
    
    def __init__(self, top, default = None, folder = False, default_file_name = ""):
        """
        Constructor for Output_selector objects.
        
        :param top: Top-most widget being used for display.
        :param default: The save location to show by default.
        :param title: A title to display around this selector.
        :param folder: Whether the output location should be a directory.
        """
        # Work out our starting directory for our browser.
        if default is None:
            # No default.
            starting_dir = None
            # Use a default file name if we have one.
            file_name = default_file_name
        
        else:
            default = Path(default)
            # We'll try and resolve our path to make it real.
            default = default.resolve()
            
            if folder:
                # Our output location is a folder.
                # We can safely assume the last part of our path is our file name.
                file_name = default.name
                starting_dir = default
                
            else:
                # Output location is a file, but default could be either a file or a directory.
                # Try and work out if default is a dir or not.
                # If the file exists this is easy, we can just check.
                # Otherwise, if the file name has an extension then we'll assume its a file, otherwise a dir.
                if default.is_dir() or len(default.suffixes) == 0:
                    # It is a directory, don't set file_name.
                    file_name = ""
                    starting_dir = default
                    
                else:
                    # Not a directory (or not real).
                    file_name = default.name
                    starting_dir = default.parents[0]
                
        # Save our default.
        self.default = default
        self.default_file_name = default_file_name
        self._updating = False
                
        # Keep track of our sub widgets.
        self.browser = Output_tree_list_box(self, starting_dir = starting_dir, folder = folder)
        self.file_name_widget = urwid.Edit(("body", "File name: "), file_name)
        self.location_widget = urwid.Edit(("body", "Location: "), str(default) if default is not None else "")
        
        body = Tab_pile([
            Pane(self.browser, "Choose file or folder"),
            ('pack', Pane(
                urwid.Pile([
                    urwid.AttrMap(self.file_name_widget, "editable"),
                    urwid.AttrMap(self.location_widget, "editable")
                ]), 
                "Output location"
            ))
        ])
        
        super().__init__(top, body)
        
        # Connect signals.
        urwid.connect_signal(self.file_name_widget, 'postchange', self.update_from_file_name_widget)
        urwid.connect_signal(self.location_widget, 'postchange', self.update_from_location_widget)
        
    def update_location_widget(self):
        """
        Update the location widget.
        """
        selected_path = self.browser.selected_directory
        
        # If nothing selected, set location as nothing.
        if selected_path is None:
            self.location_widget.set_edit_text("")
            
        else:
            # Add our file name to path.
            selected_path = Path(selected_path, self.file_name_widget.get_edit_text())
            self.location_widget.set_edit_text(str(selected_path))
            
    def update_from_browser(self):
        """
        Method called when the contents of the file browser has changed.
        """
        if not self._updating:
            try:
                self._updating = True
                self.update_location_widget()
            
            finally:
                self._updating = False
            
    def update_from_file_name_widget(self, *args):
        """
        Method called when the contents of the file name widget has changed.
        """
        if not self._updating:
            try:
                self._updating = True
                self.update_location_widget()
            
            finally:
                self._updating = False
            
    def update_from_location_widget(self, *args):
        """
        Method called when the contents of the location widget has changed.
        """
        if not self._updating:
            try:
                self._updating = True
                # Clear the browser.
                self.browser.reset()
                
                # Take the name and add to our file name widget.
                file_name = Path(self.location_widget.get_edit_text()).name
                
                self.file_name_widget.set_edit_text(file_name)
            
            finally:
                self._updating = False
        
    @property
    def value(self):
        """
        The selected output location.
        """
        return Path(self.location_widget.get_edit_text()) if self.location_widget.get_edit_text() != "" else None
        
    def on_settings_change(self, confirm):
        """
        A method that will be called when settings have been changed.
        """
        Selector_mixin.on_settings_change(self, confirm)
        
        