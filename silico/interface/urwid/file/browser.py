# Base classes for file browser widgets.
# This code is adapted from the example file browser included with urwid, originally written by Rob Lanphier.

# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.tree.base import Flaggable_tree_list_box,\
    Flaggable_tree_walker
from silico.config.configurable.option import Option
from silico.interface.urwid.file.node import Directory_node
from silico.interface.urwid.tree.selector import Enhanced_tree_selector
import pathlib


class File_browser(Flaggable_tree_list_box):
    """
    Inner widget for displaying a tree of files and directories.
    """
    
    def __init__(self, starting_dir, show_hidden = False, can_choose_folders = False, can_choose_multiple = True):
        """
        Constructor for File_browser objects.
        
        :param starting_dir: The pathlib Path of the directory to start from.
        :param show_hidden: Whether to show hidden files (those starting with a .)
        """
        starting_dir = starting_dir.resolve()
        # We store these options as a dict so we can give the same option to each node (and only change values once).
        self.options = {'show_hidden': show_hidden}
        super().__init__(Flaggable_tree_walker(Directory_node(starting_dir, starting_path = starting_dir, options = self.options)), can_choose_parents = can_choose_folders, can_choose_multiple = can_choose_multiple)


class File_selector(Enhanced_tree_selector):
    """
    A tree list box widget used to browse and select files.
    """
    
    show_hidden = Option(help = "Whether to show hidden files.", type = bool, default = False)

    def __init__(self, top, starting_dir = None, title = "File Browser", manual_widget_title = "Manual File Path", can_choose_folders = False, can_choose_multiple = True):
        """
        Constructor for File_selector objects.
        
        :param starting_dir: The starting directory that will be shown expanded.
        :param title: A title to display around this selector.
        """
        starting_dir = starting_dir if starting_dir is not None else pathlib.Path.cwd()
        browser = File_browser(starting_dir, show_hidden = self.show_hidden, can_choose_multiple = can_choose_multiple, can_choose_folders = can_choose_folders)
        browser.offset_rows = 1
        
        manual_widget = urwid.Edit(("body", "File: "))
        
        super().__init__(top, browser, manual_widget, title, manual_widget_title)
        
    def on_settings_change(self):
        """
        A method that will be called when settings have been changed.
        """
        self.browser.options['show_hidden'] = self.show_hidden
