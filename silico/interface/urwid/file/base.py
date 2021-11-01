# Base classes for file browser widgets.
# This code is adapted from the example file browser included with urwid, originally written by Rob Lanphier.

# General imports.
from pathlib import Path
import urwid

# Silico imports.
from silico.interface.urwid.tree.base import Flag_widget,\
    Flaggable_tree_list_box
from silico.interface.urwid.top import View
from silico.config.configurable.option import Option
import itertools


class File_tree_widget(Flag_widget):
    """
    Widget for individual files.
    """
    
    def get_display_text(self):
        """
        The text we'll use for display.
        
        Our key is our file path, so we can just use that.
        """
        return self.get_node().get_key()


class Empty_widget(urwid.TreeWidget):
    """
    A marker for expanded directories with no contents.
    """
    
    def get_display_text(self):
        return '(empty directory)'


class Error_widget(urwid.TreeWidget):
    """
    A marker for errors reading directories.
    """

    def get_display_text(self):
        return ('warningnode', "(error: {})".format(self.get_node().get_value()))


class Directory_widget(Flag_widget):
    """
    Widget for a directory.
    """
    
    def __init__(self, node):
        super().__init__(node)
        path = node.get_value()
        
        # Work out whether we should be expanded or not.
        self.expanded = path in itertools.chain(node.starting_path.parents, [node.starting_path])
        
        self.update_expanded_icon()        

    def get_display_text(self):
        """
        The text we'll use for display.
        
        In most cases we can use our key, which is our directory path, unless we're the root directory.
        """
        node = self.get_node()
        if node.get_depth() == 0:
            return "/"
        else:
            return node.get_key()


class File_node(urwid.TreeNode):
    """
    Metadata storage for individual files.
    """

    def __init__(self, path, starting_path, parent = None, show_hidden = False):
        """
        Constructor for File_node objects.
        
        :param path: Pathlib path of this node.
        :param path: The pathlib Path of the path that should be started at.
        :param parent: The parent of this node (next directory up).
        :param show_hidden: Whether to show hidden files.
        """
        self.show_hidden = show_hidden
        self.starting_path = starting_path
        depth = len(path.parts) -1
        key = path.name
        super().__init__(path, key = key, parent = parent, depth = depth)

    def load_parent(self):
        """
        Load the parent node of this node.
        """
        path = self.get_value()
        parent = Directory_node(path.parent, starting_path = self.starting_path, show_hidden = self.show_hidden)
        parent.set_child_node(self.get_key(), self)
        return parent

    def load_widget(self):
        """
        Load the widget we'll use for display.
        """
        return File_tree_widget(self)


class Empty_node(urwid.TreeNode):
    """
    A logical node to represent an empty directory.
    """
    
    def __init__(self, parent = None, key = None, depth = None):
        urwid.TreeNode.__init__(self, None, parent=parent, key=key, depth=depth)
    
    def load_widget(self):
        """
        Load the widget we'll use for display.
        """
        return Empty_widget(self)


class Error_node(urwid.TreeNode):
    """
    A logical node to represent a directory that failed to load.
    """
    
    def __init__(self, exception, parent=None, key=None, depth=None):
        """
        Constructor for Error_node objects.
        """
        urwid.TreeNode.__init__(self, exception, parent=parent, key=key, depth=depth)
    
    def load_widget(self):
        return Error_widget(self)


class Directory_node(urwid.ParentNode):
    """
    A logical node to represent a directory.
    """

    def __init__(self, path, starting_path, parent = None, show_hidden = False):
        """
        Constructor for Directory_node objects.
        
        :param path: The pathlib Path of this node.
        :param path: The pathlib Path of the path that should be started at.
        :param parent: The parent node, if any.
        :param show_hidden: Whether to show hidden files and directories.
        """
        self.show_hidden = show_hidden
        self.starting_path = starting_path
        
        if path == Path("/"):
            depth = 0
            key = None
            
        else:
            depth = len(path.parts) -1
            key = path.name
            
        # The number of child directories we have.
        self.num_directories = 0
            
        urwid.ParentNode.__init__(self, path, key = key, parent = parent, depth = depth)

    def load_parent(self):
        """
        Load the parent node of this node.
        """
        path = self.get_value()
        parent = Directory_node(path.parent, starting_path = self.starting_path, show_hidden = self.show_hidden)
        parent.set_child_node(self.get_key(), self)
        return parent

    def load_child_keys(self):
        """
        Load the child keys of this node.
        """
        # We want to keep our files and directories separate, because we will display directories first, then files afterwards.
        directories = []
        files = []
        
        try:
            path = self.get_value()
            # separate dirs and files
            for child in path.iterdir():
                # Ignore hidden files unless we've been told not to.
                if child.name[:1] == "." and not self.show_hidden:
                    continue
                
                if child.is_dir():
                    directories.append(child.name)
                
                elif child.is_file():
                    files.append(child.name)
                    
                else:
                    # A child can be neither a file nor a directory if it doesn't exist (or is a broken symlink).
                    pass
                
        except Exception as error:
            # Remove any directories or files we loaded before.
            directories = []
            files = []
            files.append(error)


        # Sort our files and directories.
        directories.sort()
        files.sort()
        
        # Keep track of how many directories we have.
        self.num_directories = len(directories)
        
        # Combine our keys together.
        keys = directories + files
        
        # If we've got no keys, add a single dummy key to trigger the empty node.
        if len(keys) == 0:
            keys = [None]
        
        return keys

    def load_child_node(self, key):
        """
        Load a Node for a given child key.
        
        :param key: The key of the child to load a node for.
        """
        index = self.get_child_index(key)
        if key is None:
            # An empty Node.
            return Empty_node(parent = self, key = key, depth = self.get_depth() +1)
        
        elif isinstance(key, Exception):
            # An error, use an error node.
            return Error_node(key, parent = self, key = key, depth = self.get_depth() +1)
            
        else:
            path = Path(self.get_value(), key)
            if index < self.num_directories:
                return Directory_node(path, self.starting_path, parent = self, show_hidden = self.show_hidden)
            
            else:
                return File_node(path, self.starting_path, parent = self, show_hidden = self.show_hidden)

    def load_widget(self):
        """
        Load the inner widget we'll use for display.
        """
        return Directory_widget(self)


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
        super().__init__(urwid.TreeWalker(Directory_node(starting_dir, starting_path = starting_dir, show_hidden = show_hidden)), can_choose_parents = can_choose_folders, can_choose_multiple = can_choose_multiple)


class File_selector(View):
    """
    A tree list box widget used to browse and select files.
    """

    def __init__(self, starting_dir, title = "File Browser"):
        """
        Constructor for File_selector objects.
        
        :param starting_dir: The starting directory that will be shown expanded.
        :param title: A title to display around this selector.
        """
        self.browser = File_browser(starting_dir)
        self.browser.offset_rows = 1
        
        super().__init__(self.browser, title)
        
class Coord_selector(File_selector):
    """
    A file selector for loading coordinate files.
    """

    charge = Option(help = "Forcibly set the molecular charge (an integer) of newly loaded files. If blank, the charge will be determined from each loaded file.", type = int, default = None)
    multiplicity = Option(help = "Forcibly set the molecular multiplicity (an integer) of newly loaded files. If blank, the multiplicity will be determined from each loaded file.", type = int, default = None)
    generate_3D = Option(help = "Whether to convert 2D coordinates to 3D by performing a rapid optimisation with molecular mechanics. Even if True, 3D coordinates will only be generated if it can be safely determined that the coordinates are not already in 3D.", type = bool, default = True)
