# Logical storage for file browser widgets.

# General imports.
import urwid
from pathlib import Path

# Silico imports.
from silico.interface.urwid.file.widget import File_tree_widget, Empty_widget,\
    Error_widget, Directory_widget


class File_node(urwid.TreeNode):
    """
    Metadata storage for individual files.
    """

    def __init__(self, path, starting_path, parent = None, options = None):
        """
        Constructor for File_node objects.
        
        :param path: Pathlib path of this node.
        :param path: The pathlib Path of the path that should be started at.
        :param parent: The parent of this node (next directory up).
        :param show_hidden: Whether to show hidden files.
        """
        self.options = options if options is not None else {}
        self.starting_path = starting_path
        depth = len(path.parts) -1
        key = path.name
        super().__init__(path, key = key, parent = parent, depth = depth)

    def load_parent(self):
        """
        Load the parent node of this node.
        """
        path = self.get_value()
        parent = Directory_node(path.parent, starting_path = self.starting_path, options = self.options)
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

    def __init__(self, path, starting_path, parent = None, options = None):
        """
        Constructor for Directory_node objects.
        
        :param path: The pathlib Path of this node.
        :param starting_path: The pathlib Path of the path that should be started at.
        :param parent: The parent node, if any.
        :param options: A dictionary of options.
        """
        self.options = options if options is not None else {}
        # Whether we are currently showing hidden children.
        self._showing_hidden = None
        
        # Path to the file we are showing at startup.
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
        
    @property
    def show_hidden(self):
        """
        Whether to show hidden files and directories.
        """
        return self.options.get('show_hidden')

    def load_parent(self):
        """
        Load the parent node of this node.
        """
        path = self.get_value()
        parent = Directory_node(path.parent, starting_path = self.starting_path, options = self.options)
        parent.set_child_node(self.get_key(), self)
        return parent
    
    def get_child_keys(self, reload=False):
        """
        Return a possibly ordered list of child keys
        """
        if self._child_keys is None or self.show_hidden != self._showing_hidden or reload:
            self._child_keys = self.load_child_keys()
        return self._child_keys

    def load_child_keys(self):
        """
        Load the child keys of this node.
        """
        # We want to keep our files and directories separate, because we will display directories first, then files afterwards.
        directories = []
        files = []
        
        self._showing_hidden = self.show_hidden
        
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
                return Directory_node(path, self.starting_path, parent = self, options = self.options)
            
            else:
                return File_node(path, self.starting_path, parent = self, options = self.options)

    def load_widget(self):
        """
        Load the inner widget we'll use for display.
        """
        return Directory_widget(self)