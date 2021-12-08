# Display widgets for the file browser.

# General imports.
import urwid
import itertools

# Silico imports.
from silico.interface.urwid.tree.base import Flag_widget


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
        
    @property
    def expanded(self):
        """
        Whether this directory is currently expanded.
        """
        return self._expanded
    
    @expanded.setter
    def expanded(self, expand):
        """
        Whether this directory is currently expanded.
        """
        if not expand:
            # We are collapsing.
            # If we've been ask to refresh when toggled, delete our child keys so we'll reload them.
            node = self.get_node()
            if node.options.get('refresh'):
                node._child_keys = None
        self._expanded = expand

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