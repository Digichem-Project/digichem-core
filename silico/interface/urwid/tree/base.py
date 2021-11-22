import urwid


class Flag_widget(urwid.TreeWidget):
    """
    A tree widget which can be 'flagged' to show it has been selected.
    """
    
    unexpanded_icon = urwid.Text('+')
    expanded_icon = urwid.Text('-')
    # Attributes to use for display.
    flagged_attr = "node--flagged"
    flagged_focus_attr = "node--flagged--focus"
    normal_attr = "body"
    focus_attr = "node--focus"

    def __init__(self, node):
        """
        Constructor for Flag_widget objects.
        
        :param node: The logical node we are rendering.
        """
        super().__init__(node)
        
        # We'll use an attr wrap to change formatting.
        self.wrapper = urwid.AttrWrap(self._w, None)
        self._w = self.wrapper
        self.flagged = False
        self.update_attr()

    def selectable(self):
        return True

    def update_attr(self):
        """
        Update the attributes of self.wrapper based on self.flagged.
        """
        if self.flagged:
            self.wrapper.attr = self.flagged_attr
            self.wrapper.focus_attr = self.flagged_focus_attr
        else:
            self.wrapper.attr = self.normal_attr
            self.wrapper.focus_attr = self.focus_attr
            
            
    def keypress(self, size, key):
        """
        Handle expand & collapse requests (non-leaf nodes).
        """
        if self.is_leaf:
            return key
 
        if key == 'right':
            # Expand if we're not already expanded.
            if not self.expanded:
                self.expanded = True
                self.update_expanded_icon()
                 
            else:
                # Already expanded, propagate up.
                return key
             
        elif key == "left":
            # Collapse if we're not already collapsed.
            if self.expanded:
                self.expanded = False
                self.update_expanded_icon()
                 
            else:
                # Already collapsed, propagate up.
                return key
             
        else:
            return super().keypress(size, key)

class Flaggable_tree_walker(urwid.TreeWalker):
    """
    Tree walker object for Flaggable tree list boxes.
    """
    
    def find_real_parent(self, node):
        """
        Find the first parent of a node that has not been removed.
        
        If the given node has not been removed, then it is returned.
        :param node: The node to start from.
        :returns: A node, possibly the same.
        """        
        while True:
            parent = node.get_parent()
            
            if parent is None or node.get_key() in parent.get_child_keys():
                # No more parents, or our focus is a real child of its parent.
                return node
            
            # Go again.
            node = parent
    
    def get_focus(self):
        """
        """
        # Get the current focus.
        focus = self.focus
        
        # Check if it's real.
        real = self.find_real_parent(focus)
            
        # If the focus has changed, update.
        if real != focus:
            self.focus = real
                
        return urwid.TreeWalker.get_focus(self)

class Flaggable_tree_list_box(urwid.TreeListBox):
    """
    Widget for displaying a tree that can be selected.
    """
    
    def __init__(self, walker, can_choose_parents = False, can_choose_multiple = True):
        """
        Constructor for File_browser objects.
        
        :param walker: A tree walker class to use as our body.
        :param can_choose_parents: Whether the user can select parent nodes or only children.
        :param can_choose_multiple: Whether the user can select only one item or multiple.
        """
        super().__init__(walker)
        
        self.can_choose_parents = can_choose_parents
        self.can_choose_multiple = can_choose_multiple
        
        # The nodes that have been selected.
        self.selected_nodes = []
        
    @property
    def selected(self):
        """
        A property resolving to the values that are currently selected.
        """
        return [node.get_value() for node in self.selected_nodes]
        
    def is_selectable(self, node):
        """
        Determine whether a given node is selectable.
        """
        return not hasattr(node, 'has_children') or self.can_choose_parents
        
    def keypress(self, size, key):
        """
        Handle keypress events.
        """
        if key in [' ', 'enter']:
            # The user wants to select a node, see if we can.
            focus_node = self.body.focus
            focus_widget = focus_node.get_widget()
            
            if self.is_selectable(focus_node):
                self.select(focus_node, focus_widget)
                
        elif key == "right":
            # Check to see if the widget is expanded.
            if self.body.focus.get_widget().expanded:
                # Move to child.
                self.move_focus_to_child()
            
            else:
                return super().keypress(size, key)
            
        elif key == "left":
            # Check to see if the widget is collapsed.
            if not self.body.focus.get_widget().expanded:
                # Move to parent.
                self.move_focus_to_parent(size)
                
            else:
                return super().keypress(size, key)

        else:
            return super().keypress(size, key)

    def move_focus_to_child(self):
        """
        Move focus to the first child of the widget in focus.
        """

        focus_node = self.body.focus

        # Try and get the first child of the node in focus.
        # This could fail because: 1) The node is a leaf, 2) the node is a parent with no children.
        try:
            child_node = focus_node.get_first_child()
            
        except Exception:
            return
        
        if child_node.get_widget().selectable():
            self.set_focus(child_node, "above")

    def select(self, focus_node, focus_widget):
        """
        Select (or deselect) the node currently in focus.
        """
        # Check if the node is in our list.        
        if focus_node not in self.selected_nodes:
            # Append to our list.
            self.selected_nodes.append(focus_node)
            focus_widget.flagged = True
            
        else:
            # Already in our list, remove it.
            self.selected_nodes.pop(self.selected_nodes.index(focus_node))
            focus_widget.flagged = False
            
        # Toggle the attribute for the child.
        focus_widget.update_attr()
        
    def reset(self):
        """
        Unselect any nodes that have already been selected.
        """        
        for node in self.selected_nodes:
            # Remove the flagged status.
            node.get_widget().flagged = False
            node.get_widget().update_attr()
            
        self.selected_nodes = []

