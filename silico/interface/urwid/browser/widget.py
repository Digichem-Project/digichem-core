# Widgets that control the appearance of the calculation browser.

import urwid
from urwid import TreeWidget

class Loader_widget(TreeWidget):
    """
    Mixin class for widgets that represent configurable loaders.
    """
    
    indent_cols = 4
    unexpanded_icon = urwid.Text('+')
    expanded_icon = urwid.Text('-')
    
    def load_inner_widget(self):
        """
        Load the inner widget we use for displaying text (not including the expanded icon etc).
        """
        node = self.get_node()
        loader = node.get_value()[-1]
        
        return Node_widget(self.get_display_text(), bold = not loader.partial)
    
    def get_display_text(self):
        """
        Get the text to display at this node.
        
        For leaf/half-leaf nodes, we show both the index of the leaf and the TAG/ALIAS of this loader.
        For parent nodes, we just show the alias.
        """
        node = self.get_node()
        loader = node.get_value()[-1]
        
        if not loader.partial:
            # We have an index, get it.
            
            # We need a loader path from our bottom node to the top node.
            path = node.build_loader_path()[-1]
            
            # Now get the index of that path.
            index = path[0].index_of_path(path)
            return str("[{}] {}".format(index, loader.ALIAS))
        
        else:
            return str(loader.ALIAS)


class Node_widget(urwid.WidgetWrap):
    """
    The inner widget which displays information about each node.
    """
    
    # The attribute we normally use for display when not in focus.
    normal_attr = 'body'
    # The attribute we use for display if we are a bold node.
    bold_attr = 'boldnode'
    # The attribute we use for display when we are in focus.
    focus_attr = 'focus'
    
    def __init__(self, display_text, bold = True):
        """
        Constructor for Node_widget objects.
        
        :param display_text: The text to display.
        """
        super().__init__(urwid.AttrMap(urwid.Text(display_text), self.normal_attr if not bold else self.bold_attr, self.focus_attr))
    
    def selectable(self):
        """
        Node widgets are always selectable.
        """
        return True
    
    def keypress(self, size, key):
        """
        Dummy keypress function; Node_widgets do not handle keypress themselves.
        """
        return key


class Loader_leaf_widget(Loader_widget):
    """
    Display widget for leaf nodes.
    """
   
    def selectable(self):
        """
        We are always selectable.
        """
        return True


class Loader_parent_widget(Loader_widget):
    """
    Display widget for nodes.
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # We start unexpanded so we don't have to generate the entire tree at the start.
        self.expanded = False
        self.update_expanded_icon()


class Loader_top_widget(Loader_widget):
    """
    Widget for displaying the top-most node.
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # We start unexpanded so we don't have to generate the entire tree at the start.
        self.expanded = False
        self.update_expanded_icon()
    
    def get_display_text(self):
        return "Methods"