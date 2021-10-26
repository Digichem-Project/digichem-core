# Widgets that control the appearance of the calculation browser.

# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.tree.base import Flag_widget


class Loader_widget(Flag_widget):
    """
    Mixin class for widgets that represent configurable loaders.
    """
    
    # Attributes to use for display.
    warning_normal_attr = "node--warning"
    warning_focus_attr = "node--warning--focus"
    
    
    def __init__(self, node):
        self._text = None
        super().__init__(node)
        
    @property
    def warnings(self):
        """
        A list of warnings specified for the loaders that we represent.
        """
        warnings = []
        for loader in self.get_node().get_value():
            warning = loader.config.get("warning", None) 
            if warning is not None:
                warnings.append(warning) 
        return warnings
        
    def update_attr(self):
        """
        Update the attributes of self.wrapper based on self.flagged.
        """
        if self.flagged:
            self.wrapper.attr = self.flagged_attr
            self.wrapper.focus_attr = self.flagged_focus_attr
        else:
            if len(self.warnings) > 0:
                self.wrapper.attr = self.warning_normal_attr
                self.wrapper.focus_attr = self.warning_focus_attr
            else:
                self.wrapper.attr = self.normal_attr
                self.wrapper.focus_attr = self.focus_attr
    
    def load_inner_widget(self):
        """
        Load the inner widget we use for displaying text (not including the expanded icon etc).
        """        
        return Node_widget(self.get_display_text())
    
    def get_display_text(self, reload = False):
        """
        Get the text to display at this node.
        
        For leaf/half-leaf nodes, we show both the index of the leaf and the TAG/ALIAS of this loader.
        For parent nodes, we just show the alias.
        """
        if self._text is None or reload:
            self._text = self.load_display_text()
            
        return self._text
        
    def load_display_text(self):
        """
        Load the text to display at this node.
        """
        node = self.get_node()
        loader = node.get_value()[-1]
        
        if not loader.partial:
            # We have an index, get it.
            
            # We need a loader path from our bottom node to the top node.
            path = node.build_loader_path()[-1]
            
            # Now get the index of that path.
            index = path[0].index_of_path(path)
            
            # If we are also a destination, we will resolve ourself so we can get our status.
            if loader.TYPE == "destination":
                destination = node.resolve()
                
                try:
                    status = destination.status()
                    
                except Exception:
                    # Could be because we couldn't retrieve status, or because this destination doesn't support status.
                    status = None
            
            else:
                status = None
                
            text = "[{}] {}".format(index, loader.ALIAS)
            if status is not None:
                text += " ({})".format(status)
            
            return text
        
        else:
            return str(loader.ALIAS)


class Node_widget(urwid.WidgetWrap):
    """
    The inner widget which displays information about each node.
    """
            
    def __init__(self, display_text):
        """
        Constructor for Node_widget objects.
        
        :param display_text: The text to display.
        """
        super().__init__(urwid.Text(display_text))
    
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
        self.expanded = False
        self.update_expanded_icon()
    
    def get_display_text(self):
        return "Methods"