# General imports.
import urwid


class Setedit_widget_parent_mixin():
    """
    Mixin for setting editor-like widgets that contain other widgets.
    """
    
    @property
    def child_widgets(self):
        """
        A shortcut to the property where the child widgets are stored.
        
        This must be defined by the inheriting class because it depends on the type of the inheriting class, for example self.body for listboxes.
        """
        raise NotImplementedError("Implement in subclass")
    
    @property
    def child_setedit_widgets(self):
        """
        A short-cut property to provide access to the data containing edit widgets that are children of this widget.
        """
        return {child_widget.setedit.title: child_widget for child_widget in self.child_widgets if hasattr(child_widget, "value")}
    
    
    def load_child_widgets(self, child_setedits):
        """
        Create a list of child widgets from a list of child setedits.
        """
        child_widgets = []
        
        for child_setedit in child_setedits:
            # TODO: This is weird, surely the setedit should give us its child (via .get_widget() or similar?)
            #fields.append(self.class_from_type(child_setedit.vtype)(child_setedit, has_divider = len(fields) != len(children)-1))
            child_widgets.append(child_setedit.get_widget())
            # Also add a divider if we are not at the end
            if len(child_widgets) < (len(child_setedits)*2)-1:
                child_widgets.append(urwid.Divider("-"))
            
        return child_widgets
    
    def discard(self):
        """
        Discard any changes made to the child widgets.
        """
        for child_widget in self.child_setedit_widgets.values():
            child_widget.discard()
            