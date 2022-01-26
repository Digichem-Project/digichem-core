# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.swap.swappable import Swappable
from silico.interface.urwid.layout import Pane


class Two_part_tree_selector(Swappable):
    """
    ABC for tree widgets that show an additional widget.
    """
    
    def __init__(self, top, browser, secondary_widget, browser_title, secondary_title):
        """
        """
        self.browser = browser
        self.secondary_widget = secondary_widget
        self.secondary_title = secondary_title
        
        super().__init__(top, self.browser, title = browser_title)
        
    def _get_body(self):
        """
        Get the widget we'll use for display.
        """
        main_body = super()._get_body()
        
        # We will wrap our body in a Pile and add our manual widget underneath in.
        return Tab_pile([
            main_body,
            (3, Pane(urwid.Filler(urwid.AttrMap(self.secondary_widget, "editable")), self.secondary_title))
        ])
        

class Enhanced_tree_selector(Two_part_tree_selector):
    """
    ABC for widgets that contain both a tree widget for selecting something and also a manual text-like widget
    """
        
    @property
    def manual_widget(self):
        """
        Map for manual_widget -> secondary_widget
        """
        return self.secondary_widget
    
    @manual_widget.setter
    def manual_widget(self, value):
        """
        Map for manual_widget -> secondary_widget
        """
        self.secondary_widget = value
        
    @property
    def selected(self):
        """
        A property resolving to the values that are currently selected.
        """
        selected = self.browser.selected
        
        # Also add the value from our manual edit widget (if not empty).
        if self.manual_widget.get_edit_text() != "":
            selected.append(self.manual_widget.get_edit_text())
            
        return selected
    
    def reset(self):
        """
        Reset this selector so nothing is currently selected.
        """
        self.browser.reset()
        
        # Also reset our manual widget.
        self.manual_widget.set_edit_text("")