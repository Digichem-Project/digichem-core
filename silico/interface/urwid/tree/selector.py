from silico.interface.urwid.top import View
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.section import Section
import urwid


class Enhanced_tree_selector(View):
    """
    ABC for widgets that contain both a tree widget for selecting something and also a manual text-like widget
    """
    
    def __init__(self, top, browser, manual_widget, title, manual_widget_title):
        """
        Constructor for Enhanced_tree_selector objects.
        """
        self.browser = browser
        self.manual_widget = manual_widget
        self.manual_widget_title = manual_widget_title
        
        super().__init__(top, self.browser, title)
        
    def _get_body(self):
        """
        Get the widget we'll use for display.
        """
        main_body = super()._get_body()
        
        # We will wrap our body in a Pile and add our manual widget underneath in.
        return Tab_pile([
            main_body,
            (3, Section(urwid.Filler(urwid.AttrMap(self.manual_widget, "editable")), self.manual_widget_title))
        ])
        
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