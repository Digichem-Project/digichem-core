# General imports.
import urwid

# Silico imports.


class Section(urwid.AttrMap):
    """
    A possibly selectable sub-section.
    """
    
    def __init__(self, body, title = "", focusable = True):
        """
        Constructor for urwid sections.
        
        :param body: The main widget we will display.
        :param title: The title of the section.
        :param focusable: Whether this section can take focus.
        """
        attrs = ["section"]
        if focusable:
            attrs.append("focus--section")
        
        
        self._inner = urwid.AttrMap(body, "body")
        linebox = urwid.LineBox(self._inner, title, title_align = "left")
            
        super().__init__(linebox, *attrs)
        
    @property
    def inner_body(self):
        """
        Get the inner widget that we are wrapping.
        """
        return self.base_widget
    
    @inner_body.setter
    def inner_body(self, value):
        """
        Change the inner widget that we are wrapping.
        """
        self._inner.original_widget = value