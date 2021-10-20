# General imports.
import urwid

# Silico imports.


class Section(urwid.AttrMap):
    """
    A possibly selectable sub-section.
    """
    
    def __init__(self, body, title, focusable = True):
        """
        Constructor for urwid sections.
        
        :param body: The main widget we will display.
        :param title: The title of the section.
        :param focusable: Whether this section can take focus.
        """
        attrs = ["section"]
        if focusable:
            attrs.append("focus--section")
            
        linebox = urwid.LineBox(urwid.AttrMap(body, "body"), title, title_align = "left")
            
        super().__init__(linebox, *attrs)

        