# General imports.
import urwid

# Silico imports.


class Section(urwid.AttrMap):
    """
    A possibly selectable sub-section.
    """
    
    tlcorner = '┌'
    tline = '─'
    lline = '│'
    trcorner = '┐'
    blcorner = '└'
    rline = '│'
    bline = '─'
    brcorner = '┘'
    title_align = "left"
    
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
        linebox = urwid.LineBox(
            self._inner, title,
            title_align = self.title_align,
            tlcorner = self.tlcorner,
            tline = self.tline,
            lline = self.lline,
            trcorner = self.trcorner,
            blcorner = self.blcorner,
            rline = self.rline,
            bline = self.bline,
            brcorner = self.brcorner
        )
            
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
        
class Sub_section(Section):
    """
    A section with a less intrusive border.
    """
    
    tlcorner = '─'
    tline = '─'
    lline = ''
    trcorner = '─'
    blcorner = ''
    rline = ''
    bline = ''
    brcorner = ''
    title_align = "center"
    