# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.top import Top

class Window(urwid.Frame):
    """
    Container widget for high-level widgets that emulates a window.
    """
    
    def __init__(self, body = None, title = "", help = ""):
        """
        Constructor for Swappable widgets.
        
        :param body: The body of the window, must be a Swappable class. The body can be accessed (and changed) at self.top.original_widget. If None is give, the body must be set before rendering.
        :param title: Title of this window.
        :param help: Help text to display this window.
        """
        self.header_text = urwid.Text(title, align = "center")
        self.footer_text = urwid.Text(help, align = "center")
        
        self.top = Top(body)
            
        super().__init__(urwid.AttrMap(self.top, "body"), header = urwid.AttrMap(self.header_text, "header"), footer = urwid.AttrMap(self.footer_text, "footer"))
    
    def unhandled_input(self, key):
        if key in ('q','Q'):
            raise urwid.ExitMainLoop()
    
    def run_loop(self, palette):
        """
        Run an urwid loop using this Window as the top-most frame.
        
        :param palette: The palette to display with.
        """
        self.loop = urwid.MainLoop(self, palette, unhandled_input=self.unhandled_input)
        self.loop.run()
        
        
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
        