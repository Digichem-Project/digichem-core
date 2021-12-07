# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.swap.top import Top


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
        self.loop = None
            
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