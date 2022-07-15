"""General widget classes for laying out content."""

import urwid


class Window(urwid.Frame):
    """
    A container widget that contains all on-screen content.
    
    A window contains a header and footer, which are generally static, and an inner body which can swap between various other widgets.
    """
    
    def __init__(self, top = None, title = "", help = ""):
        """
        Constructor for Window widgets.
        
        :param top: A Top widget to use as the main display widget.
        :param title: Title of this window.
        :param help: Help text to display this window.
        """
        self.header_text = urwid.Text(title, align = "center")
        self.footer_text = urwid.Text(help, align = "center")
        
        self.top = top
        self.loop = None
            
        super().__init__(urwid.AttrMap(self.top, "body"), header = urwid.AttrMap(self.header_text, "header"), footer = urwid.AttrMap(self.footer_text, "footer"))
    
#     def unhandled_input(self, key):
#         if key in ('q','Q'):
#             raise urwid.ExitMainLoop()
    
    def run_loop(self, palette):
        """
        Run an urwid loop using this Window as the top-most frame.
        
        :param palette: The palette to display with.
        """
        #self.loop = urwid.MainLoop(self, palette, unhandled_input=self.unhandled_input)
        self.loop = urwid.MainLoop(self, palette)
        self.loop.run()
        

class Pane(urwid.AttrMap):
    """
    A container widget that contains one 'section' of content within a wider window.
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
        :param title: The title of the pane.
        :param focusable: Whether this pane can take focus.
        """
        attrs = ["pane"]
        if focusable:
            attrs.append("pane--focus")
        
        
        self._pane_inner_body = urwid.AttrMap(body, "body")
        linebox = urwid.LineBox(
            urwid.AttrMap(body, "body"), title,
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
        return self._pane_inner_body.original_widget
    
    @inner_body.setter
    def inner_body(self, value):
        """
        Change the inner widget that we are wrapping.
        """
        self.original_widget.original_widget.original_widget = value


class Sub_pane(Pane):
    """
    A container widget similar to a Pane but with a more subtle border.
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
    