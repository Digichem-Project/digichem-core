# General imports.
import urwid
import functools

# Silico imports.
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.layout import Pane


class Pages(Tab_pile):
    """
    A widget that can switch between a number of different pages.
    """
    
    def __init__(self, pages, title = "Page Selector"):
        """
        
        :param title: The title to display.
        :param pages: A list of pages to switch between, where each 'page' is a tuple of (title, widget).
        """
        self.title = title
        self.pages = pages
        self.top = urwid.WidgetPlaceholder(pages[0][1])
        
        controls = urwid.Columns(self.get_controls(), dividechars = 1)
        super().__init__([
            ("pack", Pane(controls, title = self.title)),
            self.top
        ])
        
    
    def get_controls(self):
        """
        Get a list of button controls that we'll use to switch between different pages.
        
        :returns: A list of buttons.
        """
        controls = []
        
        for page_index, (page_title, widget) in enumerate(self.pages):
            switch_button = urwid.Button(page_title, functools.partial(self.switch_page, index = page_index))
            controls.append(urwid.AttrMap(switch_button, "button--small", "button--small--focus"))
            
        return controls
        
        
    def switch_page(self, button, index):
        """
        Switch the currently visibly page.
        
        :param index: The index of the page to switch to.
        """
        self.top.original_widget = self.pages[index][1]
        
    
        