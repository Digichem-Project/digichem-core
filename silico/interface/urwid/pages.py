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
    
    def __init__(self, pages, title = "Page Selector ", max_width = 50):
        """
        
        :param title: The title to display.
        :param pages: An (ordered) dict of pages to switch between, where each item is a page to switch between and each key is the title of the page.
        """
        self.title = title
        self.pages = pages
        # Set our first page as the one shown by default.
        #self.top = urwid.WidgetPlaceholder(self.pages[list(self.pages.keys())[0]])
        self.page = urwid.WidgetPlaceholder(None)
        
        # Buttons to change page.
        self.buttons = self.get_controls()
        cols = self.calculate_min_cell_width(self.buttons, max_width)
        
        controls = urwid.GridFlow(self.buttons, cell_width = cols, h_sep = 1, v_sep = 0, align = "center")
        super().__init__([
            ("pack", Pane(controls, title = self.title)),
            self.page
        ])
        
        self.switch_page(list(self.pages.keys())[0])
    
    def calculate_min_cell_width(self, items, max_width):
        """
        """
        # We need to work out how many columns to give each of our page buttons.
        # We can do this by repeatedly calling rows() on each of the buttons until they all report that they only need a single row.
        # This is probably v wasteful tho.
        cols = 0
        while cols < max_width:
            cols += 1
            too_small = False
            for item in items:
                if item.rows((cols,)) > 1:
                    too_small = True
                    break
                
            if not too_small: 
                # If we get this far we can stop.
                return cols
    
    def get_controls(self):
        """
        Get a list of button controls that we'll use to switch between different pages.
        
        :returns: A list of buttons.
        """
        controls = []
        
        def button_callback(button, title):
            self.switch_page(title)
        
        for page_title, widget in self.pages.items():
            placeholder = urwid.WidgetPlaceholder(None)
            attrmap = urwid.AttrMap(placeholder, "button--small", "button--small--focus")
            switch_button = urwid.Button(page_title, functools.partial(button_callback, title = page_title))
            placeholder.original_widget = switch_button
            controls.append(attrmap)
            
        return controls
    
    def set_button_selected(self, button_attrmap):
        """
        Set one of the page selector buttons as selected.
        
        :parma button_attrmap: The selected 'button', wrapped by an attrmap.
        """
        # First, reset all buttons.
        for button in self.buttons:
            button.set_attr_map({None: 'button--small'})
            button.set_focus_map({None: 'button--small--focus'})
            
        # Now set our button.
        button_attrmap.set_attr_map({None: 'button--small--selected'})
        button_attrmap.set_focus_map({None: 'button--small--selected--focus'})
    
    @property
    def current_page(self):
        """
        Get the title of the page that is currently visible.
        """
        page_index = list(self.pages.values()).index(self.page.original_widget)
        return list(self.pages.keys())[page_index]
        
    def switch_page(self, title):
        """
        Switch the currently visibly page.
        
        :param title: The title of the page to switch to.
        """
        self.page.original_widget = self.pages[title]
        # Change the attr
        self.set_button_selected(self.buttons[list(self.pages).index(title)])
        
    
        