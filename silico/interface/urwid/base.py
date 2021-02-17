import urwid

class Top(urwid.WidgetPlaceholder):
    """
    A placeholder widget used for switching the top-most widget being rendered by urwid.
    """
    
    def __init__(self, original_widget):
        """
        Constructor for Top widgets.
        
        :param original_widget: The main widget to wrap around. Must be given (use a placeholder if you require).
        """
        super().__init__(original_widget)
        self.stack = []
        
    def swap(self, original_widget):
        """
        Set a new widget as the top-most.
        """
        self.stack.append(self.original_widget)
        self.original_widget = original_widget
        
    def back(self):
        """
        Remove the current top-most widget and set the last in its place.
        
        If there are no more widgets to go back to, urwid.ExitMainLoop will be raised.
        
        :return: The widget just removed, for convenience.
        """
        try:
            self.original_widget = self.stack.pop()
        except IndexError:
            # The stack is empty.
            raise urwid.ExitMainLoop()
        
    def keypress(self, size, key):
        """
        Handler for keypress events.
        """
        if key == "esc":
            self.back()
        else:
            return super().keypress(size, key)
        
class Swappable():
    """
    Mix-in class for widgets that can be swapped with Top.
    """
    
    def __init__(self):
        self.top = None
        
class Window(urwid.Frame):
    """
    Container widget for high-level widgets that emulates a window.
    """
    
    def __init__(self, body, title = None, header_attr = "header", help = None, footer_attr = "footer"):
        """
        Constructor for Swappable widgets.
        
        :param title: Optional title to display when swapped to this 
        """
        if title is not None:
            header = urwid.AttrMap(urwid.Text(title, align = "center"), header_attr)
        else:
            header = None
        
        if help is not None:
            footer = urwid.AttrMap(urwid.Text(help, align = "center"), footer_attr)
        else:
            footer = None
            
        super().__init__(body, header = header, footer = footer)
        
        
        
    