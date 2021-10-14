import urwid
from urwid.graphics import LineBox
from silico.interface.urwid.misc import Tab_pile

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

    def swap_into_window(self, original_widget, cancel_callback = None, submit_callback = None):
        """
        Wrap a widget with cancel and submit buttons and then set it as the top-most widget.
        
        :param cancel_callback: A function to call when the cancel button is pressed.
        :param submit_callback: A function to call when the submit button is pressed.
        """
        self.swap(Swapping_window(original_widget, top = self, cancel_callback = cancel_callback, submit_callback = submit_callback))

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
        # Allow children to intercept first.
        key = super().keypress(size, key)
        
        if key == "esc":
            self.back()
        
        return key
        
class Swappable():
    """
    Mix-in class for widgets that can be swapped with Top.
    """
    
    def __init__(self):
        self.top = None
        
class Swapping_window(Tab_pile):
    """
    A window used when swapping widgets, provides controls to return to the previous window.
    """
    
    def __init__(self, body, top, cancel_callback = None, submit_callback = None):
        self.body = body
        self.top = top
        self.submit_callback = submit_callback
        self.cancel_callback = cancel_callback
        
        # A value indicating which of our two buttons was pressed.
        self.retval = None
        
        controls = urwid.Columns([
            urwid.AttrMap(urwid.Button("Back", lambda button: self.cancel()), "button--normal", "button--focus"),
            urwid.AttrMap(urwid.Button("Submit", lambda button: self.submit()), "button--good", "button--focus")
        ])
        
    def keypress(self, size, key):
        """
        Handler for keypress events.
        """
        # Allow children to intercept first.
        key = super().keypress(size, key)
        
        if key == "esc":
            self.cancel()
        
        return key
        
    def cancel(self):
        """
        Function called when the back button is pressed.
        """
        self.retval = False
        if self.cancel_callback is not None:
            self.cancel_callback()
        self.top.back()
        
    def submit(self):
        """
        Function called when the submit button is pressed.
        """
        self.retval = True
        if self.submit_callback is not None:
            self.submit_callback()
        self.top.back()


class Window(urwid.Frame):
    """
    Container widget for high-level widgets that emulates a window.
    """
    
    def __init__(self, body, title = "", help = ""):
        """
        Constructor for Swappable widgets.
        
        :param body: The body of the window, must be a Swappable class.
        :param title: Title of this window.
        :param help: Help text to display this window.
        """
        self.header_text = urwid.Text(title, align = "center")
        self.footer_text = urwid.Text(help, align = "center")
        
        self.swappable_top = Top(body)
        body.top = self.swappable_top
            
        super().__init__(urwid.AttrMap(self.swappable_top, "body"), header = urwid.AttrMap(self.header_text, "header"), footer = urwid.AttrMap(self.footer_text, "footer"))
    
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
    