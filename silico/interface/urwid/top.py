# Silico imports.
from silico.interface.urwid.wrapper import Confirm_or_cancel

# General imports.
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

    def swap_into_window(self, original_widget, cancel_callback = None, submit_callback = None):
        """
        Wrap a widget with cancel and submit buttons and then set it as the top-most widget.
        
        :param cancel_callback: A function to call when the cancel button is pressed.
        :param submit_callback: A function to call when the submit button is pressed.
        """
        self.swap(Swapping_window(original_widget, top = self, cancel_callback = cancel_callback, submit_callback = submit_callback))

    def back(self, number = 1):
        """
        Remove the current top-most widget and set the last in its place.
        
        If there are no more widgets to go back to, urwid.ExitMainLoop will be raised.
        
        :param number: The number of times to go back.
        :return: The widget just removed, for convenience.
        """
        for iter in range(0, number):
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
        
class Swapping_window(Confirm_or_cancel):
    """
    A window used when swapping widgets, provides controls to return to the previous window.
    """