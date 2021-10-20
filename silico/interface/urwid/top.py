# Silico imports.
from silico.interface.urwid.wrapper import Confirm_or_cancel, Confirm_settings_cancel, Confirm

# General imports.
import urwid
from silico.interface.urwid.base import Section
from silico.config.configurable.base import Configurable
from silico.interface.urwid.configurable import Configurable_editor


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
        Wrap a widget buttons and then set it as the top-most widget.
        
        :param cancel_callback: A function to call when the cancel button is pressed.
        :param submit_callback: A function to call when the submit button is pressed.
        """
        # First decide which kind of wrapper to use.
        # If we have options, use a Confirm_settings_cancel.
        if isinstance(original_widget, View) and original_widget.has_settings:
            window = Confirm_settings_cancel(original_widget, top = self, settings_editor = original_widget.get_settings_editor(), cancel_callback = cancel_callback, submit_callback = submit_callback)
            
        elif submit_callback is not None:
            window = Confirm_or_cancel(original_widget, top = self, cancel_callback = cancel_callback, submit_callback = submit_callback)
            
        else:
            window = Confirm(original_widget, top = self, submit_callback = cancel_callback)
        
        self.swap(window)

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


class View(Section, Configurable):
    """
    A widget designed to be shown inside a swapping window.
    """
    
    def __init__(self, body, title, focusable = True):
        self._settings_editor = None
        Section.__init__(self, body, title, focusable = focusable)
        Configurable.__init__(self, True)
    
    @property
    def has_settings(self):
        """
        Does this View have any configurable options set on it?
        """
        return len(self.OPTIONS) != 0
    
    def get_settings_editor(self, reload = False):
        """
        Return the widget we use to change our settings.
        """
        if self._settings_editor is None or reload:
            self._settings_editor = self.load_settings_editor()
            
        return self._settings_editor
    
    def load_settings_editor(self):
        """
        Load the browser widget we'll use to change our settings.
        """
        return Configurable_editor(self, "Settings")
    