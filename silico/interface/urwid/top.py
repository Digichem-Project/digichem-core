# Silico imports.
from silico.interface.urwid.wrapper import Confirm_or_cancel, Confirm_settings_cancel, Confirm
from silico.interface.urwid.section import Section, Sub_section
from silico.config.configurable.base import Configurable
from silico.interface.urwid.setedit.base import Settings_editor
from silico.interface.urwid.setedit.view import View_browser
from silico.interface.urwid.dialogue import Output_dialogue

# General imports.
import urwid
from urwid.widget import WidgetWrap
from urwid.decoration import WidgetPlaceholder


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
        
        self.current_top = WidgetPlaceholder(original_widget)
        
        # We keep track of two separate widget stacks, one for normal widgets and one for dialogues (which appear on top of normal widgets).
        self.stack = []
        self.dialogue_stack = []
        
        # A widget we'll use to popup and show program output.
        self.output_widget = Output_dialogue(self)
        
    def set_top(self, original_widget):
        """
        Set the currently visible widget.
        """
        self.original_widget = original_widget
        self.current_top.original_widget = original_widget
        
    def output(self, message, error = False):
        """
        Add text to the output widget and show it on top.
        
        :param message: The text to add.
        :param error: Whether to show alternative formatting indicating an error.
        """
        # Urwid doesn't render tab characters, so we'll replace those.
        message = message.replace('\t', '  ')
        
        # First, add our message to our widget.
        self.output_widget.output(message, error)
        
        # Then swap to it if it's not already the top.
        if len(self.dialogue_stack) == 0 or self.dialogue_stack[-1].top_w != self.output_widget:
            self.popup(self.output_widget, width = ('relative', 100), height = ('relative', 100), left = 1, top = 1, right = 1, bottom = 2)
            
    def update_view(self):
        """
        Update the currently visible widget.
        """
        if len(self.dialogue_stack) > 0:
            self.original_widget = self.dialogue_stack[-1]
        
        else:
            self.original_widget = self.stack[-1]

    def swap(self, original_widget):
        """
        Set a new widget as the top-most.
        """
        # First, add our new body to the stack.
        self.stack.append(original_widget)
        
        # Also add as the topmost.
        self.current_top.original_widget = original_widget
        
        # Update.
        self.update_view()
        
    #def popup(self, dialogue, align = "center", width = ('relative', 80), valign = "middle", height = ('relative', 80), **kwargs):
    def popup(self, dialogue, align = "center", width = 50, valign = "middle", height = 10, **kwargs):
        """
        Set a new dialogue popup as the top-most widget.
        
        In addition to the widget to show, this function takes the same arguments as the constructor for urwid Overlay widgets.
        
        :param dialogue: The widget to display as a popup.
        """
        # First, get an overlay widget to display both our top-most normal widget and our popup.
        # We will use another placeholder widget for the underneath widget, so we can swap it out.
        overlay = urwid.Overlay(dialogue, self.current_top, align = align, width = width, valign = valign, height = height, **kwargs)
        
        # Add the overlay to our dialogue stack.
        self.dialogue_stack.append(overlay)
        
        # Update.
        self.update_view()

    def swap_into_window(self, original_widget, cancel_callback = None, submit_callback = None):
        """
        Wrap a widget buttons and then set it as the top-most widget.
        
        :param cancel_callback: A function to call when the cancel button is pressed.
        :param submit_callback: A function to call when the submit button is pressed.
        """
        # If our widget is a View and a submit_callback has not been given, we can use a default.
        if isinstance(original_widget, View) and submit_callback is None:
            submit_callback = original_widget.submit
        
        # First decide which kind of wrapper to use.
        # If we have options, use a Confirm_settings_cancel.
        if isinstance(original_widget, View) and original_widget.has_settings:
            window = Confirm_settings_cancel(original_widget, top = self, settings_editor = original_widget.get_settings_editor(), cancel_callback = cancel_callback, submit_callback = submit_callback)
            
        elif submit_callback is not None:
            window = Confirm_or_cancel(original_widget, top = self, cancel_callback = cancel_callback, submit_callback = submit_callback)
            
        else:
            window = Confirm(original_widget, top = self, submit_callback = cancel_callback)
        
        self.swap(window)

    def back(self):
        """
        Remove the current top-most widget and set the last in its place.
        
        If there are no more widgets to go back to, urwid.ExitMainLoop will be raised.
        
        :param number: The number of times to go back.
        :return: The widget just removed, for convenience.
        """
        if len(self.dialogue_stack) > 0:
            # Close the top-most dialogue.
            self.close_popup()
            
        else:
            # Close the top-most real widget.
            self.close_widget()
            
    def close_widget(self):
        """
        Close the top-most currently viewable non-popup.
        """
        try:
            # Remove the top-most real widget.
            self.stack.pop()
            
            # Update our current_top.
            self.current_top.original_widget = self.stack[-1]
            
            # Update.
            self.update_view()
            
        except IndexError:
            # The stack is empty.
            raise urwid.ExitMainLoop()
            
    def close_popup(self, popup = None):
        """
        Close the top-most currently viewable popup.
        """
        # If we've been given a widget to close, close that one.
        if popup is not None:
            try:
                self.dialogue_stack.pop([True if overlay.top_w == popup else False for overlay in self.dialogue_stack].index(True))
                
            except IndexError:
                # Give a slightly more descriptive error.
                raise IndexError("close_popup() called but widget '{}' is not currently visible".format(popup))
            
        else:
            try:
                # Remove the top-most dialogue widget.
                self.dialogue_stack.pop()
                
            except IndexError:
                # Give a slightly more descriptive error.
                raise IndexError("close_popup() called but no popups are currently visible")
        
        # Update our view.
        self.update_view()
        
    def keypress(self, size, key):
        """
        Handler for keypress events.
        """
        # Allow children to intercept first.
        key = super().keypress(size, key)
        
        if key == "esc":
            if len(self.dialogue_stack) > 0:
                self.close_popup()
                
            else:
                self.back()
        
        else:
            return key


class View(WidgetWrap, Configurable):
    """
    A widget designed to be shown inside a swapping window.
    """
    
    def __init__(self, body, title = "", border = True, focusable = True):
        self._settings_editor = None
        inner_cls = Section if border else Sub_section
        self.inner_body = inner_cls(body, title, focusable)
        
        WidgetWrap.__init__(self, self.inner_body)
        Configurable.__init__(self, True)
        
    def submit(self):
        """
        A method to call when this View's submit button is pressed.
        """
        raise NotImplementedError("This view does not have a default submit() defined.")
    
    @property
    def has_settings(self):
        """
        Does this View have any configurable options set on it?
        """
        return len(self.OPTIONS) != 0
    
    def on_settings_change(self):
        """
        A method that will be called when settings have been changed.
        """
        # This default implementation does nothing.
    
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
        #return Settings_editor(View_browser(self), "Settings for {}".format(self.original_widget._title))
        return Settings_editor(View_browser(self), "Settings")
    