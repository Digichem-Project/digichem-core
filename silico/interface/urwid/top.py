# Silico imports.
from silico.interface.urwid.wrapper import Confirm_or_cancel, Confirm_settings_cancel, Confirm
from silico.interface.urwid.section import Section, Sub_section
from silico.config.configurable.base import Configurable
from silico.interface.urwid.setedit.base import Settings_editor
from silico.interface.urwid.setedit.view import View_browser

# General imports.
import urwid
from urwid.widget import WidgetWrap
from silico.interface.urwid.dialogue import Dialogue
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
        self.output_widget = Output_widget(self)
        
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
        # First, add our message to our widget.
        self.output_widget.output(message, error)
        
        # Then swap to it if it's not already the top.
        if len(self.dialogue_stack) == 0 or self.dialogue_stack[-1].top_w != self.output_widget:
            self.popup(self.output_widget)
            
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
        
    def popup(self, dialogue, align = "center", width = ('relative', 80), valign = "middle", height = ('relative', 80)):
        """
        Set a new dialogue popup as the top-most widget.
        """        
        # First, get an overlay widget to display both our top-most normal widget and our popup.
        # We will use another placeholder widget for the underneath widget, so we can swap it out.
        overlay = urwid.Overlay(dialogue, self.current_top, align = align, width = width, valign = valign, height = height)
        
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

    def back(self, number = 1):
        """
        Remove the current top-most widget and set the last in its place.
        
        If there are no more widgets to go back to, urwid.ExitMainLoop will be raised.
        
        :param number: The number of times to go back.
        :param force: Whether to force switching back even if the current top is the output_widget (which is normally sticky).
        :return: The widget just removed, for convenience.
        """
        for iter in range(0, number):
            try:
                # Remove the top-most real widget.
                self.stack.pop()
                
                # Update.
                self.update_view()
                
            except IndexError:
                # The stack is empty.
                raise urwid.ExitMainLoop()
            
    def close_popup(self):
        """
        Close the top-most currently viewable popup.
        """
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
        
        return key


class Output_widget(Dialogue):
    """
    A widget that pops up to display stdout and stderr.
    """
    
    def __init__(self, top, max_length = 100):
        """
        Constructor for Confirm dialogue boxes.
        
        :param top: The top widget used for display.
        """
        self.max_length = max_length            
        super().__init__(title = "Program Output", top = top)
        
    def output(self, message, error = False):
        """
        Add output text to this widget.
        
        :param message: The output text to add.
        :param error: Whether to show alternative formatting to indicate an error.
        """
        # Decide on formatting.
        if error:
            # Error formatting.
            attr = "section--dialogue--error"
        
        else:
            # Normal formatting.
            attr = "dialogue"
        
        # Add text.
        for sub_message in message.splitlines():
            new_text = urwid.Text((attr, sub_message))
            self.list.body.append(new_text)
        
        # If we've gone over max, remove some.
        while (len(self.list.body) > self.max_length):
            self.list.body.pop(0)
        
        # Scroll to bottom.
        self.list.body.focus = len(self.list.body) -1
        
    def close(self):
        """
        """
        self.topmost.close_popup()
        return False
        
    def get_body(self):
        """
        Get the outer widget we'll use for display.
        """
        return Confirm(self.get_inner(), self.topmost, submit_callback = self.close)


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
    