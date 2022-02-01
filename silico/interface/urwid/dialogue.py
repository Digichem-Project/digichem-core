# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.wrapper import Confirm_or_cancel, Confirm
from silico.interface.urwid.layout import Pane


class Dialogue_mixin():
    """
    A mixin class for dialogue objects.
    """
    
    def setup(self, title, body, error = False):
        """
        Constructor for Dialogues, this should be called by real classes.
        This is not called __init__ to avoid problems with multiple inheritance.
        """
        self.dialogue_body = body
        self.title = title
        self.error = error
    
    def get_body(self):
        """
        Get the inner body of this confirmation dialogue.
        """
        # Decide on formatting.
        if self.error:
            # Error formatting.
            attrs = {"body": "pane--dialogue--error", "pane": "pane--dialogue--error"}
        
        else:
            # Normal formatting.
            attrs = {"body": "dialogue", "pane": "pane--dialogue"}
            
        return urwid.AttrMap(Pane(
            self.dialogue_body,
            title = self.title,
            focusable = False
        ), attrs)
        
    def back(self):
        """
        Function called to close this wrapper.
        """
        self.top.close_popup(self)
    
    def convert_callback(self, callback):
        """
        Convert certain values which can be interpreted as a function into a real function.
        
        The interpretation is as follows:
            - int: A number of steps to call back() on the top widget.
            - otherwise: A function to call.
            
        :returns: The possibly new callback function.
        """
        def back_callback():
            # Close ourself first.
            self.back()
            
            # Then close any widgets.
            for iter in range(0, callback -1):
                self.top.close_widget()
            
            # Stop going back any further.
            return False
        
        if isinstance(callback, int):
            return back_callback
        
        else:
            return callback
        

class Widget_dialogue(Confirm_or_cancel, Dialogue_mixin):
    """
    Class for dialogues that display other widgets.
    """
    
    def __init__(self, title, widget, top, error = False, cancel_callback = None, submit_callback = None):
        """
        Constructor for Confirm dialogue boxes.
        
        :param title: Text markup to use as the title.
        :param widget: The widget to wrap around, this should probably be inside a list box or similar...
        :param top: The top widget used for display.
        :param error: Whether to use alternative formatting to indicate an error.
        :param cancel_callback: Function to call when the cancel button is clicked.
        :param submit_callback: Function to call when the submit button is clicked.
        """
        Dialogue_mixin.setup(self, title, widget, error)
        super().__init__(self.get_body(), top, cancel_callback, submit_callback)

    def back(self, *args, **kwargs):
        return Dialogue_mixin.back(self, *args, **kwargs)
    
    def convert_callback(self, *args, **kwargs):
        return Dialogue_mixin.convert_callback(self, *args, **kwargs)
    
    
        
class Text_dialogue(Dialogue_mixin):
    """
    Mixin class for dialogues that display text
    """
    
    def setup(self, title, message, error = False):
        """
        Constructor for Dialogues, this should be called by real classes.
        """
        super().setup(title, urwid.ListBox(urwid.SimpleListWalker([urwid.Text(('body', message))])), error)


class Confirm_or_cancel_dialogue(Confirm_or_cancel, Text_dialogue):
    """
    A dialogue box with two buttons.
    """
    
    def __init__(self, title, message, top, error = False, cancel_callback = None, submit_callback = None):
        """
        Constructor for Confirm dialogue boxes.
        
        :param title: Text markup to use as the title.
        :param message: Text markup to use as the message body.
        :param top: The top widget used for display.
        :param error: Whether to use alternative formatting to indicate an error.
        :param cancel_callback: Function to call when the cancel button is clicked.
        :param submit_callback: Function to call when the submit button is clicked.
        """
        Text_dialogue.setup(self, title, message, error)
        super().__init__(self.get_body(), top, cancel_callback, submit_callback)

    def back(self, *args, **kwargs):
        return Text_dialogue.back(self, *args, **kwargs)
    
    def convert_callback(self, *args, **kwargs):
        return Text_dialogue.convert_callback(self, *args, **kwargs)


class Confirm_dialogue(Confirm, Text_dialogue):
    """
    A dialogue box with one button.
    """
    
    def __init__(self, title, message, top, error = False, submit_callback = None):
        """
        Constructor for confirmation dialogues.
        
        :param title: Text markup to use as the title.
        :param message: Text markup to use as the message body.
        :param top: The top widget used for display.
        :param error: Whether to use alternative formatting to indicate an error.
        :param submit_callback: Function to call when the button is clicked.
        """
        Text_dialogue.setup(self, title, message, error)
        super().__init__(self.get_body(), top, submit_callback = submit_callback)

    def back(self, *args, **kwargs):
        return Text_dialogue.back(self, *args, **kwargs)
    
    def convert_callback(self, *args, **kwargs):
        return Text_dialogue.convert_callback(self, *args, **kwargs)


class Output_dialogue(Confirm, Dialogue_mixin):
    """
    A widget that pops up to display stdout and stderr.
    """
    
    def __init__(self, top, max_length = 500):
        """
        Constructor for Confirm dialogue boxes.
        
        :param top: The top widget used for display.
        :param max_length: The maximum number of items to display. If this value is exceeded, older values will be deleted.
        """
        self.max_length = max_length
        self.list = urwid.ListBox(urwid.SimpleListWalker([]))
        
        Dialogue_mixin.setup(self, "Program Output", self.list, error = False)
        super().__init__(self.get_body(), top)

    def back(self, *args, **kwargs):
        # Empty our text so we always display something new.
        self.list.body.clear()
        return Dialogue_mixin.back(self, *args, **kwargs)
    
    def convert_callback(self, *args, **kwargs):
        return Dialogue_mixin.convert_callback(self, *args, **kwargs)
        
    def output(self, message, error = False):
        """
        Add output text to this widget.
        
        :param message: The output text to add.
        :param error: Whether to show alternative formatting to indicate an error.
        """
        # Decide on formatting.
        if error:
            # Error formatting.
            attr = "pane--dialogue--error"
        
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
        return False
        