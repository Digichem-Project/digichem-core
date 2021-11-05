# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.section import Section
from silico.interface.urwid.wrapper import Confirm_or_cancel, Confirm


class Dialogue(urwid.AttrMap):
    """
    A confirmation box, designed to be used as a pop-up. Presents the user with a message and buttons.
    """
    
    def __init__(self, title, top, initial = None, error = False):
        """
        General constructor for Dialogue objects.
        
        :param title: The title of this dialogue.
        :param top: The top-most widget to use for display.
        :param initial: initial list of widgets to display in our body.
        """
        self.topmost = top
        self.title = title
        self.list = urwid.ListBox(urwid.SimpleListWalker(initial if initial != None else []))
        self.error = error
        
        # Decide on formatting.
        if self.error:
            # Error formatting.
            attrs = {"body": "dialogue", "section": "section--dialogue--error"}
        
        else:
            # Normal formatting.
            attrs = {"body": "dialogue", "section": "section--dialogue"}
        
        # Call parent.
        super().__init__(self.get_body(), attrs)
                
    def get_inner(self):
        """
        Get the inner widget we'll use for display.
        """
        
        return Section(
            self.list,
            title = self.title,
            focusable = False
        )
        
    def get_body(self):
        """
        Get the outer widget we'll use for display.
        """
        return self.get_inner()


class Confirm_or_cancel_dialogue(Dialogue):
    """
    A dialogue box with two buttons.
    """
    
    def __init__(self, title, top, message = "", error = False, cancel_callback = None, submit_callback = None):
        """
        Constructor for Confirm dialogue boxes.
        
        :param title: Text markup to use as the title.
        :param message: Text markup to use as the message body.
        :param top: The top widget used for display.
        :param error: Whether to use alternative formatting to indicate an error.
        :param cancel_callback: Function to call when the cancel button is clicked.
        :param submit_callback: Function to call when the submit button is clicked.
        """
        self.cancel_callback = cancel_callback
        self.submit_callback = submit_callback
        
        super().__init__(title, top, [urwid.Text(message, wrap="any")], error)
        
    def get_body(self):
        """
        Get the outer widget we'll use for display.
        """
        return Confirm_or_cancel(self.get_inner(), self.topmost, cancel_callback = self.cancel_callback, submit_callback = self.submit_callback)


class Confirm_dialogue(Dialogue):
    """
    A dialogue box with one button.
    """
    
    def __init__(self, title, top, message = "",  error = False, submit_callback = None):
        """
        Constructor for Confirm dialogue boxes.
        
        :param title: Text markup to use as the title.
        :param message: Text markup to use as the message body.
        :param top: The top widget used for display.
        :param error: Whether to use alternative formatting to indicate an error.
        :param submit_callback: Function to call when the button is clicked.
        """
        self.submit_callback = submit_callback
                        
        super().__init__(title, top, [urwid.Text(message, wrap="any")], error)
        
    def get_body(self):
        """
        Get the outer widget we'll use for display.
        """
        return Confirm(self.get_inner(), self.topmost, submit_callback = self.submit_callback)