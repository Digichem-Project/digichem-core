# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.base import Section
from silico.interface.urwid.wrapper import Confirm_or_cancel, Confirm


class Dialogue(urwid.Overlay):
    """
    A confirmation box, designed to be used as a pop-up. Presents the user with a message and buttons.
    """
    
                
    def get_inner(self):
        """
        Get the inner widget we'll use for display.
        """
        return Section(
            urwid.ListBox(urwid.SimpleListWalker([urwid.Text(self.message, wrap="any")])),
            title = self.title,
            focusable = False
        )
        
class Confirm_or_cancel_dialogue(Dialogue):
    """
    A dialogue box with two buttons.
    """
    
    def __init__(self, title, message, top, cancel_callback = None, submit_callback = None):
        """
        Constructor for Confirm dialogue boxes.
        
        :param title: Text markup to use as the title.
        :param message: Text markup to use as the message body.
        :param top: The top widget used for display.
        :param cancel_callback: Function to call when the cancel button is clicked.
        :param submit_callback: Function to call when the submit button is clicked.
        """
        self.title = title
        self.message = message
        
        body = urwid.AttrMap(Confirm_or_cancel(self.get_inner(), top, cancel_callback = cancel_callback, submit_callback = submit_callback), {"body": "dialogue", "section": "section--dialogue"})
        
        # Call parent.
        super().__init__(
            body,
            top.original_widget,
            align = "center",
            width = ('relative', 80),
            valign = "middle",
            height = ('relative', 80),
        )


class Confirm_dialogue(Dialogue):
    """
    A dialogue box with one button.
    """
    
    def __init__(self, title, message, top, error = False, submit_callback = None):
        """
        Constructor for Confirm dialogue boxes.
        
        :param title: Text markup to use as the title.
        :param message: Text markup to use as the message body.
        :param top: The top widget used for display.
        :param error: Whether to use alternative formatting to indicate an error.
        :param submit_callback: Function to call when the button is clicked.
        """
        self.title = title
        self.message = message
        
        # Decide on formatting.
        if error:
            # Error formatting.
            attrs = {"body": "dialogue", "section": "section--dialogue"}
        
        else:
            # Normal formatting.
            attrs = {"body": "dialogue", "section": "section--dialogue--error"}
            
        
        body = urwid.AttrMap(Confirm(self.get_inner(), top, submit_callback = submit_callback), attrs)
                
        # Call parent.
        super().__init__(
             body,
            top.original_widget,
            align = "center",
            width = ('relative', 80),
            valign = "middle",
            height = ('relative', 80),
        )