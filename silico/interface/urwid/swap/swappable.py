# General imports.
import urwid

# Silico imports.
from silico.config.configurable.base import Configurable


class Swappable(urwid.WidgetWrap, Configurable):
    """
    A widget designed to be shown inside a swapping window.
    """
    
    def __init__(self, top, body):
        self.top = top
        self._settings_editor = None
        
        urwid.WidgetWrap.__init__(self, body)
        Configurable.__init__(self, True)
        
    def submit(self):
        """
        A method to call when this Swappable's submit button is pressed.
        """
        raise NotImplementedError("This Swappable object does not have a default submit() defined.")
    
    @property
    def settings_title(self):
        """
        A title to use for the settings widget that can be used to edit the configurable options of this Swappable.
        """
        return "Settings"
    
    @property
    def has_settings(self):
        """
        Does this Swappable have any configurable options set on it?
        """
        return len(self.OPTIONS) != 0
    
    @property
    def additional_option_pages(self):
        """
        A dict of additional 'pages' of options to edit.
        """
        return {}
    
    def on_settings_change(self):
        """
        A method that will be called when settings have been changed.
        """
        # This default implementation does nothing.