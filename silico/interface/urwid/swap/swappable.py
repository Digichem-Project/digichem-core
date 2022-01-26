# General imports.
import urwid

# Silico imports.
from silico.config.configurable.base import Configurable
from silico.interface.urwid.layout import Sub_pane, Pane

class Swappable(urwid.WidgetWrap, Configurable):
    """
    A widget designed to be shown inside a swapping window.
    """
    
    def __init__(self, top, body, title = "", border = True, focusable = True):
        self.top = top
        self._settings_editor = None
        self.inner = body
        self.title = title
        self.border = border
        self.focusable = focusable
        
        urwid.WidgetWrap.__init__(self, self._get_body())
        Configurable.__init__(self, True)
    
    def _get_body(self):
        """
        Get the widget we'll use for display.
        """
        # TODO: That the body is wrapped in a Pane object is weird and odd and should probably be removed. 
        inner_cls = Pane if self.border else Sub_pane
        inner_body = inner_cls(self.inner, self.title, self.focusable)
        return inner_body
        
    def submit(self):
        """
        A method to call when this Swappable's submit button is pressed.
        """
        raise NotImplementedError("This view does not have a default submit() defined.")
    
    @property
    def has_settings(self):
        """
        Does this Swappable have any configurable options set on it?
        """
        return len(self.OPTIONS) != 0
    
    def on_settings_change(self):
        """
        A method that will be called when settings have been changed.
        """
        # This default implementation does nothing.
