import urwid

from silico.interface.urwid.base import Window
from silico.interface.urwid.misc import Tab_pile
from .editors import Option_editor, Bool_editor, Text_editor
from silico.interface.urwid.editors import Options_editor

class Configurable_editor(Options_editor):
    """
    Urwid widget used to edit a configurable.
    """
    
    def __init__(self, configurable, *, top):
        """
        Constructor for Configurable_editor objects.
        """
        self.top = top
        
        # Save our configurable.
        self.configurable = configurable
                
        # Build our editing widgets.
        self.fields = []
        #editors = Options_editor.get_widgets(self, self.configurable, self.configurable)
        editors = super().get_widgets(self.configurable, self.configurable)
            
        # Remove the last widget (unneeded Divider).
        editors.pop()
        
        edit_list = urwid.ListBox(
            urwid.SimpleFocusListWalker(editors)
        )
        
        # Buttons.
        cancel = urwid.Button("Cancel", self.cancel)
        confirm = urwid.Button("Confirm", self.confirm)
        
        # Call our parent constructor to get our box.
        Tab_pile.__init__(self, [
            ('weight', 1, urwid.LineBox(edit_list)),
            ('pack', urwid.Columns([
                urwid.AttrMap(urwid.Padding(cancel, align = "center", width = 10), 'bad_button'),
                urwid.AttrMap(urwid.Padding(confirm, align = "center", width = 11), 'good_button')
            ]))
        ])
        
    @property
    def title(self):
        """
        """
        return "Editing {} '{}' with ID: {}".format(self.configurable.TYPE, self.configurable.CLASS, self.configurable.ID)
    
    def window(self):
        """
        Get a window containing this editor.
        """
        return Window(
            urwid.AttrMap(self, 'body'),
            #title = "Silico Configurable Editer"
            title = self.title
        )
    
    def cancel(self, *args, **kwargs):
        """
        Discard any changes and close the window.
        """
        self.top.back()
        
    
    def confirm(self, *args, **kwargs):
        """
        Confirm any changes and close the window.
        """
        self.save()
        self.top.back()
        
        
        