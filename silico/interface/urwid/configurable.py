# Code for editing configurables.

# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.setedit.base import Settings_editor, Setedit_browser
from silico.interface.urwid.pages import Pages


class Configurable_editor(Settings_editor):
    """
    A high-level widget for changing the settings of a configurable.
    """
    
    def __init__(self, configurable):
        self.configurable = configurable
        
        super().__init__(Setedit_browser.from_configurable(configurable), "Editing {} with ID: {}".format(configurable.TYPE, configurable.index()))
    

class Method_editor(Pages):
    """
    A widget for editing different parts of a method.
    """
    
    def __init__(self, destination, program, calculation):
        """
        """
        # We show in reverse order, because typically the calculation is what the user wants to change anyway.
        super().__init__([
            ("Calculation", Configurable_editor(calculation)),
            ("Program", Configurable_editor(program)),
            ("Destination", Configurable_editor(destination))
        ])