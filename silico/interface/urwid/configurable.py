# Code for editing configurables.

# General imports.

# Silico imports.
from silico.interface.urwid.setedit.base import Settings_editor, Setedit_browser,\
    Setedit
from silico.interface.urwid.pages import Pages


class Configurable_browser(Setedit_browser):
    """
    A widget that permits viewing and editing lists of options.
    """
            
    def __init__(self, configurable):
        """
        Construct a Setedit_browser from a configurable object.
        
        :param configurable: A configurable object to construct from.
        """
        self.configurable = configurable
        setedits = [Setedit.from_configurable_option(configurable, option) for option in configurable.OPTIONS.values() if option.no_edit is False]
        super().__init__(setedits)
        
    def save(self):
        """
        Update the configurable we are editing with the current values.
        """
        options = self.configurable.OPTIONS
        
        for setedit in self.body:
            value = setedit.get_value()
            options[setedit.title].__set__(self.configurable, value)
            
            # Also update the 'default' value of the setedit.
            setedit.default_value = setedit.value_from_configurable_option(self.configurable, options[setedit.title])


class Configurable_editor(Settings_editor):
    """
    A high-level widget for changing the settings of a configurable.
    """
    
    def __init__(self, configurable, title):
        self.configurable = configurable
        super().__init__(Configurable_browser(configurable), title)
    

class Method_editor(Pages):
    """
    A widget for editing different parts of a method.
    """
    
    def __init__(self, destination, program, calculation):
        """
        """
        # We show in reverse order, because typically the calculation is what the user wants to change anyway.
        super().__init__([
            ("Calculation", Configurable_editor(calculation, self.get_title(calculation))),
            ("Program", Configurable_editor(program, self.get_title(program))),
            ("Destination", Configurable_editor(destination, self.get_title(destination)))
        ])
        
    def get_title(self, method_target):
        """
        Get the title to display for a given method_target.
        """
        return "Editing {} with ID: {}".format(method_target.TYPE, method_target.index())