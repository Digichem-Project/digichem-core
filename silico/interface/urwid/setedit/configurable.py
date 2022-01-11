# Code for editing configurables.

# General imports.

# Silico imports.
import silico.logging
from silico.interface.urwid.setedit.base import Setedit_browser, Option_setedit


class Configurable_browser(Setedit_browser):
    """
    A widget that permits viewing and editing lists of options.
    """
            
    def __init__(self, top, configurable, on_change_callback = None):
        """
        Construct a Setedit_browser from a configurable object.
        
        :param top: Top-most widget being used for display.
        :param configurable: A configurable object to construct from.
        :param on_change_callback: A function to call when settings are saved.
        """
        self.configurable = configurable
        setedits = [Option_setedit.from_configurable_option(top, configurable, option) for option in configurable.OPTIONS.values() if option.no_edit is False]
        super().__init__(setedits, on_change_callback = on_change_callback)
        
    def save(self):
        """
        Update the configurable we are editing with the current values.
        """
        options = self.configurable.OPTIONS
        child_widgets = list(self.child_setedit_widgets.values())
        
        for child_widget in child_widgets:
            value = child_widget.value
            options[child_widget.setedit.title].__set__(self.configurable, value)
            # Also update the 'default' value of the setedit.
            #setedit.default_value = setedit.value_from_configurable_option(self.configurable, options[setedit.title])
        
        # Check the values we just set are allowed.
        try:
            self.configurable.validate()
        
        except Exception:
            silico.logging.get_logger().error("Failed to save changes; an option has an invalid value", exc_info = True)
            return False
        
        # All good, confirm changes.
        for child_widget in child_widgets:
            child_widget.setedit.confirm()
        
        