# Base code for Settings_browser widgets.

# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.setedit.widget import Setedit_widget
from silico.interface.urwid.base import Section


class Setedit():
    """
    Logical storage for setting editors.
    """
    
    def __init__(self, title, default_value, vtype, help = None):
        """
        
        :param title: The title/name of this setedit.
        :param default_value: The default value.
        :param vtype: The type of value.
        :param help: Optional help text to display.
        """
        self.title = title
        self.default_value = default_value
        self.vtype = vtype
        self.help = help
        self._widget = None
        
    def get_children(self):
        """
        Get child setedits of this setedit.
        
        This method should only be called if self.default_value contains values of the form (default_value, vtype, help). This is likely only to be the case if self.vtype == "Options".
        """
        return [type(self)(name, default_value, vtype, help) for name, (default_value, vtype, help) in self.default_value.items()]
    
    def get_widget(self, reload = False):
        """
        Get the widget we'll use for display.
        """
        if self._widget is None or reload:
            self._widget = self.load_widget()
            
        return self._widget
    
    def load_widget(self):
        """
        Load the widget we'll use for display.
        """
        return Setedit_widget.class_from_type(self.vtype)(self)
    
    def get_value(self):
        """
        Get the possibly edited value of this Setedit.
        """
        return self.get_widget().value
    
    @classmethod
    def vtype_from_configurable_option(self, option):
        """
        Get the type of value expected by a configurable option.
        
        :param option: The configurable option.
        :returns: The value type (a string).
        """
        if hasattr(option, "OPTIONS"):
            return "Options"
        
        try:
            return option.type.__name__
        
        except Exception:
            return "str"
    
    @classmethod
    def value_from_configurable_option(self, owning_obj, option):
        """
        Get a dictionary of values from a Configurable options object.
        
        :param owning_obj: The owning object on which the Option object is set as a class attribute.
        :param option: The configurable option.
        """
        # First, keep an eye out for nested options.
        if hasattr(option, "OPTIONS"):
            values = {}
            for sub_option in option.OPTIONS.values():
                values[sub_option.name] = (self.value_from_configurable_option(owning_obj, sub_option), self.vtype_from_configurable_option(sub_option), sub_option.help)    
                    
            return values
        
        else:
            # A normal option, just get the value.
            return option.__get__(owning_obj)
    
    @classmethod
    def from_configurable_option(self, owning_obj, option):
        """
        Construct a new Setedit object from a configurable option.
        
        :param owning_obj: The owning object on which the Option object is set as a class attribute.
        :param option: The configurable option.
        :returns: The Setedit object.
        """
        return self(option.name, self.value_from_configurable_option(owning_obj, option), self.vtype_from_configurable_option(option), option.help)    


class Setedit_walker(urwid.SimpleFocusListWalker):
    """
    ListWalker-compatible class for displaying Settings.

    Positions are Setedit objects, while actual display is handled by Setedit_widgets.
    """

    def __init__(self, data):
        """
        Constructor for Row_walker objects.
        
        :param data: List of data (Row_items) to display.
        """
        super().__init__(data)
            
    def get_focus(self):
        """
        Get the current focus.
        
        :returns: A tuple (widget, position) of the current focus.
        """
        try:
            focus = self.focus
            return self[focus].get_widget(), focus
        except (IndexError, KeyError, TypeError):
            return None, None
        
    def get_next(self, position):
        """
        Get the next.
        
        :returns: A tuple (widget, position) of position.
        """
        try:
            position = self.next_position(position)
            return self[position].get_widget(), position
        except (IndexError, KeyError):
            return None, None

    def get_prev(self, position):
        """
        Get the previous.
        
        :returns: A tuple (widget, position) of position.
        """
        try:
            position = self.prev_position(position)
            return self[position].get_widget(), position
        except (IndexError, KeyError):
            return None, None


class Setedit_browser(urwid.ListBox):
    """
    A widget that permits viewing and editing lists of options.
    """
    
    def __init__(self, setedits):
        super().__init__(Setedit_walker(setedits))
        
    @classmethod
    def from_configurable(self, configurable):
        """
        Construct a Setedit_browser from a configurable object.
        
        :param configurable: A configurable object to construct from.
        :returns: A setedit_browser object.
        """        
        setedits = [Setedit.from_configurable_option(configurable, option) for option in configurable.OPTIONS.values() if option.no_edit is False]
        return self(setedits)
    
    def discard(self):
        """
        Discard any changes made without saving.
        """
        for setedit in self.body:
            setedit.get_widget().reset()
            


class Settings_editor(Section):
    """
    A high-level widget for editing a list of settings.
    """
    
    def __init__(self, browser, title):
        """
        Constructor for Settings_editor objects.
        
        :param browser: A Setedit_browser to use as our body.
        :param title: The title to show.
        """
        self.browser = browser
        
        super().__init__(self.browser, title)
        
        
        
    
    