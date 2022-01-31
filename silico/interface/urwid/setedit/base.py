# Base code for Settings_browser widgets.

# General imports.
import urwid
from urwid.listbox import SimpleFocusListWalker

# Silico imports.
from silico.interface.urwid.setedit.widget import Setedit_widget
from silico.interface.urwid.setedit.common import Setedit_widget_parent_mixin
from silico.interface.urwid.pages import Pages
from silico.logging.base import get_logger
from silico.interface.urwid.layout import Pane


class Setedit():
    """
    Logical storage for setting editors.
    """
    
    def __init__(self, top, title, starting_value, vtype, help = None, choices = None):
        """
        
        :param top: Top-most widget being used for display.
        :param title: The title/name of this setedit.
        :param starting_value: The default value.
        :param vtype: The type of value.
        :param help: Optional help text to display.
        :param choices: An optional list of choices which restricts the values this setedit can take.
        """
        self.top = top
        self.title = title
        self.previous_value = starting_value
        self.vtype = vtype
        self.help = help
        self.choices = choices
        self._widget = None
        
    def confirm(self):
        """
        Confirm the changes made to this Setedit, so that future rollbacks will return to the current value.
        """
        self.previous_value = self.get_widget().value
        
    def get_children(self):
        """
        Get child setedits of this setedit.
        
        This method should only be called if self.previous_value contains values of the form (previous_value, vtype, help). This is likely only to be the case if self.vtype == "Options".
        """
        return [type(self)(name, previous_value, vtype, help) for name, (previous_value, vtype, help) in self.previous_value.items()]
    
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
        
    @classmethod
    def vtype_from_configurable_option(self, option):
        """
        Get the type of value expected by a configurable option.
        
        :param option: The configurable option.
        :returns: The value type (a string).
        """
        if hasattr(option, "OPTIONS"):
            return "Options"
        
        elif len(option.choices) > 0:
            return "choices"
        
        elif option.list_type is not None:
            return "list"
        
        else: 
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
            # This doesn't work if option is a sub option.
            return option.__get__(owning_obj)


class Setedit_browser(urwid.ListBox, Setedit_widget_parent_mixin):
    """
    A widget that permits viewing and editing lists of options.
    """
    
    def __init__(self, setedits, on_change_callback = None):
        """
        Constructor for Setedit_browser objects.
        
        :param setedits: A list of  Setedit objects to browse through.
        :param on_change_callback: A function to call (with no arguments) when new settings are saved.
        """
        super().__init__(SimpleFocusListWalker(self.load_child_widgets(setedits)))
        self.on_change_callback = on_change_callback
        
    @property
    def child_widgets(self):
        """
        A shortcut to the property where the child widgets are stored.
        
        This must be defined by the inheriting class because it depends on the type of the inheriting class, for example self.body for listboxes.
        """
        return self.body
            
    def save(self):
        """
        Save the changes made.
        """
        for child_widget in self.child_setedit_widgets.values():
            child_widget.setedit.confirm()
            
    def confirm_callback(self):
        """
        Method called when settings have been changed.
        """
        try:
            retval = self.save()
            
        except Exception:
            get_logger().error("Failed to save changes", exc_info = True)
            return False
        
        if self.on_change_callback is not None:
            self.on_change_callback()
            
        return retval


# TOOD: Should inherit from some ABC?
class Paginated_settings_browser(Pages):
    """
    A widget for editing multiple pages of settings.
    """
    
    def __init__(self, *args, on_change_callback = None, **kwargs):
        """
        Constructor for Paginated_settings_browser objects.
        
        :param pages: An (ordered) dict of Setedit_browser objects (or similar), where each item is the browser object and each key is the name/title of that browser.
        :param title: The title of this paginated browser.
        """
        super().__init__(*args, **kwargs)
        self.on_change_callback = on_change_callback
        
    def confirm_callback(self):
        """
        Method called when settings have been changed.
        """
        try:
            # TODO: Need to be smarter about this, at the very least we should swap to the page where the bad option is.
            self.save()
        
        except Exception:
            get_logger().error("Failed to save changes", exc_info = True)
            return False
    
    def save(self):
        """
        Save any changes made.
        """
        # Save each of the individual browsers that we encapsulate.
        for page in self.pages.values():
            page.base_widget.save()
            
        self.validate()
        self.on_change_callback()
            
    def validate(self):
        """
        Validate any recent changes made.
        
        If any of the changes are invalid, an exception will be raised.
        """
        # This default implementation does nothing.
    
    def discard(self):
        """
        Discard any changes made.
        """
        for page in self.pages.values():
            page.base_widget.discard()
