# Base code for Settings_browser widgets.

# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.setedit.widget import Setedit_widget
from silico.interface.urwid.section import Section
from silico.interface.urwid.setedit.common import Setedit_widget_parent_mixin
from urwid.listbox import SimpleFocusListWalker


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


class Option_setedit(Setedit):
    """
    A setedit for editing a Configurable Option.
    """
    
    def __init__(self, top, owning_obj, dict_obj, option):
        """
        :param top: Top-most widget being used for display.
        :param owning_obj: The owning object on which the Option object is set as a class attribute.
        :param dict_obj: The dict in which the value of this Option is stored.
        :param value: The new value to set.
        """
        # Note that we don't call our parent's constructor
        self.top = top
        self.owning_obj = owning_obj
        self.dict_obj = dict_obj
        self.option = option
        self._widget = None
        
        self.title = option.name
        self.vtype = self.vtype_from_configurable_option(option)
        self.help = option.help
        self.choices = option.choices
        # This is the last value we saved, if we're asked to reset we'll roll back to this.
        self.previous_value = self.option.get_from_dict(self.owning_obj, self.dict_obj)
        
    @property
    def default_value(self):
        """
        The default value of this setedit. If an equivalent value is given by the user it will not be saved, instead the true default will be used.
        """
        # This can raise an attribute error.
        return self.option.default(self.owning_obj)
        
    def confirm(self):
        """
        Confirm the changes made to this Setedit, so that future rollbacks will return to the current value.
        """
        # First update the value of our widget from the value of the configurable option (because it could have changed, for example from type conversion etc).
        self.get_widget().value = self.option.get_from_dict(self.owning_obj, self.dict_obj)
        super().confirm()
    
    def get_children(self, reload = False):
        """
        Get a list of Option_setedit object that represent the options that are children of this option.
        """
        raise NotImplementedError("Option objects don't have children")
    
    @classmethod
    def from_configurable_option(self, top, owning_obj, option):
        """
        Construct a new Setedit object from a configurable option.
        
        :param top: The top-most widget being used for display.
        :param owning_obj: The owning object on which the Option object is set as a class attribute.
        :param option: The configurable option.
        :returns: The Setedit object.
        """
        if hasattr(option, "OPTIONS"):
            cls = Options_setedit
            
        else:
            cls = self
            
        return cls(top, owning_obj, owning_obj._configurable_options, option)


class Options_setedit(Option_setedit):
    """
    A setedit for editing a group of Configurable Options.
    """
    
    def __init__(self, top, owning_obj, dict_obj, option):
        super().__init__(top, owning_obj, dict_obj, option)
        self._children = None
        
    def confirm(self):
        """
        Confirm the changes made to this Setedit, so that future rollbacks will return to the current value.
        """
        for child in self.get_children():
            child.confirm()
    
    @property
    def previous_value(self):
        raise NotImplementedError("Options_setedit does not contain values")
    
    @previous_value.setter
    def previous_value(self, value):
        # Do nothing.
        pass
    
    def load_children(self):
        """
        Load the child setedits of this one.
        """
        mapping = self.option.get_from_dict(self.owning_obj, self.dict_obj)
        
        children = []
        
        for sub_option in self.option.OPTIONS.values():
            if hasattr(sub_option, "OPTIONS"):
                # This sub option has sub options of its own.
                children.append(Options_setedit(self.top, self.owning_obj, mapping, sub_option))
                
            else:
                children.append(Option_setedit(self.top, self.owning_obj, mapping, sub_option))
                
        return children
    
    def get_children(self, reload = False):
        """
        Get a list of Option_setedit object that represent the options that are children of this option.
        """
        if self._children is None or reload:
            self._children = self.load_children()
        
        return self._children


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
        
    def _save(self):
        """
        Wrapper for save. This method should not be modified in child classes.
        """
        retval = self.save()
        if retval != False and self.on_change_callback is not None:
            self.on_change_callback()
        return retval


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
        
        
        
    
    