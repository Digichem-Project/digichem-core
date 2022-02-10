# Code for editing configurables.

# General imports.

# Silico imports.
from silico.interface.urwid.setedit.base import Setedit_browser, Setedit,\
    Paginated_settings_browser
from silico.exception.base import Silico_exception
from silico.interface.urwid.setedit.widget import Solo_sub_editor
from silico.interface.urwid.layout import Pane
from silico.interface.urwid.dialogue import Confirm_or_cancel_dialogue
import urwid
from silico.exception.configurable import Missing_option_exception


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
        try:
            #self.previous_value = self.option.get_from_dict(self.owning_obj, self.dict_obj)
            self.previous_value = self.option.dump(self.owning_obj, self.dict_obj)
        
        except Missing_option_exception:
            self.previous_value = None
            
    @property
    def required(self):
        """
        Whether this option must have a value set.
        """
        return self.option.required
        
    def reset(self):
        """
        Reset the value of this option back to its default (if it has one).
        """
        try:
            default = self.option.default(self.owning_obj)
            
        except AttributeError:
            # No default, do nothing.
            return
        
        self.get_widget().value = default
        
    def refresh(self):
        """
        Refresh the current edit value of the edit widgets of this setedit (in case the underlying configurable option value has changed).
        """
        self.previous_value = self.option.get_from_dict(self.owning_obj, self.dict_obj)
        self.get_widget().discard()
        
    def confirm(self):
        """
        Confirm the changes made to this Setedit, so that future rollbacks will return to the current value.
        """
        # First update the value of our widget from the value of the configurable option (because it could have changed, for example from type conversion etc).
        #self.get_widget().value = self.option.get_from_dict(self.owning_obj, self.dict_obj)
        self.get_widget().value = self.option.dump(self.owning_obj, self.dict_obj)
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
    
    def reset(self):
        """
        Reset the value of this option back to its default (if it has one).
        """
        for child in self.get_children():
            child.reset()
    
    def refresh(self):
        """
        Refresh the current edit value of the edit widgets of this setedit (in case the underlying configurable option value has changed).
        """
        for child in self.get_children():
            child.refresh()
        
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
    

class Options_solo_setedit(Options_setedit):
    """
    An alternative setedit class for editing configurable Options (plural) settings.
    
    This setedit acts like a browser, allowing a single set of sub-settings to be edited. It is best used in conjunction with a paginated settings browser.
    """
    
    def __init__(self, top, owning_obj, dict_obj, option):
        super().__init__(top, owning_obj, dict_obj, option)
    
    def load_widget(self):
        """
        Load the widget we'll use for display.
        """
        return Solo_sub_editor(self)


class Configurable_browser(Setedit_browser):
    """
    A widget that permits viewing and editing lists of options.
    """
            
    def __init__(self, setedits, top, configurable, on_change_callback = None, can_reset = True):
        """
        Construct a Setedit_browser from a configurable object.
        
        :param setedits: A list of setedit objects to edit.
        :param top: Top-most widget being used for display.
        :param configurable: The configurable object which is being edited.
        :param on_change_callback: A function to call when settings are saved.
        """
        self.configurable = configurable
        
        # Get our list of edit widgets.
        widgets = self.load_child_widgets(setedits)
        
        if can_reset:
            # Add a button at the bottom which will reset our options.
            widgets.append(urwid.Divider())
            widgets.append(urwid.AttrMap(urwid.Button("Reset to defaults", lambda button: self.popup_reset_dialogue()), 'button--settings', 'button--settings--focused'))
        
        super().__init__(widgets, top, on_change_callback = on_change_callback)    
        
    def popup_reset_dialogue(self):
        """
        Show the popup dialogue which prompts the user whether they are sure they want to reset to defaults.
        """
        text = ""
        self.top.popup(Confirm_or_cancel_dialogue("Reset to defaults", "Are you sure you want to reset all options in this list?\n\nNote: After reseting, changes will not be saved until the 'save' button is pressed.", self.top, submit_callback = self.reset))
        
    def reset(self):
        """
        Reset each of the settings in this browser to their default (if they have one).
        """
        for widget in self.child_setedit_widgets.values():
            widget.setedit.reset()
        
    @classmethod
    def from_configurable(self, top, configurable, on_change_callback = None, can_reset = True):
        """
        Alternative constructor that creates a Configurable_browser object with all configurable options of a Configurable a its body.
        
        :param top: Top-most widget being used for display.
        :param configurable: A configurable object to construct from.
        :param on_change_callback: A function to call when settings are saved.
        """
        setedits = [Option_setedit.from_configurable_option(top, configurable, option) for option in configurable.OPTIONS.values() if option.no_edit is False]
        return self(setedits, top, configurable, on_change_callback = on_change_callback, can_reset = can_reset)
        
    def refresh(self):
        """
        Refresh the current edit value of each of the child widgets of this browser (in case the underlying configurable option value has changed).
        """
        for child_widget in self.child_setedit_widgets.values():
            child_widget.setedit.refresh()

    def save(self, validate = True):
        """
        Update the configurable we are editing with the current values.
        
        :param validate: Whether to validate the changes made.
        """
        options = self.configurable.OPTIONS
        child_widgets = list(self.child_setedit_widgets.values())
        
        for child_widget in child_widgets:
            value = child_widget.value
            options[child_widget.setedit.title].__set__(self.configurable, value)
            
        if validate:
            self.validate_setedits()
            # TODO: This should probably not be called from here...
            self.confirm()
        
            
    def validate_setedits(self):
        """
        Check the currently set values of our configurable are valid.
        
        Note that this function validates all configurable options of the configurable; not just the ones that have been changed by this editor.
        """
        child_widgets = list(self.child_setedit_widgets.values())
        # Check the values we just set are allowed.
        try:
            self.configurable.validate()
        
        except Exception as e:
            raise Silico_exception("An option has an invalid value") from e
        
        ## All good, confirm changes.
        #for child_widget in child_widgets:
        #    child_widget.setedit.confirm()
            
    def confirm(self):
        """
        Confirm the currently set values, so future resets() will rollback to the state as it is now.
        """
        child_widgets = list(self.child_setedit_widgets.values())
        # All good, confirm changes.
        for child_widget in child_widgets:
            child_widget.setedit.confirm()
        


#########################
# Convenience functions #
#########################
def make_paginated_configurable_browser(configurable, top, general_page_name = "general", on_change_callback = None, page_selector_title = None):
    """
    Create a paginated setting browser that has one page for each sub_option of a configurable object, plus one additional page for all top-level options.
    
    :param configurable: A configurable to create a settings browser for.
    :param top: A top-level widget to use for switching the current widget.
    :param general_page_name: The name to use for the general settings page where top-level settings are located.
    :param on_change_callback: A function to call when settings are changed.
    :param page_selector_title: The title to use for the widget which allows changing page.
    """
    # First, split our options based on those that do and do not have children.
    options_without_children = []
    options_with_children = []

    for name, option in configurable.OPTIONS.items():
        if not option.no_edit:
            # Split based on children.
            if option.num_child_options == 0:
                options_without_children.append(option)
            
            else:
                options_with_children.append((name, option))
    
    
    # Next, make a page for top-level (no children) settings.
    pages = {general_page_name: Pane(Configurable_browser([Option_setedit(top, configurable, configurable._configurable_options, option) for option in options_without_children], top, configurable), general_page_name)}
    
    # Next, add one page for each sub_option.
    pages.update({name: Pane(Configurable_browser([Options_solo_setedit(top, configurable, configurable._configurable_options, option)], top, configurable), name) for name, option in options_with_children})
    
    browser = Paginated_settings_browser(pages, on_change_callback = on_change_callback, title = page_selector_title)
    
    # Done.
    return browser


def make_settings_page_from_configurable_option(top, configurable, configurable_options, option):
    """
    Create a single page that might be used in a paginated_configurable_browser to edit the settings of an option from a configurable.
    
    :param top: A top-level widget to use for switching the current widget.
    :param configurable: The configurable that contains the option.
    :param configurable_options: The dict where the value of option is stored.
    :param option: The configurable option to edit.
    :returns: A tuple of (name, page).
    """
    return (option.name, Pane(Configurable_browser([Options_solo_setedit(top, configurable, configurable_options, option)], top, configurable), option.name))
        
        