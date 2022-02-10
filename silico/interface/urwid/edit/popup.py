# Edit widgets that appear as buttons, then 'popup'/switch to a knew window to allow editing.

# General imports.
import urwid

# Silico imports.
import silico.interface.urwid.file.browser
import silico.interface.urwid.file.output
from silico.interface.urwid.dialogue import Widget_dialogue
from silico.interface.urwid.method.browser import Method_selector


class Popup_edit(urwid.Button):
    """
    ABC for editors that show a popup.
    """
    
    def __init__(self, top, initial = None):
        """
        Constructor for Popup_edit widgets.
        
        :param top: The top-most widget being used for display.
        """
        self.top = top
        self._popup = None
        self._value = None
        self.initial = initial
        
        super().__init__("", lambda button: self.open_popup())
        self.value = initial
    
    def open_popup(self):
        """
        Method called to show the popup the user can interact with.
        """
        self.top.popup(self.get_popup())
        
    def close_popup(self):
        """
        Method called to hide the popup.
        """
        self.top.close_popup(self.get_popup())
        
    def get_popup(self, reload = False):
        """
        Get the popup widget.
        """
        if self._popup is None or reload:
            self._popup = self.load_popup()
            
        return self._popup
        
    def load_popup(self):
        """
        Load/create the popup.
        """
        raise NotImplementedError("Implement in subclass")
    
    @property
    def value(self):
        return self._value
    
    @value.setter
    def value(self, val):
        self._value = val
        self.update_label()
    
    def reset(self):
        """
        Reset the current value back to the default value.
        """
        self._value = self.setedit.previous_value
        self.update_label()
        
    def update(self):
        """
        Update the value of this widget.
        """
        raise NotImplementedError("Implement in subclass")
        
    def update_label(self):
        """
        Update the label of this widget with the current value.
        """
        self.set_label(str(self.value) if self.value is not None else "")


class Choice(urwid.AttrMap):
    """
    Widget used to display a choice in a choices picker object.
    """
    
    body_attr = "body"
    focus_attr = "editable"
    
    def __init__(self, value, picker):
        """
        """
        self.button = urwid.Button("", lambda button: self.submit())
        super().__init__(self.button, self.body_attr, self.focus_attr)
        
        self.value = value
        self.picker = picker
        
    def submit(self):
        """
        Method called when our choice is chosen.
        """
        self.picker.edit.update()
        self.picker.edit.close_popup()
        
    @property
    def value(self):
        """
        The value of this choice.
        """
        return self._value
    
    @value.setter
    def value(self, value):
        """
        Change the value of this choice.
        """
        self._value = value
        self.button.set_label(str(value) if value is not None else "(None)")


class Choices_picker(urwid.ListBox):
    """
    Widget that allows picking from a number of choices.
    """
    
    def __init__(self, edit):
        """
        """
        super().__init__(urwid.SimpleFocusListWalker([]))
        
        # Keep our widget for later.
        self.edit = edit
        
        for choice in edit.choices:
            self.body.append(self.get_widget(choice))
            
        self.set_choice(self.edit.value)
        
    def get_widget(self, value):
        """
        Get one of the choices widget (a button) we will display.
        """
        return Choice(value, self)
    
    def set_choice(self, value):
        """
        Set a choice as focus.
        
        :param value: The value of the choice to set.
        """
        try:
            match = [choice.value for choice in self.body].index(value)
            self.set_focus(match)
            
        except ValueError:
            # Couldn't find the given value, ignore?
            pass


class Choices_edit(Popup_edit):
    """
    An edit widget for selecting between a number of pre-chosen choices.
    """
    
    def __init__(self, top, choices, initial = None, title = "Select option", change_callback = None):
        """
        Constructor for Choices_edit widgets.
        """
        super().__init__(top, initial)
        self.choices = choices
        self.title = title
        self.picker = Choices_picker(self)
        self.change_callback = change_callback
    
    def load_popup(self):
        return Widget_dialogue(self.title, self.picker, self.top, submit_callback = self.update)
        
    def update(self):
        """
        Update the value of this widget.
        """
        self.value = self.picker.focus.value
        self.update_label()
        
        # Also update out picked widget to show the currently selected as having default focus.
        self.picker.set_choice(self.value)
        
        if self.change_callback is not None:
            self.change_callback()


class File_edit(Popup_edit):
    """
    An edit widget for picking (existing) file locations.
    """
    
    def __init__(self, top, initial = None, title = "Select file", can_choose_folders = True):
        """
        Constructor for File_edit widgets.
        """
        self.file_selector = silico.interface.urwid.file.browser.File_selector(top, can_choose_folders = can_choose_folders, can_choose_multiple = False)
        super().__init__(top, initial)
        
    def open_popup(self):
        """
        Method called to show the popup the user can interactive with.
        """
        self.top.swap_into_window(self.get_popup(), submit_callback = self.update)
    
    def load_popup(self):
        return self.file_selector
    
    def update(self):
        """
        Update the value of this widget.
        """
        # Get whatever has been selected.
        selected_files = self.file_selector.selected
        self.value = selected_files[-1] if len(selected_files) > 0 else None
        self.file_selector.reset()
        
        self.update_label()
        
        
class Output_edit(File_edit):
    """
    An edit widget for picking a location to save a new file to.
    """
    
    def __init__(self, top, initial = None, folder = False, default_file_name = ""):
        """
        Constructor for File_edit widgets.
        """
        self.file_selector = silico.interface.urwid.file.output.Output_selector(top, default = initial, folder = folder, default_file_name = default_file_name)
        Popup_edit.__init__(self, top, initial)

    def update(self):
        """
        Update the value of this widget.
        """
        # Get whatever has been selected.
        self.value = self.file_selector.value
        
        self.update_label()
        
        
class Method_target_picker(Popup_edit):
    """
    An edit widget for picking part of a method (a calculation, program or destination).
    """
    
    def __init__(self, top, method_targets, initial = None):
        """
        Constructor for File_edit widgets.
        """
        self.method_selector = Method_selector(top, method_targets, one_type_only = True, can_choose_multiple = False)
        super().__init__(top, initial)
        
    def open_popup(self):
        """
        Method called to show the popup the user can interactive with.
        """
        self.top.swap_into_window(self.get_popup(), submit_callback = self.update)
    
    def load_popup(self):
        return self.method_selector
    
    def update(self):
        """
        Update the value of this widget.
        """
        # Get whatever has been selected.
        selected_nodes = self.method_selector.browser.selected_nodes
        if len(selected_nodes) == 0:
            self.value = None
        
        else:
            path = selected_nodes[-1].build_loader_path()[-1]
            method = path[0].resolve_path(path)
            
            self.value = method
            
        self.method_selector.reset()
        
        self.update_label()
        
    def update_label(self):
        """
        Update the label of this widget with the current value.
        """
        self.set_label(str(self.value.name) if self.value is not None else "")
    
    