# General imports.
import urwid
import pathlib

# Silico imports.
import silico.interface.urwid.file.browser
from silico.interface.urwid.dialogue import Widget_dialogue


class Min_edit(urwid.Edit):
    """
    A fixed urwid edit widget that always leaves enough space for the cursor.
    """
    
    def pack(self, *args, **kwargs):
        """
        Modified pack method to fix an urwid bug
        """
        packed = super().pack(*args, **kwargs)
        # Always make sure we have enough space for our cursor.
        return (packed[0] +1, packed[1])


class Setedit_widget(urwid.Pile):
    """
    Top class for widgets that handle display of Setedit objects.
    """
    
    title_attr = "bold"
    edit_attr = "editable"
    help_attr = "body"
    
    def __init__(self, setedit, has_divider = True):
        """
        Constructor for Setedit_widget objects.
        
        :param widget_list: List of widgets to stack vertically.
        :param has_divider: Whether to show a divider after this widget.
        """
        self.setedit = setedit
        self.divider = urwid.Divider()
        
        widget_list = self.get_widgets()
        
        if has_divider:
            widget_list.append(self.divider)
            
        urwid.Pile.__init__(self, widget_list)
        
    def reset(self):
        """
        Reset the current value back to the default value.
        """
        raise NotImplementedError("Implement in subclass")
        
    def load_widgets(self):
        """
        Load the list of inner widgets we'll use for display.
        """
        raise NotImplementedError("Implement in subclass")
    
    def get_widgets(self, reload = False):
        """
        Get the list of inner widgets we'll use for display.
        """
        # No need to cache?
        return self.load_widgets()
    
    @property
    def value(self):
        """
        Get the possibly edited value of this widget.
        """
        raise NotImplementedError("Implement in subclass")
    
    def get_help_widget(self, help_text):
        """
        Return a widget used to display a help message.
        """
        return urwid.Text((self.help_attr, help_text if help_text is not None else ""))
    
    @classmethod
    def value_to_str(self, value):
        """
        Convert a value to a string.
        
        This function handles mapping of 'None' values.
        """
        return str(value) if value is not None else ""
    
    @classmethod
    def str_to_value(self, value):
        """
        Convert a string to a real value.
        
        This function handles mapping of 'None' values.
        """
        return value if value != "" else None
    
    @classmethod
    def class_from_type(self, vtype):
        """
        Get a concrete Setedit_widget suitable for a certain type.
        
        :param vtype: The name of the type (as a string). If in doubt try type(value).__name__
        :returns: A suitable class which will have a constructor of the form __init__(title, value(s), help).
        """
        if vtype == "bool":
            return Bool_editor
        
        elif vtype == "Path":
            return File_editor
        
        elif vtype == "choices":
            return Choices_editor
        
        elif vtype == "list" or vtype == "tuple":
            return List_editor
        
        elif vtype == "dict":
            return Dict_editor
        
        elif vtype == "Options":
            return Sub_setedit
        
        else:
            return Text_editor


class Single_editor(Setedit_widget):
    """
    A Setedit_widget for changing a single value.
    """
    
    def __init__(self, setedit, has_divider = True):
        """
        Constructor for Setedit_widger objects.
        
        :param edit_widget: A widget to use for editing this settings.
        :param help: Help text to display for this editor.
        :param has_divider: Whether to show a divider after this widget.
        """
        self.edit_widget = None
        super().__init__(setedit, has_divider)
        
    def load_widgets(self, edit_widget):
        """
        Load the list of inner widgets we'll use for display.
        """
        return [
            urwid.AttrMap(edit_widget, self.edit_attr),
            self.get_help_widget(self.setedit.help)
        ]


class Text_editor(Single_editor):
    """
    Editor widget for generic text types.
    """
        
    def load_widgets(self):
        """
        Load the list of inner widgets we'll use for display.
        """
        self.edit = urwid.Edit((self.title_attr, self.setedit.title + ": "), self.value_to_str(self.setedit.starting_value))
        return super().load_widgets(self.edit)
    
    def reset(self):
        """
        Reset the current value back to the default value.
        """
        self.edit.set_edit_text(self.value_to_str(self.setedit.starting_value))
    
    @property
    def value(self):
        """
        The currently, possibly edited, value of this widget.
        """
        value = self.edit.get_edit_text()
        return self.str_to_value(value)


class Bool_editor(Single_editor):
    """
    Editer widget for bool values.
    """
    
    def load_widgets(self):
        """
        Load the list of inner widgets we'll use for display.
        """
        self.checkbox = urwid.CheckBox((self.title_attr, self.setedit.title + ": "), self.setedit.starting_value)
        
        return super().load_widgets(self.checkbox)
    
    def reset(self):
        """
        Reset the current value back to the default value.
        """
        self.checkbox.set_state(self.setedit.starting_value)
    
    @property
    def value(self):
        """
        The current, possible edited, value of this widget.
        """
        return self.checkbox.get_state()
    
class Popup_editor(Single_editor):
    """
    ABC for editors that show a popup.
    """
    
    def __init__(self, setedit, *args, **kwargs):
        """
        """
        self._popup = None
        self._value = setedit.starting_value
        
        # The button that can be used to launch our popup.
        self.button = urwid.Button("", lambda button: self.open_popup())
        self.update_label()
        
        super().__init__(setedit, *args, **kwargs)
    
    def load_widgets(self):
        """
        Load the list of inner widgets we'll use for display.
        """
        # We'll wrap this in a columns so we can add our title.
        body = urwid.Columns([
            ('pack', urwid.Text((self.title_attr, "{}: ".format(self.setedit.title)))),
            self.button
        ], dividechars = 1)
        
        return super().load_widgets(body)
    
    def open_popup(self):
        """
        Method called to show the popup the user can interactive with.
        """
        self.setedit.top.popup(self.get_popup())
        
    def close_popup(self):
        """
        Method called to hide the popup.
        """
        self.setedit.top.close_popup(self.get_popup())
        
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
    
    def reset(self):
        """
        Reset the current value back to the default value.
        """
        self._value = self.setedit.starting_value
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
        self.button.set_label(self.value_to_str(self._value))
    
    @property
    def value(self):
        """
        The currently, possibly edited, value of this widget.
        """
        return self.str_to_value(self._value)


class File_editor(Popup_editor):
    """
    An editor for picking files.
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.file_selector = silico.interface.urwid.file.browser.File_selector(self.setedit.top, title = "File for {}:".format(self.setedit.title), can_choose_folders = True, can_choose_multiple = False)
    
    def open_popup(self):
        """
        Method called to show the popup the user can interactive with.
        """
        self.setedit.top.swap_into_window(self.get_popup(), submit_callback = self.update)
    
    def load_popup(self):
        return self.file_selector

    def update(self):
        """
        Update the value of this widget.
        """
        selected_files = self.file_selector.selected
        self._value = selected_files[-1] if len(selected_files) > 0 else None
        self.file_selector.reset()
        
        self.update_label()
        
    @classmethod
    def value_to_str(self, value):
        """
        Convert a value to a string.
        
        This function handles mapping of 'None' values.
        """
        return str(value) if value is not None else ""
    
    @classmethod
    def str_to_value(self, value):
        """
        Convert a string to a real value.
        
        This function handles mapping of 'None' values.
        """
        return pathlib.Path(value) if value != "" and value is not None else None


class Choices_widget(urwid.AttrMap):
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
        self.picker.editor_widget.update()
        self.picker.editor_widget.close_popup()
        
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
        self.button.set_label(value)


class Choices_picker(urwid.ListBox):
    """
    Widget that allows picking from a number of choices.
    """
    
    def __init__(self, editor_widget):
        """
        """
        super().__init__(urwid.SimpleFocusListWalker([]))
        
        # Keep our widget for later.
        self.editor_widget = editor_widget
        
        for choice in editor_widget.setedit.choices:
            self.body.append(self.get_widget(choice))
            
        self.set_choice(self.editor_widget._value)
        
    def get_widget(self, value):
        """
        Get one of the choices widget (a button) we will display.
        """
        return Choices_widget(self.editor_widget.value_to_str(value), self)
    
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


class Choices_editor(Popup_editor):
    """
    An editor that allows the user to pick from a number of choices.
    """
    
    def __init__(self, setedit, *args, **kwargs):
        super().__init__(setedit, *args, **kwargs)
        self.picker = Choices_picker(self)
    
    def load_popup(self):
        return Widget_dialogue("Select option for {}".format(self.setedit.title), self.picker, self.setedit.top, submit_callback = self.update)
        
    def update(self):
        """
        Update the value of this widget.
        """
        self._value = self.picker.focus.value
        self.update_label()
        
        # Also update out picked widget to show the currently selected as having default focus.
        self.picker.set_choice(self._value)


class Text_list_editor(Min_edit):
    """
    Text editor that appears as part of a list.
    """
    
    # Attribute that controls appearance of the position indicator.
    position_attr = "body"
    
    def __init__(self, position, edit_text):
        """
        Constructor for Text_list_editor widgets.
        
        :param position: Starting position (index +1) of this editor.
        :param edit_text: Default value of this editor.
        """
        
        # Replace None values with actual text.
        edit_text = Setedit_widget.value_to_str(edit_text)
                
        super().__init__("", edit_text)
        
        # Set our position.
        self.set_position(position)
        
    def set_position(self, position):
        """
        Tell this widget where it appears in the list.
        """
        self.set_caption((self.position_attr, "{}) ".format(position)))
        
    @property
    def value(self):
        """
        The currently, possible edited, value of this widget.
        """
        value = self.get_edit_text()
        return Setedit_widget.str_to_value(value)


class List_editor(Setedit_widget):
    """
    Editor widget for list values.
    """
    
    def __init__(self, setedit, has_divider = True, strip_empty = True):
        """
        Constructor for Bool_editor widgets.
        
        :param title: The title of this editor.
        :param values: List of starting values.
        :param help: Help text to display for this editor.
        :param has_divider: Whether to show a divider after this widget.
        :param strip_empty: Whether to remove empty values in the list.
        """
        self.strip_empty = strip_empty
        
        super().__init__(setedit, has_divider)
        
    def reset(self):
        """
        Reset the current value back to the default value.
        """
        self.inner_pile.contents = [(field, self.inner_pile.options()) for field in self.get_fields(self.setedit.starting_value)]
        
    def load_widgets(self):
        """
        Load the list of inner widgets we'll use for display.
        """
        # An inner pile we'll use for adding/removing values.
        self.inner_pile = urwid.Pile(self.get_fields(self.setedit.starting_value))
        
        return [
            # Widget for our title.
            urwid.Text((self.title_attr, self.setedit.title)),
            # Help text.
            self.get_help_widget(self.setedit.help),
            # Actual contents.
            self.inner_pile
        ]
                
    def get_fields(self, values):
        """
        Build a list of edit fields based on a value.
        """
        fields = [self.get_field(index +1, value) for index, value in enumerate(values)]
        
        # Add one final editor so we can add more values.
        fields.append(self.get_field(len(fields) +1, None))
        
        return fields
        
    def get_field(self, position, value):
        """
        Get a new editable field widget.
        
        :param position: The position of the new field in the list.
        :param value: The starting position of value.
        :returns: The new widget
        """
        field =  Text_list_editor(position, edit_text = value)
        
        # Add a signal to watch for changes.
        urwid.connect_signal(field, 'postchange', lambda *data: self.adjust())
        
        return urwid.AttrMap(field, self.edit_attr)
    
    @property
    def value(self):
        """
        The currently, possible edited, value of this widget.
        """
        return [field.base_widget.value for field, options in self.inner_pile.contents if field.value is not None or not self.strip_empty]
        
    def adjust(self):
        """
        Adjust the length of this list editor.
        """
        # Check to see if the last field has something in it or not.
        if self.inner_pile.contents[-1][0].base_widget.value is not None:
            # Add new final editor.
            self.inner_pile.contents.append((self.get_field(len(self.contents) +1, None), self.options()))

        elif len(self.inner_pile.contents) > 1 and self.inner_pile.contents[-2][0].base_widget.value is None:
            # We have too many empty editors, remove the last.
            self.inner_pile.contents.pop()
                        
            # Check again.
            self.adjust()
            
        # Update positions.
        for index, (field, options) in enumerate(self.inner_pile.contents):
            field.base_widget.set_position(index +1)


class Dict_item_editor(urwid.Columns):
    """
    Editor that appears as part of dictionary.
    """
    
    edit_attr = "editable"
    body_attr = "body"
    
    def __init__(self, key, value, change_callback = None):
        """
        Constructor for Dict_item_editor objects.
        
        :param key: The key of this dictionary item.
        :param value: The value of this dictionary item.
        :param change_callback: Function to call when this editor changes value.
        """
        # We have two separate editable entities.
        self.key_edit = Min_edit(edit_text = Setedit_widget.value_to_str(key))
        self.value_edit = Min_edit(edit_text = Setedit_widget.value_to_str(value))
                
        # Add a signal to watch for changes.
        urwid.connect_signal(self.key_edit, 'postchange', lambda *data: change_callback())
        urwid.connect_signal(self.value_edit, 'postchange', lambda *data: change_callback())
                
        # Call parent.
        super().__init__([
            ('pack', urwid.AttrMap(self.key_edit, self.edit_attr)),
            ('pack', urwid.Text(": ")),
            ('pack', urwid.AttrMap(self.value_edit, self.edit_attr))],
            min_width = 5
        )
    
    @property
    def value(self):
        """
        Get the possibly edited value of this field.
        """
        return (Setedit_widget.str_to_value(self.key_edit.get_edit_text()), Setedit_widget.str_to_value(self.value_edit.get_edit_text()))
        


class Dict_editor(List_editor):
    """
    Editor widget for dicts
    """
        
    def get_field(self, key, value):
        """
        Get a new editable field widget.
        
        :param key: The key of this dictionary item.
        :param value: The value of this dictionary item.
        """
        return Dict_item_editor(key, value, self.adjust)
    
    def get_fields(self, values):
        """
        Build a list of edit fields based on a value.
        """
        fields = [self.get_field(key, value) for key, value in values.items()]
        
        # Add one final editor so we can add more values.
        fields.append(self.get_field(None, None))
        
        return fields
                    
    @property
    def value(self):
        """
        The currently, possible edited, value of this widget.
        """
        val =  {}
        for field, options in self.inner_pile.contents:
            key, value = field.value
            if not (key is None and value is None) or not self.strip_empty:
                val[key] = value
                
        return val

    def adjust(self):
        """
        Adjust the length of this list editor.
        """
        # Check to see if the last field has something in it or not.
        if self.inner_pile.contents[-1][0].value != (None, None):
            # Add new final editor.
            self.inner_pile.contents.append((self.get_field(None, None), self.options()))

        elif len(self.inner_pile.contents) > 1 and self.inner_pile.contents[-2][0].value == (None, None):
            # We have too many empty editors, remove the last.
            self.inner_pile.contents.pop()
                        
            # Check again.
            self.adjust()


class Sub_setedit(Setedit_widget):
    """
    A Setedit object for editing sub options.
    """
    
    
    def reset(self):
        """
        Reset the current value back to the default value.
        """
        self.inner_pile.contents = [(field, self.inner_pile.options()) for field in self.get_fields()]
        
    def load_widgets(self):
        """
        Load the list of inner widgets we'll use for display.
        """
        self.inner_pile = urwid.Pile(self.get_fields())
        
        return [
            # Widget for our title.
            urwid.Text((self.title_attr, self.setedit.title)),
            # Help text.
            self.get_help_widget(self.setedit.help),
            # Actual contents.
            urwid.Padding(self.inner_pile, left = 2)
        ]
        
    def get_fields(self):
        """
        Get the child widgets to display.
        """
        fields = []
        
        # Ask our setedit for its children.
        children = self.setedit.get_children()
        for child_setedit in children:
            fields.append(self.class_from_type(child_setedit.vtype)(child_setedit, has_divider = len(fields) != len(children)-1))
            
        return fields
    
    @property
    def value(self):
        """
        The currently, possible edited, value of this widget.
        """
        val =  {}
        for field, options in self.inner_pile.contents:
            key, value = field.value
            val[key] = value
                
        return val

