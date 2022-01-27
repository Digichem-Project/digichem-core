# General imports.
import urwid

# Silico imports.
from silico.exception.base import Silico_exception
from silico.interface.urwid.edit.popup import File_edit, Choices_edit
from silico.interface.urwid.setedit.common import Setedit_widget_parent_mixin
from silico.interface.urwid.layout import Pane


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
    
    def __init__(self, setedit):
        """
        Constructor for Setedit_widget objects.
        """
        self.setedit = setedit
        
        widget_list = self.get_widgets()
        
        urwid.Pile.__init__(self, widget_list)
        
    def discard(self):
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
    
    @value.setter
    def value(self, value):
        """
        Change the edit value of this widget.
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
    
    def __init__(self, setedit):
        """
        Constructor for Setedit_widger objects.
        
        :param edit_widget: A widget to use for editing this settings.
        :param help: Help text to display for this editor.
        """
        self.edit_widget = None
        super().__init__(setedit)
        
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
        self.edit = urwid.Edit((self.title_attr, self.setedit.title + ": "), self.value_to_str(self.setedit.previous_value))
        return super().load_widgets(self.edit)
    
    def discard(self):
        """
        Reset the current value back to the default value.
        """
        self.edit.set_edit_text(self.value_to_str(self.setedit.previous_value))
    
    @property
    def value(self):
        """
        The currently, possibly edited, value of this widget.
        """
        value = self.edit.get_edit_text()
        return self.str_to_value(value)
    
    @value.setter
    def value(self, value):
        """
        Change the edit value of this widget.
        """
        self.edit.set_edit_text(self.value_to_str(value))


class Bool_editor(Single_editor):
    """
    Editer widget for bool values.
    """
    
    def load_widgets(self):
        """
        Load the list of inner widgets we'll use for display.
        """
        # Panic if our previous_value is not True or False.
        if self.setedit.previous_value is not True and self.setedit.previous_value is not False:
            raise Silico_exception("Only True or False are allowed values for checkboxes, not '{}'".format(self.setedit.previous_value))
        
        self.checkbox = urwid.CheckBox((self.title_attr, self.setedit.title + ": "), self.setedit.previous_value)
        
        return super().load_widgets(self.checkbox)
    
    def discard(self):
        """
        Reset the current value back to the default value.
        """
        self.checkbox.set_state(self.setedit.previous_value)
    
    @property
    def value(self):
        """
        The current, possible edited, value of this widget.
        """
        return self.checkbox.get_state()
    
    @value.setter
    def value(self, value):
        """
        Change the edit value of this widget.
        """
        self.checkbox.set_state(value)

    
class Popup_editor(Single_editor):
    """
    ABC for editors that show a popup.
    """
    
    def __init__(self, setedit, popup_widget, *args, **kwargs):
        """
        """
        self.popup_widget = popup_widget
        self._value = setedit.previous_value
        
        super().__init__(setedit, *args, **kwargs)
    
    def load_widgets(self):
        """
        Load the list of inner widgets we'll use for display.
        """
        # We'll wrap this in a columns so we can add our title.
        body = urwid.Columns([
            ('pack', urwid.Text((self.title_attr, "{}: ".format(self.setedit.title)))),
            self.popup_widget
        ], dividechars = 1)
        
        return super().load_widgets(body)
    
    def discard(self):
        """
        Reset the current value back to the default value.
        """
        self.popup_widget.value = self.setedit.previous_value
    
    @property
    def value(self):
        """
        The currently, possibly edited, value of this widget.
        """
        return self.str_to_value(self.popup_widget.value)
    
    @value.setter
    def value(self, value):
        """
        Change the edit value of this widget.
        """
        self.popup_widget.value = self.value_to_str(value)


class File_editor(Popup_editor):
    """
    An editor for picking files.
    """
    
    def __init__(self, setedit):
        popup_widget = File_edit(setedit.top, setedit.previous_value, "File for {}:".format(setedit.title))
        super().__init__(setedit, popup_widget = popup_widget)


class Choices_editor(Popup_editor):
    """
    An editor that allows the user to pick from a number of choices.
    """
     
    def __init__(self, setedit):
        popup_widget = Choices_edit(setedit.top, setedit.choices, setedit.previous_value, "Select option for {}".format(setedit.title))
        super().__init__(setedit, popup_widget = popup_widget)


class Text_list_editor(Min_edit):
    """
    Text editor that appears as part of a list.
    """
    
    # Attribute that controls appearance of the position indicator.
    position_attr = "body"
    
    def __init__(self, editor, position, edit_text):
        """
        Constructor for Text_list_editor widgets.
        
        :param editor: The parent editor object:
        :param position: Starting position (index +1) of this editor.
        :param edit_text: Default value of this editor.
        """
        self.editor = editor
        
        # Replace None values with actual text.
        edit_text = self.editor.value_to_str(edit_text)
                
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
        return self.editor.str_to_value(value)


class List_editor(Setedit_widget):
    """
    Editor widget for list values.
    """
    
    def __init__(self, setedit, strip_empty = True):
        """
        Constructor for Bool_editor widgets.
        
        :param title: The title of this editor.
        :param values: List of starting values.
        :param help: Help text to display for this editor.
        :param strip_empty: Whether to remove empty values in the list.
        """
        self.strip_empty = strip_empty
        
        super().__init__(setedit)
        
    def discard(self):
        """
        Reset the current value back to the default value.
        """
        self.inner_pile.contents = [(field, self.inner_pile.options()) for field in self.get_fields(self.setedit.previous_value)]
        
    def load_widgets(self):
        """
        Load the list of inner widgets we'll use for display.
        """
        # An inner pile we'll use for adding/removing values.
        self.inner_pile = urwid.Pile(self.get_fields(self.setedit.previous_value))
        
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
        if values is not None:
            fields = [self.get_field(index +1, value) for index, value in enumerate(values)]
            
        else:
            fields = []
        
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
        field =  Text_list_editor(self, position, edit_text = value)
        
        # Add a signal to watch for changes.
        urwid.connect_signal(field, 'postchange', lambda *data: self.adjust())
        
        return urwid.AttrMap(field, self.edit_attr)
    
    @property
    def value(self):
        """
        The currently, possible edited, value of this widget.
        """
        return [field.base_widget.value for field, options in self.inner_pile.contents if field.base_widget.value is not None or not self.strip_empty]
    
    @value.setter
    def value(self, value):
        """
        Change the edit value of this widget.
        """
        self.inner_pile.contents = [(field, self.inner_pile.options()) for field in self.get_fields(value)]
        
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
    
    def __init__(self, editor, key, value, change_callback = None):
        """
        Constructor for Dict_item_editor objects.
        
        :param editor: The parent Dict_editor we belong to.
        :param key: The key of this dictionary item.
        :param value: The value of this dictionary item.
        :param change_callback: Function to call when this editor changes value.
        """
        self.editor = editor
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
        return (self.editor.str_to_value(self.key_edit.get_edit_text()), self.editor.str_to_value(self.value_edit.get_edit_text()))        


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
        return Dict_item_editor(self, key, value, self.adjust)
    
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
    
    @value.setter
    def value(self, value):
        super(self.__class__, self.__class__).value.fset(self, value)

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

# TOOD: Rename?
class Sub_setedit(Setedit_widget, Setedit_widget_parent_mixin):
    """
    A Setedit object for editing sub options.
    """
    
    padding_amount = 4
    
    @property
    def child_widgets(self):
        """
        A shortcut to the property where the child widgets are stored.
        
        This must be defined by the inheriting class because it depends on the type of the inheriting class, for example self.body for listboxes.
        """
        return [widget for widget, options in self.inner_pile.contents]
            
    def load_widgets(self):
        """
        Load the list of inner widgets we'll use for display.
        """
        # Prepend a divider.
        child_widgets = [urwid.Divider("-")]
        child_widgets.extend(self.load_child_widgets(self.setedit.get_children()))
        
        self.inner_pile = urwid.Pile(child_widgets)
        
        return [
            # Widget for our title.
            urwid.Text((self.title_attr, self.setedit.title)),
            # Help text.
            self.get_help_widget(self.setedit.help),
            # Actual contents.
            urwid.Padding(self.inner_pile, left = self.padding_amount)
        ]
        
    def discard(self):
        """
        Reset the current value back to the default value.
        """
        return Setedit_widget_parent_mixin.discard(self)
    
    @property
    def value(self):
        """
        The currently, possible edited, value of this widget.
        """
        val =  {}
        for child_widget in self.child_setedit_widgets.values():
            value = child_widget.value
            val[child_widget.setedit.title] = value
                
        return val
    
    @value.setter
    def value(self, value):
        """
        Change the edit value of this widget.
        
        :param value: The new value to set; a dictionary of values mapping to the child widgets of this widget.
        """
        # We'll cache this dict temporarily because otherwise we'll be recreating it each time.
        child_widgets = self.child_setedit_widgets
        for key, sub_value in value:
            child_widgets[key].value = sub_value

      
class Solo_sub_editor(Sub_setedit):
    """
    A widget for displaying a set of sub-options as the main editing item.
    """
            
    def load_widgets(self):
        """
        Load the list of inner widgets we'll use for display.
        """
        # Prepend a divider.
        child_widgets = list(self.load_child_widgets(self.setedit.get_children()))
        
        self.inner_pile = urwid.Pile(child_widgets)
        
#         return [
#             ## Widget for our title.
#             #urwid.Text((self.title_attr, self.setedit.title)),
#             ## Help text.
#             #self.get_help_widget(self.setedit.help),
#             ## A divider.
#             #urwid.Divider('-'),
#             # Actual contents.
#             self.inner_pile
#         ]
        return [
            Pane(urwid.Pile([
                self.get_help_widget(self.setedit.help),
                urwid.Divider('-'),
                self.inner_pile
                ]),
            self.setedit.title)
        ] 

