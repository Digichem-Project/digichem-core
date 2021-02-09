import urwid
from silico.interface.urwid.misc import Tab_pile
from silico.config.configurable.options import Options

class Option_editor(urwid.Pile):
    """
    Top-level class for those that edit options in a Configurable.
    """
    
    def __init__(self, original_widget, configurable, option, dictobj = None):
        """
        """
        help_msg = option.help
        choices = option.choices(configurable)
        # If we have choices, append them.
        if choices is not None:
            if help_msg[-1] != ".":
                help_msg += "."
            help_msg += " Possible values are:\n" + "\n".join(["- '{}'".format(choice) for choice in choices])
        
        
        # Construct our Pile of editable widget and help.
        super().__init__([
            urwid.AttrMap(original_widget, 'editable'),
            urwid.Text(help_msg)
        ])
        self.configurable = configurable
        self.option = option
        self.dictobj = dictobj
        
    def save(self):
        """
        Save the current value of this editor.
        """
        self.option.__set__(self.configurable, self.value, dictobj = self.dictobj)
    
    @property
    def value(self):
        """
        The currently, possible edited, value of this widget.
        
        All subclasses should define an implementation of this property.
        """
        raise NotImplementedError()
    
    @classmethod
    def from_type(self, configurable, option, padding = 4, **kwargs):
        """
        Get an appropriate urwid widget depending on the type of an option.
        """
        if isinstance(option, Options):
            return Options_editor(configurable, option, padding)
        elif option.type is bool:
            return Bool_editor(configurable, option, **kwargs)
        elif option.type is list or option.type is tuple:
            return List_editor(configurable, option, **kwargs)
        elif option.type is dict:
            return Dict_editor(configurable, option, **kwargs)
        else:
            return Text_editor(configurable, option, **kwargs)


class Text_editor(Option_editor):
    """
    Editer widget for generic text types.
    """
    
    def __init__(self, configurable, option, **kwargs):
        """
        Constructor for Bool_editor widgets.
        """
        value = option.getraw(configurable, **kwargs)
        self.edit = urwid.Edit(('boldnode', option.name +": "), edit_text = str(value) if value is not None else "")
        super().__init__(self.edit, configurable, option, **kwargs)
    
    @property
    def value(self):
        """
        The currently, possibly edited, value of this widget.
        """
        value = self.edit.get_edit_text()
        return value if value != "" else None


class Text_list_editor(urwid.Edit):
    """
    Text editor that appears as part of a list.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Constructor for urwid.Edit widgets.
        """
        if 'edit_text' in kwargs and kwargs['edit_text'] is None:
            kwargs['edit_text'] = ""
        urwid.Edit.__init__(self, *args, **kwargs)
    
    def pack(self, *args, **kwargs):
        """
        Modified pack method to fix an urwid bug
        """
        packed = super().pack(*args, **kwargs)
        # Always make sure we have enough space for our cursor.
        return (packed[0] +1, packed[1])
    
    @property
    def value(self):
        """
        The currently, possible edited, value of this widget.
        """
        value = self.get_edit_text()
        return value if value != "" else None
        
            
class Bool_editor(Option_editor):
    """
    Editer widget for bool values.
    """
    
    def __init__(self, configurable, option, **kwargs):
        """
        Constructor for Bool_editor widgets.
        """
        self.checkbox = urwid.CheckBox(('boldnode', option.name), option.getraw(configurable, **kwargs))
        super().__init__(self.checkbox, configurable, option, **kwargs)
    
    @property
    def value(self):
        """
        The currently, possible edited, value of this widget.
        """
        return self.checkbox.get_state()
    
class List_editor(Option_editor):
    """
    Editor widget for list values.
    """
    
    def __init__(self, configurable, option, *, strip_empty = True, **kwargs):
        """
        Constructor for Bool_editor widgets.
        
        :param strip_empty: Whether to remove empty values in the list.
        """
        self.strip_empty = strip_empty
        
        # Build our list of editors.
        editors = [urwid.Text(('boldnode', option.name + ": "))]
        self.fields = []
        self.pile = urwid.Pile(editors)
        self.build_fields(configurable, option, **kwargs)
        
        # Set focus properly.
        self.pile.focus_position = len(self.pile.contents) -1
        
        super().__init__(self.pile, configurable, option, **kwargs)
        
    def build_fields(self, configurable, option, **kwargs):
        """
        Build a list of edit fields based on a value.
        """
        # Add edit fields.
        for num, value in enumerate(option.getraw(configurable, **kwargs)):
            # Create a new text edit.
            self.new_field(num, value)
            
        # Add one final editor so we can add more values.
        self.new_field(len(self.fields), None)
        
    def new_field(self, num, value):
        """
        Get a new editable field widget.
        
        :returns: The new widget, which is also appended to the fields attribute.
        """
        field =  Text_list_editor("{}) ".format(num+1), edit_text = value)
        
        # Add a signal to watch for changes.
        urwid.connect_signal(field, 'postchange', self.adjust)
        
        self.fields.append(field)
        self.pile.contents.append((field, self.pile.options()))
        return field
    
    @property
    def value(self):
        """
        The currently, possible edited, value of this widget.
        """
        return [field.value for field in self.fields if field.value is not None or not self.strip_empty]
        
    def adjust(self, *args):
        """
        Adjust the length of this list editor.
        """
        # Check to see if the last field has something in it or not.
        if self.fields[-1].value is not None:
            # Add new final editor.
            #self.append_editor()
            self.new_field(len(self.fields), None)
        elif len(self.fields) > 1 and self.fields[-2].value is None:
            # We have too many empty editors, remove the last.
            self.pile.contents.pop()
            self.fields.pop()
            
            # Check again.
            self.adjust()


class Dict_editor(List_editor):
    """
    Editor widget for dicts
    """
    
    def build_fields(self, configurable, option, **kwargs):
        """
        Build a list of edit fields based on a value.
        """
        # Add edit fields.
        value = option.getraw(configurable, **kwargs)
        for name, value in (value.items() if value is not None else {}.items()):
            # Create a new text edit.
            self.new_field(name, value)
        
#         value = option.__get__(configurable, type(configurable))
#         for name, value in (value.items() if value is not None else {}.items()):
#             # Create a new text edit.
#             self.new_field(name, value)
            
        # Add one final editor so we can add more values.
        self.new_field(None, None)
        
    def new_field(self, name, value):
        """
        Get a new editable field widget.
        
        :returns: The new widget, which is also appended to the fields attribute.
        """
        name_part = Text_list_editor(edit_text = name)
        value_part = Text_list_editor(edit_text = value)
        
        field = (name_part, value_part)
        
        # Add a signal to watch for changes.
        urwid.connect_signal(name_part, 'postchange', self.adjust)
        urwid.connect_signal(value_part, 'postchange', self.adjust)
        
        # Add to list of fields.
        self.fields.append(field)
        
        # Wrap in columns and return.
        widget = urwid.Columns([
            ('pack', name_part),
            ('pack', urwid.Text(": ")),
            ('pack', value_part)
        ],
            min_width = 5
        )
        
        self.pile.contents.append((widget, self.pile.options()))
        return widget
        
    @property
    def value(self):
        """
        The currently, possible edited, value of this widget.
        """
        val =  {
            name.value: value.value if value.value is not None else ""
            for name,value
            in self.fields
            if not (name.value is None and value.value is None) or not self.strip_empty
        }
        return val
    
    def adjust(self, *args):
        """
        Adjust the length of this list editor.
        """
        # Check to see if the last field has something in it or not.
        if self.fields[-1][0].value is not None or self.fields[-1][1].value is not None:
            # Add new final editor.
            #self.append_editor()
            self.new_field(None, None)
        elif len(self.fields) > 1 and self.fields[-2][0].value is None and self.fields[-2][1].value is None:
            # We have too many empty editors, remove the last.
            self.pile.contents.pop()
            self.fields.pop()
            
            # Check again.
            self.adjust()


class Options_editor(Tab_pile):
    """
    Editor widget for groups of options.
    """
        
    def __init__(self, configurable, option, padding = 4):
        """
        Constructor for Configurable_editor objects.
        """    
        # Save our configurable.
        self.configurable = configurable
        
        # Build our editing widgets.
        self.fields = []
        widgets = self.get_widgets(configurable, option, padding)
            
        # Remove the last widget (unneeded Divider).
        widgets.pop()
        
        # Wrap our walker in a box.
        edit_list = urwid.Pile(widgets)
                
        super().__init__([
            urwid.Text(('boldnode', option.name)),
            urwid.Text(option.help),
            urwid.Padding(edit_list, left = padding)
        ])
        
        # Save our Options (plural) object.
        if isinstance(option, Options):
            self.options_obj = option
        else:
            self.options_obj = None
        self.configurable = configurable
        
    def get_widgets(self, configurable, option_or_configurable, padding = 4):
        """
        Get a list of widgets that will make up the body of this editor.
        
        Editable widgets will also be added to the 'fields' list attribute.
        
        :return: List of widgets (some editable, some not).
        """
        # Build our editing widgets.
        widgets = []
        for name,sub_option in option_or_configurable.OPTIONS.items():
            if not sub_option.no_edit:
                # Get our edit field.
                field = Option_editor.from_type(configurable, sub_option, padding = padding + 4, dictobj = option_or_configurable.getraw(configurable) if isinstance(option_or_configurable, Options) else None)
                                
                # Add.
                self.fields.append(field)
                widgets.append(field)
                widgets.append(urwid.Divider('─'))
                
        return widgets    
        
    def save(self, *args, **kwargs):
        """
        Save any changes.
        """
        for field in self.fields:
            field.save()
            
        # Also configure our Options object.
        if getattr(self, 'options_obj', None) is not None:
            self.options_obj.configure(self.configurable)
