from silico.exception.configurable import Configurable_option_exception,\
    Missing_option_exception, Disallowed_choice_exception
import itertools

class Option():
    """
    Class for specifying an option in a configurable.
    
    Options are descriptors that perform type checking and other functionality for Configurables; they expose the options that a certain configurable expects.
    """
    
    def __init__(self, name = None, *, default = None, help = None, choices = None, validate = None, list_type = None, type = None, rawtype = None, exclude = None, required = False, no_edit = False):
        """
        Constructor for Configurable Option objects.
        
        :param name: The name of this option. If None is given this will be determined automatically from the name of the attribute this option is stored under.
        :param default: Default value for this option. Alternatively, default can be a callable which will be called with 2 arguments: this Option object and the owning Configurable object and should return the default value. Be careful about setting mutable options (eg, lists, dicts) as default, if one object modifies the object this will be propagated to all objects that have this option.
        :param help: Descriptive help string
        :param choices: An optional iterable of valid choices for this option.
        :param validate: Function called to check that the given value is valid. The function will be called with 3 arguments: this Option object, the owning Configurable object and the value being set, and should return True or False as appropriate.
        :param list_type: If given, indicates that this option should accept a number of elements which will be stored in an object of type list_type (commonly list or tuple, but anything with a similar interface should work). If list_type and type are both given, then type will be used to set the type of each element in the list.
        :param type: A callable that is used to set the type of value.
        :param exclude: A list of strings of the names of attributes that this option is mutually exclusive with.
        :param required: Whether this option is required or not.
        :param no_edit: Flag to indicate that this option shouldn't be edited.
        """
        self.name = name
        # TODO: Need to be smarter about storing mutable objects (most notable, list and dict, normally empty ones) as default
        # When a default is accessed, it is this actual object that is returned, hence modifications would be shared by all configurable options of the same type.
        # Perhaps a solution is to on demand copy() the default and store the resulting clone?
        self._default = default
        self.list_type = list_type
        self.type = type
        # TODO: This is probably deprecated
        self.rawtype = type if rawtype is None else rawtype
        self.help = help
        self.choices = choices if choices is not None else []
        self._validate = validate if validate is not None else self.default_validate
        self.exclude = exclude if exclude is not None else []
        self.required = required
        self.no_edit = no_edit
        
        # If we are a list_type and a default of None has been given, change the default to an empty list.
        if self.list_type is not None and self._default is None:
            self._default = []
        
        # If we are a sub-object (ie, part of a dict), this is a hierarchy of names of the Options object we are owned by.
        self.parents = []
        
        # By definition, Options that are required can have no default, so we'll delete this attribute.
        if self.required:
            del(self._default)
    
    def default_validate(self, option, configurable, value):
        """
        A function used as the default for _validate; always returns True
        """
        return True
    
    @property
    def full_name(self):
        """
        The full name/path of this option, including any parents.
        """
        return ": ".join(itertools.chain([parent.name for parent in self.parents], (self.name,)))
            
    def add_parent(self, parent):
        """
        Add an owning parent Options object to this Option object.
        
        This method is called by the parent Options object when this Option is added to it.
        """
        self.parents.insert(0, parent)

    def __set_name__(self, cls, name):
        """
        Called automatically during class creation, allows us to know the attribute name we are stored under.
        """
        self.name = name if self.name is None else self.name


    def __get__(self, owning_obj, cls = None):
        """
        Compute/retrieve the value of this option.
        """
        if owning_obj is None:
            return self
        
        return self.get_from_dict(owning_obj, owning_obj._configurable_options)


    def get_from_dict(self, owning_obj, dict_obj):
        """
        Compute/retrieve the value of this option which is stored in a given dictionary.
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param dict_obj: The dict in which the value of this Option is stored. In most cases, the value of this option is evaluated simply as dict_obj[self.name]
        """
        try:            
            return dict_obj[self.name]
        except KeyError:
            # No value set, return our default value (if we have one).
            try:
                return self.default(owning_obj)
            except AttributeError:
                # No value set and no default, panic.
                raise Missing_option_exception(owning_obj, self) from None


    def __set__(self, owning_obj, value):
        """
        Set the value of this option.
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param value: The new value to set.
        """
        self.set_into_dict(owning_obj, owning_obj._configurable_options, value)


    def set_into_dict(self, owning_obj, dict_obj, value):
        """
        Set the value of this option into a specified dict object.
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param dict_obj: The dict in which the value of this Option is stored. In most cases, the value of this option is evaluated simply as dict_obj[self.name]
        :param value: The new value to set.
        """
        dict_obj[self.name] = value


    def __delete__(self, owning_obj):
        """
        Delete the explicit value of this option, resorting to the default (if one is given).
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        """
        self.set_default(owning_obj, owning_obj._configurable_options)


    def set_default(self, owning_obj, dict_obj):
        """
        Reset this option to default.
        
        Note that this method is an alias for calling del() on this attribute.
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param dict_obj: The dict in which the value of this Option is stored. In most cases, the value of this option is evaluated simply as dict_obj[self.name]
        """
        dict_obj.pop(self.name, None)


    def default(self, owning_obj):
        """
        Get the default value of this Option.
        
        :raises AttributeError: If this Option object is required.
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        """
        if not callable(self._default):
            return self._default
        else:
            return self._default(self, owning_obj)


    def is_default(self, dict_obj):
        """
        Whether the value of this option is currently the default or not.
        
        :param dict_obj: The dict in which the value of this Option is stored. In most cases, the value of this option is evaluated simply as dict_obj[self.name]
        """
        return not self.name in dict_obj


    def validate(self, owning_obj, dict_obj = None):
        """
        Validate the value of this option.
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param dict_obj: The dict in which the value of this Option is stored. In most cases, the value of this option is evaluated simply as dict_obj[self.name]
        """
        if dict_obj is None:
            dict_obj = owning_obj._configurable_options
        
        value = self.get_from_dict(owning_obj, dict_obj)
        
        # Try and set the type.
        if not self.is_default(dict_obj) and value is not None:
            # If we are a list type, convert the type of value and also each element.
            if self.list_type is not None:
                # First convert the list itself.
                try:
                    value = self.list_type(value)
                except (TypeError, ValueError) as e:
                    raise Configurable_option_exception(owning_obj, self, "value '{}' of type '{}' could not be converted to the list-like type '{}'".format(value, type(value).__name__, self.list_type)) from e
                
                # Now check each element.
                if self.type is not None:
                    for index, element in enumerate(value):
                        try:
                            value[index] = self.type(element)
                            
                        except (TypeError, ValueError) as e:
                            raise Configurable_option_exception(owning_obj, self, "item '{}') '{}' of type '{}' is of invalid type".format(index, element, type(element).__name__)) from e
            else:
                # Not a list type, just convert.
                try:
                    value = self.type(value) if self.type is not None else value
                    self.set_into_dict(owning_obj, dict_obj, value)
                except (TypeError, ValueError) as e:
                    raise Configurable_option_exception(owning_obj, self, "value '{}' of type '{}' is of invalid type".format(value, type(value).__name__)) from e
        
        # If a list of possible choices has been set, check those now.
        if len(self.choices) != 0:
            # If we are a list type, we'll check each item in value (rather than value itself).
            if self.list_type is not None:
                # Check each of our values, storing each in a new list in case they get changed.
                values = value
                new_values = []
                
                for sub_value in values:
                    new_values.append(self.validate_choices(sub_value, owning_obj, dict_obj))
                    
                value = new_values
                    
            else:
                # Not a list type, only a single value.
                value = self.validate_choices(value, owning_obj, dict_obj)
                
            # Now we need to set our value again incase it changed from validation.
            self.set_into_dict(owning_obj, dict_obj, value)
                
        # Check the value is valid.
        if not self._validate(self, owning_obj, value):
            # Invalid.
            raise Configurable_option_exception(owning_obj, self, "value '{}' of type '{}' is invalid".format(value, type(value).__name__))
        
        # Finally, if our value is equivalent to our default, we'll actually delete our value and use the default instead.
        if value == self.default(owning_obj):
            self.set_default(owning_obj, dict_obj)


    def validate_choices(self, value, owning_obj, dict_obj = None):
        """
        Check whether the value of this option is one of the allowed choices.
        
        This method is called automatically by validate()
        
        :param value: The value of this option.
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param dict_obj: The dict in which the value of this Option is stored. In most cases, the value of this option is evaluated simply as dict_obj[self.name]
        """
        for choice in self.choices:
            if value == choice:
                # Found a match, all ok.
                return value
            
            elif isinstance(value, str) and isinstance(choice, str) and value.lower() == choice.lower():
                # Found a match, but with different cAsInG, convert to the correct case.
                return choice
            
        # If we get here, there was no match.
        raise Disallowed_choice_exception(owning_obj, self, value)
        
            