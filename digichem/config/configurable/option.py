from silico.exception.configurable import Configurable_option_exception,\
    Missing_option_exception

class Option():
    """
    Class for specifying an option in a configurable.
    
    Options are descriptors that perform type checking and other functionality for Configurables; they expose the options that a certain configurable expects.
    """
    
    def __init__(self, name = None, *, default = None, help = None, choices = None, validate = None, type = None, rawtype = None, exclude = None, required = False, no_edit = False):
        """
        Constructor for Configurable Option objects.
        
        :param name: The name of this option. If None is given this will be determined automatically from the name of the attribute this option is stored under.
        :param default: Default value for this option. Alternatively, default can be a callable which will be called with 2 arguments: this Option object and the owning Configurable object and should return the default value.
        :param help: Descriptive help string
        :param choices: An iterable of valid choices for this option. Alternatively, choices can be a callable which will be called with 2 arguments: this Option object and the owning Configurable object and should return the list of options.
        :param validate: Function called to check that the given value is valid. The function will be called with 3 arguments: this Option object, the owning Configurable object and the value being set, and should return True or False as appropriate.
        :param type: A callable that is used to set the type of value.
        :param exclude: A list of strings of the names of attributes that this option is mutually exclusive with.
        :param required: Whether this option is required or not.
        :param no_edit: Flag to indicate that this option shouldn't be edited.
        """
        self.name = name
        self._default = default
        self.type = type
        self.rawtype = type if rawtype is None else rawtype
        self.help = help
        self._choices = choices
        self._validate = validate if validate is not None else lambda option, configurable, value: True
        self.exclude = exclude if exclude is not None else []
        self.required = required
        self.no_edit = no_edit
        
        # By definition, Options that are required can have no default, so we'll delete this attribute.
        if self.required:
            del(self._default)


    def choices(self, owning_obj):
        """
        Get the list of allowed options for this option.
        
        This property will evaluate self._choices if it is a callable.
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :return: The list of choices, or None if no choices.
        """
        if not callable(self._choices):
            return self._choices
        else:
            return self._choices(self, owning_obj)


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
                raise Missing_option_exception(owning_obj, self.name)


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
            try:
                value = self.type(value) if self.type is not None else value
                self.set_into_dict(owning_obj, dict_obj, value)
            except (TypeError, ValueError) as e:
                raise Configurable_option_exception(owning_obj, self, "value '{}' of type '{}' is of invalid type".format(value, type(value).__name__)) from e
        
        # If we have a list of options, check we chose one.
        choices = self.choices(owning_obj)
        
        if choices is not None:
            # If we are a list type, we'll check each item in value (rather than value itself).
            try:
                if issubclass(self.type, list) or issubclass(self.type, tuple):
                    values = value
                else:
                    values = [value]
            except TypeError:
                # issubclass raises this all the time...
                values = [value]
            
            for subvalue in values:
                if subvalue not in choices:
                    raise Configurable_option_exception(owning_obj, self, "value '{}' is not one of the allowed choices".format(subvalue))
            
        # Check the value is valid.
        if not self._validate(self, owning_obj, value):
            # Invalid.
            raise Configurable_option_exception(owning_obj, self, "value '{}' of type '{}' is invalid".format(value, type(value).__name__))
            