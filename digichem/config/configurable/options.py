# General imports.
from collections.abc import MutableMapping

# Silico imports.
from silico.config.configurable.option import Option
from silico.exception.base import Silico_exception
from silico.exception.configurable import Configurable_option_exception
from silico.config.configurable.base import Options_mixin


class Options_mapping(MutableMapping):
    """
    A class that 'binds' an Options object with a parent Configurable.
    """
    
    def __init__(self, options_obj, owner_obj, dict_obj):
        """
        Constructor for Options_binder objects.
        
        :param options_obj: An Options object to bind.
        :param owning_obj: A Configurable object which contains Options as a class attribute.
        :param dict_obj: The dict in which the dict corresponding to options_obj is stored.
        """
        self.options_obj = options_obj
        self.owning_obj = owner_obj
        self.dict_obj = dict_obj


    @property
    def sub_dict_obj(self):
        """
        Access to the dict object which holds values for the sub options of the Options object we map are stored.
        
        Options objects are essentially a dict in disguise. This makes access slightly more confusing to understand (for humans), because there are two dicts to consider:
          - The first is the 'normal' dict_obj in which the value of the Options object is stored (the 'value' here is the second dict). This dict_obj is probably the _configurable_options attribute of our owning object, but this is not guaranteed (for example, if this Options object is a child of another Options object).
          - The second is the dict in which our child options will store their values, which can be thought of as the dict our Options object is mapping. This second dict is what is referred to here as the sub_dict_obj.
        
        Access to the second, sub_dict_obj is provided through this property because there is no guarantee that it actually exists.
        """
        # TODO: This property can probably be replaced by the get_sub_dict() method of the Options class.
        try:
            return self.dict_obj[self.options_obj.name]
        
        except KeyError:
            if self.options_obj.name not in self.dict_obj:
                self.dict_obj[self.options_obj.name] = {}
                return self.dict_obj[self.options_obj.name]
            
            else:
                raise


    def get_sub_option(self, name):
        """
        Retrieve a sub option of the Options object we are mapping.
        """
        try:
            return self.options_obj.OPTIONS[name]
        
        except KeyError:
            # Unrecognised sub option?
            raise Configurable_option_exception(self.owning_obj, self.options_obj, "'{}' is not recognised as a valid sub option".format(name))
        
    def __len__(self):
        """
        The number of options in the Options object we are mapping.
        """
        return len(self.options_obj.OPTIONS)
    
    def __iter__(self):
        """
        Iteration magic method.
        """
        yield from self.options_obj.OPTIONS

    def __getitem__(self, key):
        """
        Fetch an item from the Options object we are mapping.
        
        The underlying Option will be evaluated to return its value, rather than returning the Option object itself.
        
        :param key: The name of an Option contained within the Options object that we are mapping.
        """
        return self.get_sub_option(key).get_from_dict(self.owning_obj, self.sub_dict_obj)


    def __setitem__(self, key, value):
        """
        Set the value of one of the Options contained in the Options object we are mapping.
        
        :param key: The name of the Option to set.
        :param value: The value of the Option to set.
        """
        self.get_sub_option(key).set_into_dict(self.owning_obj, self.sub_dict_obj, value)


    def __delitem__(self, key):
        """
        Reset the value of one of the Options contained in the Options object we are mapping to its default value.
        
        :param key: The key to delete.
        """
        self.get_sub_option(key).set_default(self.sub_dict_obj)


class Options(Option, Options_mixin):
    """
    A type of option that expects more options (another dict).
    """
    
    def __init__(self, *args, name = None, help = None, exclude = None, **kwargs):
        """
        """        
        # Use parent constructor.
        super().__init__(name = name, help = help, rawtype = dict, exclude = exclude)
        
        # Go through args and add to kwargs.
        for arg in args:
            if arg.name is None:
                raise Silico_exception("Configurable option given as positional argument to Options must have a name")
            kwargs[arg.name] = arg
        
        self.OPTIONS = {}
        
        # Set names of all kwargs.
        for argname in kwargs:
            kwargs[argname].name = argname
            self.OPTIONS[argname] = kwargs[argname]
            
            # Set ourselves as parent.
            kwargs[argname].add_parent(self)
            
    def add_parent(self, parent):
        """
        Add an owning parent Options object to this Option object.
        
        This method is called by the parent Options object when this Option is added to it.
        """
        #for child in self.OPTIONS:
        for child in self.OPTIONS.values():
            child.add_parent(parent)

            
    @property
    def OPTIONS(self):
        """
        A dictionary of child Option objects that we contain.
        
        This has to be accessed via a property because it is a property in the parent class
        """
        return self._OPTIONS


    @OPTIONS.setter
    def OPTIONS(self, value):
        self._OPTIONS = value


    def __get__(self, owning_obj, cls = None):
        """
        Get the values of this dictionary of Options.
        
        Options objects do not store any values themselves, instead a mapping object is returned (which supports dict[key] style access) which can return the true values of any contained Options objects.
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        """
        if owning_obj is None:
            return self
        
        return self.get_from_dict(owning_obj, owning_obj._configurable_options)


    def get_from_dict(self, owning_obj, dict_obj):
        """
        Get the values of this dictionary of Options from a specific dict.
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param dict_obj: The dict in which the value of this Option is stored. In most cases, the value of this option is evaluated simply as dict_obj[self.name]
        """
        return Options_mapping(self, owning_obj, dict_obj)


    def __set__(self, owning_obj, value):
        """
        Set the value of this Options object.
        
        Options objects cannot be set directly (because they are not actually dictionaries), but if the given value is iterable it will be used to set our options.
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param value: The new value to upate from. This should be a dict-like object that supports iteration via items().
        """
        self.set_into_dict(owning_obj, owning_obj._configurable_options, value)
        
    def get_sub_dict(self, dict_obj):
        """
        Get access to the dict object which holds values for our sub options.
        
        Options objects are essentially a dict in disguise. This makes access slightly more confusing to understand (for humans), because there are two dicts to consider:
          - The first is the 'normal' dict_obj in which the value of the Options object is stored (the 'value' here is the second dict). This dict_obj is probably the _configurable_options attribute of our owning object, but this is not guaranteed (for example, if this Options object is a child of another Options object).
          - The second is the dict in which our child options will store their values, which can be thought of as the dict our Options object is mapping. This second dict is what is referred to here as the sub_dict_obj.
        
        Access to the second, sub_dict_obj is provided through this function because there is no guarantee that it actually exists.
        """
        try:
            return dict_obj[self.name]
        
        except KeyError:
            if self.name not in dict_obj:
                dict_obj[self.name] = {}
                return dict_obj[self.name]
            
            else:
                raise


    def set_into_dict(self, owning_obj, dict_obj, value):
        """
        Set the value of this Options object.
        
        Options objects cannot be set directly (because they are not actually dictionaries), but if the given value is iterable it will be used to set our options.
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param dict_obj: The dict in which the value of this Option is stored. In most cases, the value of this option is evaluated simply as dict_obj[self.name]
        :param value: The new value to upate from. This should be a dict-like object that supports iteration via items().
        """
        try:
            for key, sub_value in value.items():
                try:
                    self.OPTIONS[key].set_into_dict(owning_obj, self.get_sub_dict(dict_obj), sub_value)
                
                except KeyError:
                    # Unknown sub option?
                    if key not in self.OPTIONS:
                        raise Configurable_option_exception(owning_obj, self, "Cannot update values from iterable, sub option '{}' is not recognised".format(key))
                    
                    else:
                        raise
            
        except Exception:
            # Check to see if value has an items() method.
            if not callable(getattr(value, "items", None)):
                raise Configurable_option_exception(owning_obj, self, "Cannot update values; the given value does not support the items() method")
            
            else:
                raise


    def set_default(self, owning_obj, dict_obj):
        """
        Reset this option to default.
        
        Note that this method is an alias for calling del() on this attribute.
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param dict_obj: The dict in which the value of this Option is stored. In most cases, the value of this option is evaluated simply as dict_obj[self.name]
        """
        for sub_option in self.OPTIONS.values():
            sub_option.set_default(owning_obj, self.get_sub_dict(dict_obj))


    def is_default(self, dict_obj):
        """
        Whether the value of this option is currently the default or not.
        
        Options objects are considered as default only if all their children are default.
        
        :param dict_obj: The dict in which the value of this Option is stored. In most cases, the value of this option is evaluated simply as dict_obj[self.name]
        """
        # This is safe, because all([]) == True.
        return all([sub_option.is_default(self.get_sub_dict(dict_obj)) for sub_option in self.OPTIONS.values()])


    def validate(self, owning_obj, dict_obj = None):
        """
        Validate the options contained within this Options object.
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param dict_obj: The dict in which the value of this Option is stored.
        """
        if dict_obj is None:
            dict_obj = owning_obj._configurable_options
        
        # Our children will find their values in a different dict to where we find ourself.
        sub_dict_obj = dict_obj.get(self.name, {})
#         try:
#             sub_dict_obj = dict_obj[self.name]
#             
#         except KeyError:
#             print("blagh")
#             raise   
        
        # Validate each of our sub options.
        self.validate_children(owning_obj, sub_dict_obj)        

