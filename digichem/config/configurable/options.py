# General imports.
from collections.abc import MutableMapping
from itertools import chain
import deepmerge

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
    
    def get_sub_option(self, name, _mro = None):
        """
        Retrieve a sub option of the Options object we are mapping.
        """
        try:
            return self.options_obj.OPTIONS[name]
        
        except KeyError:                
            # Give up and panic.
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
        self.get_sub_option(key).set_default(self.owning_obj, self.sub_dict_obj)


class Options(Option, Options_mixin):
    """
    A type of option that expects more options (another dict).
    """
    
    def __init__(self, *args, name = None, help = None, exclude = None, no_edit = False, **kwargs):
        """
        """        
        # Use parent constructor.
        super().__init__(name = name, help = help, exclude = exclude, no_edit = no_edit)
        
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
            
    @property
    def num_child_options(self):
        """
        The number of child/sub options contained within this one. For 'normal' options, this is always 0.
        """
        return len(self.OPTIONS)
            
    def add_parent(self, parent):
        """
        Add an owning parent Options object to this Option object.
        
        This method is called by the parent Options object when this Option is added to it.
        """
        super().add_parent(parent)
        #for child in self.OPTIONS:
        for child in self.OPTIONS.values():
            child.add_parent(parent)

    def get_inherited_options(self, owning_obj, _mro = None):
        """
        Build a dictionary of child Option objects that we contain, including those we inherit from any base classes of our owning class.
        """
        # The full 'access' path to this option.
        # This will consist of an attribute access as the first item,
        # eg: obj.option
        # followed by a number of dict like accesses to a specific option,
        # eg: obj.option['sub1']['sub2']
        resolve_path = [part.name for part in chain(self.parents, (self,))]
        
        # The 'parent' object of our owning object.
        # Decide which parent class to look at.
        # We do this by walking up the method resolution order of our owning_obj,
        # which we keep track of via the recursive argument _mro.
        if _mro is None:
            _mro = list(type(owning_obj).__mro__[1:])
        
        # Here, we are not actually interested in the value of the option,
        # but rather the Option(s) object itself.
        # Hence we access via the base class itself, rather than using super().
        parent_options = None
        while parent_options is None:
            parent_cls = _mro.pop(0)
            
            try:
                current = getattr(parent_cls, resolve_path[0])
            
            except AttributeError:
                if not hasattr(parent_cls, resolve_path[0]):
                    # The parent object doesn't have any configurable options.
                    parent_options = {}
                    break
                else:
                    raise
                
            # Walk up the nested Options object to find ourself.
            try:
                for resolve_part in resolve_path[1:]:
                    current = current._options[resolve_part]
            
            except KeyError:
                # Couldn't walk all the way up the path, try again from the next base class.
                # Doing nothing will loop us around again.
                pass
            
            else:
                # Got our option, stop for now.
                # Get the inherited options of our parent Options object.
                parent_options = current.get_inherited_options(owning_obj, _mro)
        
        
        # Merge the returned options with our own.
        return deepmerge.always_merger.merge(parent_options, self._options)
    
    def get_options(self, owning_obj):
        """
        Get a dict of all the configurable options of this object.
        
        The key of each item is the name of the corresponding option.
        """
        return self.get_inherited_options(owning_obj)
    

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
            
    def dump(self, owning_obj, dict_obj, explicit = False):
        """
        Dump the value of this option so it can be serialised (for example, to yaml).
        
        :param explicit: If True, all values will be dumped. If False, only non-default values will be dumped.
        :returns: A dumped version of this option's value.
        """
        dump = {}
        
        for option in self.OPTIONS.values():
            if explicit or not option.is_default(self.get_sub_dict(dict_obj)):
                dump[option.name] = option.dump(owning_obj, self.get_sub_dict(dict_obj), explicit = explicit)
                
        return dump

    def set_default(self, owning_obj, dict_obj):
        """
        Reset this option to default.
        
        Note that this method is an alias for calling del() on this attribute.
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param dict_obj: The dict in which the value of this Option is stored. In most cases, the value of this option is evaluated simply as dict_obj[self.name]
        """
        for sub_option in self.OPTIONS.values():
            sub_option.set_default(owning_obj, self.get_sub_dict(dict_obj))


    def is_default(self, owning_obj, dict_obj):
        """
        Whether the value of this option is currently the default or not.
        
        Options objects are considered as default only if all their children are default.
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param dict_obj: The dict in which the value of this Option is stored. In most cases, the value of this option is evaluated simply as dict_obj[self.name]
        """
        # This is safe, because all([]) == True.
        return all([sub_option.is_default(owning_obj, self.get_sub_dict(dict_obj)) for sub_option in self.get_options(owning_obj).values()])


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

