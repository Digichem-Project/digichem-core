# General imports.
from collections.abc import MutableMapping
import deepmerge

# Silico imports.
from silico.config.configurable.option import Option, InheritedAttrError
from silico.exception.base import Silico_exception
from silico.exception.configurable import Configurable_option_exception,\
    Configurable_exception
from silico.misc import Default


class Options_mixin():
    """
    Mixin class for those that contain configurable options.
    """
    
    @classmethod
    def get_cls_option(cls, name):
        """
        Get a specific sub option.
        """
        try:
            return cls.get_cls_options()[name]
        
        except KeyError:
            raise ValueError("No class option with name '{}'".format(name))
    
    @classmethod
    def get_cls_options(cls):
        """
        Get a dict of all the configurable options that belong directly to (are defined on) this class.
        
        Unlike get_options(), this method will not attempt to resolve inheritance.
        
        The key of each item is the name of the corresponding option.
        """
        # TODO: This function gets called quite a lot (anecdotally), might be worth optimising (caching perhaps?).
        # We apply a custom sort here which ignores case (otherwise all uppercase options appear before all lowercase which is annoying).
        return dict(sorted({getattr(cls, attr).name: getattr(cls, attr) for attr in dir(cls) if isinstance(getattr(cls, attr), Option)}.items(), key = lambda v: v[0].upper()))
    
    def get_options(self, owning_obj = None):
        """
        Get a dict of all the configurable options of this object.
        
        The key of each item is the name of the corresponding option.
        """
        return self.get_cls_options()
    
    def prune(self, owning_obj, dict_obj):
        """
        Prune none and empty values, removing them.
        """
        for key, value in tuple(dict_obj.items()):
            # First, if the value is another dict, prune that.
            if isinstance(value, dict):
                self.prune(owning_obj, value)
                
                # Once we've pruned, check the length of the dict.
                if len(value) == 0:
                    # This dict is empty and can be removed.
                    del(dict_obj[key])
                    
            # If None, delete.
            #if value is None:
            #    del(dict_obj[key])
    
    def validate_children(self, owning_obj, dict_obj):
        """
        Validate the child Options of this object.
        
        :param owning_obj: The owning object which contains these Options.
        :param dict_obj: The dict in which the values of the child Options are stored.
        """
        # First, prune empty values.
        self.prune(owning_obj, dict_obj)
        
        options = self.get_options(owning_obj)
        
        # Next, validate each of our known options.
        for option in options.values():
            
            option.validate(owning_obj, dict_obj)
        
        # Next, check for any exclusions.
        # NOTE: We have to do this in a different loop to the one where we call validate() above.
        # This is because validate() can set an option to default, which can change whether
        # the option is considered to be clashing or not.
        for option in options.values():
            # Also check for exclusions.
            # TODO: Extend exclusions so they can support nested options.
            for exclusion in option.exclude:
                # This option has an exclusion, check at least one of it and the exclusion is not set.
                try:
                    if not options[exclusion].is_default(owning_obj, dict_obj) and not option.is_default(owning_obj, dict_obj):
                        raise Configurable_exception(owning_obj, "options '{}' and '{}' cannot be set at the same time (mutually exclusive)".format(option.name, exclusion))
                
                except KeyError:
                    # One of the given options cannot be found.
                    raise Configurable_option_exception(owning_obj, option, "The option '{}' in exclude cannot be found".format(exclusion)) from None
            
        # We also need to make sure there are no unexpected options.
        for unexpected_key in set(dict_obj).difference(options):
            # Although this looks like a loop, we will obviously only raise the first exception.
            msg = "unrecognised option '{}' with value '{}'".format(unexpected_key, dict_obj[unexpected_key])
            if hasattr(self, "is_configurable"):
                raise Configurable_exception(owning_obj, msg)
            
            else:
                raise Configurable_option_exception(owning_obj, self, msg)


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
            return self.options_obj.get_options(self.owning_obj)[name]
        
        except KeyError:                
            # Give up and panic.
            raise Configurable_option_exception(self.owning_obj, self.options_obj, "'{}' is not recognised as a valid sub option".format(name))
        
    def __len__(self):
        """
        The number of options in the Options object we are mapping.
        """
        return len(self.options_obj.get_options(self.owning_obj))
    
    def __iter__(self):
        """
        Iteration magic method.
        """
        yield from self.options_obj.get_options(self.owning_obj)

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
    
    def __init__(self, *args, name = None, help = Default(None), validate = Default(None), exclude = Default(None), no_edit = Default(False), **kwargs):
        """
        """
        consumed_kwargs = {
            "name": name,
            "help": help,
            "validate" : validate,
            "exclude": exclude,
            "no_edit": no_edit,
        }
        # Check that none of the arguments consumed by this constructor are Option objects.
        # This can be an easy mistake, where an Options object is created with a sub option
        # with the same name as one of our arguments ("name", "help", "exclude", "no_edit" etc).
        for kwarg_name, kwarg_value in consumed_kwargs.items():
            if isinstance(kwarg_value, Option):
                raise TypeError("Constructor keyword argument '{}' cannot be an Option object.".format(kwarg_name))
        
        # Certain constructor arguments can be inherited from a parent Options object if they're not given.
        # Because of this, we check whether they are in kwargs rather than specifying them explicitly, and
        # set a flag not to inherit them if they're not given.
        # Args to inherit from parent.
        self._inherit = []
        
        #TODO: Can we not use the inheritance mechanism in Option?
        for arg in ("help", "validate", "exclude", "no_edit"):
            if isinstance(consumed_kwargs[arg], Default):
                self._inherit.append(arg)
        
        # Use parent constructor.
        super().__init__(name = name, help = help, validate = validate, exclude = exclude, no_edit = no_edit)
        
        # Go through args and add to kwargs.
        for arg in args:
            if arg.name is None:
                raise Silico_exception("Configurable option given as positional argument to Options must have a name")
            kwargs[arg.name] = arg
        
        # Dict of child options.
        self._options = {}
        
        # Set names of all kwargs.
        for argname in kwargs:
            kwargs[argname].name = argname
            self._options[argname] = kwargs[argname]
            
            # Set ourselves as parent.
            kwargs[argname].add_parent(self)
    
    
    def __set_name__(self, owning_cls, name):
        """
        Called automatically during class creation, allows us to know the attribute name we are stored under.
        """
        super().__set_name__(owning_cls, name)
        
        # Also call set name on our children.
        # This is not particularly useful for setting the name, but it is for passing the owning class name to children.
        for child_name, child in self._options.items():
            child.__set_name__(owning_cls, child_name)
    
    
    @property
    def num_child_options(self):
        """
        The number of child/sub options contained within this one. For 'normal' options, this is always 0.
        """
        # WARNING: This only counts direct children of this Options.
        return len(self._options)
    
            
    def add_parent(self, parent):
        """
        Add an owning parent Options object to this Option object.
        
        This method is called by the parent Options object when this Option is added to it.
        """
        super().add_parent(parent)
        #for child in self.OPTIONS:
        for child in self._options.values():
            child.add_parent(parent)
    

    def get_inherited_options(self, owning_obj, _mro = None):
        """
        Build a dictionary of child Option objects that we contain, including those we inherit from any base classes of our owning class.
        """
        # Get our base options.
        try:
            parent_obj, _mro = self.get_base_option(type(owning_obj), _mro)
            parent_options = parent_obj.get_inherited_options(owning_obj, _mro)
            
        except InheritedAttrError:
            parent_options = {}
            
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
        options = self.get_options(owning_obj)
        try:
            for key, sub_value in value.items():
                try:
                    options[key].set_into_dict(owning_obj, self.get_sub_dict(dict_obj), sub_value)
                
                except KeyError:
                    # Unknown sub option?
                    if key not in options:
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
        
        for option in self.get_options(owning_obj).values():
            if explicit or not option.is_default(owning_obj, self.get_sub_dict(dict_obj)):
                dump[option.name] = option.dump(owning_obj, self.get_sub_dict(dict_obj), explicit = explicit)
                
        return dump

    def set_default(self, owning_obj, dict_obj):
        """
        Reset this option to default.
        
        Note that this method is an alias for calling del() on this attribute.
        
        :param owning_obj: The owning object on which this Option object is set as a class attribute.
        :param dict_obj: The dict in which the value of this Option is stored. In most cases, the value of this option is evaluated simply as dict_obj[self.name]
        """
        for sub_option in self.get_options(owning_obj).values():
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
        
        # Validate each of our sub options.
        self.validate_children(owning_obj, sub_dict_obj)
        
        # Finally, call our custom validate function if given.
        value = self.get_from_dict(owning_obj, dict_obj)
        if not self._validate(self, owning_obj, value):
            # Invalid.
            # This won't get raised often, because most validation functions will throw their own exceptions.
            raise Configurable_option_exception(owning_obj, self, "Validation for Options object failed")

