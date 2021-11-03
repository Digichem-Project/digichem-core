from silico.exception import Configurable_exception
from silico.misc import Dynamic_parent
from silico.config.configurable.option import Option


class Options_mixin():
    """
    Mixin class for those that contain configurable options.
    """
    
    @property
    def OPTIONS(self):
        """
        Get a list of all Configurable Options of this object.
        """
        return {getattr(type(self), attr).name: getattr(type(self), attr) for attr in dir(type(self)) if isinstance(getattr(type(self), attr), Option)}
    
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
            if value is None:
                del(dict_obj[key])
    
    def validate_children(self, owning_obj, dict_obj):
        """
        Validate the child Options of this object.
        
        :param owning_obj: The owning object which contains these Options.
        :param dict_obj: The dict in which the values of the child Options are stored.
        """
        # First, prune empty values.
        self.prune(owning_obj, dict_obj)
        
        # Next, validate each of our known options.
        for option in self.OPTIONS.values():
            option.validate(owning_obj, dict_obj)
            
            # Also check for exclusions.
            for exclusion in option.exclude:
                # This option has an exclusion, check at least one of it and the exclusion is not set.
                if not self.OPTIONS[exclusion].is_default(dict_obj) and not option.is_default(dict_obj):
                    raise Configurable_exception(self, "options '{}' and '{}' cannot be set at the same time (mutually exclusive)".format(option.name, exclusion))
            
        # We also need to make sure there are no unexpected options.
        for unexpected_key in set(dict_obj).difference(self.OPTIONS):
            # Although this looks like a loop, we will obviously only raise the first exception.
            raise Configurable_exception(owning_obj, "unrecognised option '{}'".format(unexpected_key))


class Configurable(Options_mixin):
    """
    Class that represents a Configurable.
    
    Configurables are classes that specify how certain object attributes should be set by defining a number of Option objects as class attributes.
    Each Option object maps a certain attribute on the owning configurable object and defines, for example, an allowed type, a default value, a help string, a list of allowed values etc.
    """
    
    def __new__(cls, *args, validate_now = True, **kwargs):
        instance = super().__new__(cls)
        
        # This is where the actual values for configurable options are stored.
        # We set this here so configurable options can be used in __init__ (particularly by subclasses).
        instance._configurable_options = {}
        
        # Set all our configurable options.
        # Look through kwargs for options that we recognise.
        values = {}
        for option in instance.OPTIONS.values():
            try:
                values[option.name] = kwargs[option.name]
            
            except KeyError:
                pass
            
        # Setting like this might be unsafe because we're not deep copying...
        instance._configurable_options = values
        
        return instance
    
    def __init__(self, validate_now = True, **kwargs):
        """
        Constructor for Configurable objects.
        
        :param validate_now: If True, the given options will be validated before this constructor returns. Validation can also be performed at any time by calling validate().
        :param **kwargs: Initial values for the Options of this configurable.
        """
        # If we've been asked to, validate.
        if validate_now:
            self.validate()

    def validate(self):
        """
        Check that all the configurable options of this configurable have been set appropriately.
        
        :raises Exception: If one of the Options of this configurable is invalid.
        """
        self.validate_children(self, self._configurable_options)
    
    @property    
    def description(self):
        """
        """
        return str(type(self))


class Configurable_class_target(Dynamic_parent, Configurable):
    """
    A configurable object which specifies which type of class it is.
    """
        
    # Configurable options.
    #TYPE = Option(help = "The type of this Configurable", choices = ('destination', 'program', 'calculation', 'basis_set'), required = True, default = 'calculation', type = str, no_edit = True)
    TYPE = Option(help = "The parent class of this target, the class we will be replaced as will be a child class of this.", required = True, type = str, no_edit = True)
    class_name = Option(help = "The name of a class that we will be replaced as.", required = True, no_edit = True)
    name = Option(help = "The unique name of this configurable target", type = str, required = True)


    def __init__(self, loader_list = None, validate_now = True, **kwargs):
        """
        Constructor for Configurable objects.
        
        :param loader_list: If this configurable was loaded from a (number of) configurable loaders, an ordered list of those loaders. 
        :param validate_now: If True, the given options will be validated before this constructor returns. Validation can also be performed at any time by calling validate(). 
        """
        self.inner_cls = None
        self.loader_list = loader_list if loader_list is not None else []
        
        self.configure_auto_name()
        
        Configurable.__init__(self, validate_now =validate_now, **kwargs)
    
    @property
    def file_names(self):
        """
        A list of files from which this configurable was loaded.
        """
        return [loader.file_name for loader in self.loader_list if loader.file_name is not None]
        
    @property
    def file_name(self):
        """
        The file names from which this configurable was loaded from, formatted as a single string (will return None if not loaded from any files).
        """
        if len(self.file_names) == 0:
            return None
        
        elif len(self.file_names) == 1:
            return self.file_names[0]
        
        else:
            return "\n".join(self.file_names)
        
    @property
    def tag_hierarchy(self):
        """
        An ordered list of the TAGs of each of the configurable loaders that were combined to generate this configurable.
        """
        return [loader.ALIAS for loader in self.loader_list if not loader.pseudo and loader.ALIAS is not None]
    
    def index(self):
        """
        Get the index of this configurable.
        
        Note that a configurable can only have an index if it was loaded from a (number of) configurable loaders, otherwise this method will throw an index error.
        """
        return self.loader_list[0].index_of_path(self.loader_list)
    
    def configure_auto_name(self):
        """
        Setup automatic names for this configurable.
         
        Automatic names are ones that are generated automatically from some other information about the configurable.
        """
        tag_hierarchy = self.tag_hierarchy
        # If our name is empty, set from our tag hierarchy.
        if not hasattr(self, "name") and len(tag_hierarchy) > 0:
            # No name, set from the category.
            self.name = " ".join(tag_hierarchy)
    
    @property
    def description(self):
        """
        A string that describes this Configurable object.
        """
        # Get our ID if we've got one.
        try:
            desc = "[{}] ".format(self.index())
        except Exception:
            desc = ""
        
        # Start with the name.
        try:
            desc += self.name
        except AttributeError:
            # No name set, this is pretty unusual but try and continue.
            desc += "NO-NAME"
            
        # Add our type.
        try:
            desc += " ({})".format(self.TYPE)
        except AttributeError:
            pass
            
        return desc
                                
    
    ############################
    # Class creation mechanism #
    ############################
    
    def classify(self):
        """
        Create a new class from this Configurable object.
        
        The new class will inherit this object's attributes as class-level attributes, so all new objects created from the class will 'share' the attributes of this object.
        """
        cls = type(type(self).__name__ + "_actual", (self._actual, type(self)), vars(self))
        cls.__module__ = '__main__'
        return cls
    
    def finalize(self, force = True):
        """
        Finalize this configurable, indicating that no further changes are going to be made.
        
        Calling finalize creates a new class which is stored at this object's inner_cls attribute, which can be used to create new objects, using this object as a template.
        Alternatively, once finalize() has been called, new objects from this inner class can be created by calling this object.
        
        Once called, all objects created from this class template will have the same type and share class-level attributes.
        Further modifications to this object will not be reflected in children objects, unless finalize is called again in which case changes will propagate to newly created children only.
        Objects created before and after a call to finalize() will not share the same type (although they will have the same type name...)
        
        :param force: Whether to finalize again if this method has already been called. If False and finalize() has been previously called, nothing will happen.
        """
        if force or self.inner_cls is None:
            self.inner_cls = self.classify()
        
    def __call__(self, *args, **kwargs):
        """
        Create a new object using this object as a class template.
        
        :raises TypeError: If finalize has not yet been called.
        """
        try:
            return self.inner_cls(*args, **kwargs)
        except TypeError:
            # Type Error, possibly because inner_cls is None.
            if self.inner_cls is None:
                raise TypeError("Attempted to create object before calling finalize() (inner_cls is None)")
            else:
                raise
        
    class _actual():
        """
        The inner class that is used as the base for the classes returned by classify().
        
        Objects of this class will have attributes of the outer-level object as class-level attributes.
        
        This default implementation does nothing.
        """
        
        
        def __new__(cls, *args, **kwargs):
            return object.__new__(cls)
            
        
        