import re

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
    
    def validate_children(self, owning_obj, dict_obj):
        """
        Validate the child Options of this object.
        
        :param owning_obj: The owning object which contains these Options.
        :param dict_obj: The dict in which the values of the child Options are stored.
        """
        # First, validate each of our known options.
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


class Configurable(Dynamic_parent, Options_mixin):
    """
    Class that represents a Configurable.
    
    Configurables are classes that specify how certain object attributes should be set by defining a number of Option objects as class attributes.
    Each Option object maps a certain attribute on the owning configurable object and defines, for example, an allowed type, a default value, a help string, a list of allowed values etc.
    """
    
    def __init__(self, file_name = None, validate_now = False, **kwargs):
        """
        Constructor for Configurable objects.
        
        :param file_name: If this configurable was loaded from a file, the full file path to that file.
        :param relative_file_name: If this configurable was loaded from a file, the file path to that file relative to the top directory for storing config files.
        :param validate_now: If True, the given options will be validated before this constructor returns. Validation can also be performed at any time by calling validate(). 
        """
        self._configurable_options = {}
        self.inner_cls = None
        self.file_name = file_name
        
        # Set all our configurable options.
        # Setting like this might be unsafe because we're not deep copying...
        self._configurable_options = kwargs
        
        # If we've been asked to, validate.
        if validate_now:
            self.validate()
        
    # Configurable options.
    NAME = Option(help = "The unique name of this Configurable", type = str, required = True)
    ALIASES = Option(help = "A list of alternative names for this Configurable", default = [], type = list, no_edit = True)
    TAG_HIERARCHY = Option(help = "A hierarchical list of tags that were combined to form this configurable", default = [], type = list, no_edit = True)
    CATEGORY = Option(help = "A hierarchical list of categories that this configurable belongs to", default = [], type = list)
    AUTO_CATEGORY = Option(help = "If no category is explicitly set and this option is true, the category will be automatically set from the path from which this configurable was loaded", type = bool, default = True)
    UPDATE_CATEGORY = Option(help = "The name of a category to update. If given, this configurable will be used to override the options set by all other configurable that have UPDATE_CATEGORY as one of their categories", type = str, default = None, no_edit = True)
    TYPE = Option(help = "The type of this Configurable", choices = ('method', 'program', 'calculation', 'basis_set'), required = True, default = 'calculation', type = str, no_edit = True)
    CLASS = Option(
        help = "The specific sub-type of this Configurable",
        choices = lambda option, configurable: [handle for cls in Configurable.from_class_handle(configurable.TYPE).recursive_subclasses() if hasattr(cls, 'CLASS_HANDLE') for handle in cls.CLASS_HANDLE],
        required = True,
        no_edit = True
    )
    hidden = Option(help = "If True, this configurable will not appear in lists (but can still be specified by the user). Useful for configurables that should not be used naively", type = bool, default = False)
    warning = Option(help = "A warning message to display when this configurable is chosen", default = None, type = str)

    # A regex used to extract useful parts of a directory name.
    DIR_TO_CATEGORY_REGEX = re.compile(r"[0-9]* (.*)")

    @classmethod
    def dir_name_to_category(self, dir_name):
        """
        Convert a directory name to an appropriate category name.
        """
        match = self.DIR_TO_CATEGORY_REGEX.match(dir_name)
        if match:
            return match.group(1)
        else:
            return dir_name
    
    def configure_auto_name(self):
        """
        Setup automatic names for this configurable.
         
        Automatic names are ones that are generated automatically from some other information about the configurable.
        Currently, this includes the CATEGORY and NAME attributes, which can be set automatically based on the file path from which the configurable was loaded. 
        """
        # First, set CATEGORY, as we may need this later for our NAME.
        if len(self.CATEGORY) == 0 and self.AUTO_CATEGORY and self.TAG_HIERARCHY is not None:
            # No category and we're able to set one.
            # We will modify the folder names we're going to use to build our category in the following ways:
            # - Numerical characters from the start of each category name will be stripped (as these are often used to set ordering).
            #self.CATEGORY = [self.dir_name_to_category(category_name) for category_name in self.relative_file_path.parts]
            self.CATEGORY = self.TAG_HIERARCHY
             
        # Now set our NAME if empty.
        if not hasattr(self, "NAME") and len(self.CATEGORY) > 0:
            # No NAME, set from the category.
            self.NAME = " ".join(self.CATEGORY)

    def validate(self):
        """
        Check that all the configurable options of this configurable have been set appropriately.
        
        :raises Exception: If one of the Options of this configurable is invalid.
        """
        self.validate_children(self, self._configurable_options)

    
    def match(self, identifier):
        """
        Determine whether a given identifier matches this Configurable object.
        
        Identifiers match (case insensitively) to either NAME or any of ALIASES.
        
        :param identifier: A string to match against.
        :returns: True or False
        """
        return identifier.lower() in [name.lower() for name in self.NAMES]
            
    @property
    def NAMES(self):
        """
        A list of all the unique names this Configurable is known by (NAME + ALIASES).
        """
        names = []
        
        # Add NAME if we have one.
        if self.NAME is not None:
            names.append(self.NAME)
        
        # Add any ALIAS.
        names.extend(self.ALIASES)
        
        # Done
        return names
    
    @property
    def description(self):
        """
        A string that describes this Configurable object.
        """
        # Start with the name.
        try:
            desc = self.NAME
        except AttributeError:
            # No name set, this is pretty unusual but try and continue.
            desc = "NO-NAME"
        
        # Add alias if present.
        if len(self.ALIASES) > 0:
            desc += " ({})".format(" ,".join(self.ALIASES))
            
        return desc
                    
    def grouped(self, grouped_dict, *, names = None, number):
        # TODO: Unused?
        names = self.GROUP if names is None else names
        
        # This is the name we will be adding to in grouped_dict.
        # If it is None, we will be adding to a list (we have no more names to go through).
        curname = names[0] if len(names) > 0 else None
            
        if curname is None:
            # Create a new list if one doesn't exist.
            if curname not in grouped_dict:
                grouped_dict[curname] = []
            
            # Add to list
            grouped_dict[curname].append((number, self))
        else:
            # We have more names to work through.
            # Create an empty dict under the name if not existing.
            if curname not in grouped_dict:
                grouped_dict[curname] = {}
                
            # Recurse.
            self.group(grouped_dict[curname], names[1:])
            
        return grouped_dict
        
    
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
        pass
            
        
        