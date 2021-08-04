#
# Classes that represent one or more configurables. Loaders construct real configurable objects from one or more dict objects on demand.
# Loaders are essentially a speed hack that means we don't have to construct 1000s of configurable objects each time we start up.
#

# General imports.
import deepmerge

# Silico imports.
from silico.config.configurable import Configurable
from silico.exception import Silico_exception

class Configurable_loader():
    """
    Abstract top-level class for configurable loaders.
    """
    
    @property
    def TAG(self):
        """
        An identifying tag.
        """
        return self.config.get('TAG', None)
    
    @property
    def TOP(self):
        """
        Whether this loader is a top level loader.
        
        TOP defaults to TRUE unless this loader is listed as a child of another loader, in which case it defaults to FALSE.
        This value can also be overriden by the value given in the config.
        """
        return self.config.get('TOP', self._TOP)
    
    @TOP.setter
    def TOP(self, value):
        self._TOP = value
    
    def __init__(self, file_name, TYPE, config, pseudo = False):
        """
        Default constructor for Configurable_loader objects.
        
        :param file_name: The path to the file from which we were parsed from. If not parsed from a file, file_name should be None.
        :param TYPE: The TYPE of the configurables we represent, a string identifying a class that is a parent to all the configurables we represent. Typically, Method, Program, Calculation or Basis_set.
        :param config: The config dict associated with this loader.
        :param pseudo: Whether this partial is a pseudo configurable. Pseudo configs act as placeholders only and won't appear as separate options in list etc (instead, the configurables in NEXT will).
        """
        # The file from which this loader was parsed from (if any).
        self.file_name = file_name
        
        # Save our TYPE.
        self.TYPE = TYPE
        
        # The config options at this node.
        self.config = config
        
        # Get our TYPE class.
        self.type_class = Configurable.from_class_handle(self.TYPE)
        
        # A list of the next loaders in the chain. For single loaders, len(self.NEXT) == 0.
        self.NEXT = []
        
        self._TOP = True
        self.pseudo = pseudo
    
    def resolve(self, index, *, parent_offset = 0, parent_config = None):
        """
        Resolve (one of) the configurable we represent, given by the unique ID/offest/index index.
        
        :param index: The index of the configurable to resolve. This value may be out of bounds, in which case an IndexException should be raised.
        :param parent_offset: The index of this loader in the parent loader, used when called recursively from hierarchical loaders.
        :param parent_config: The partially re-constructed config as processed by the previous loader in the chain, used when called recursively from hierarchical loaders.
        """
        raise NotImplementedError()
    
    def merge_with_parent(self, parent_config):
        """
        Merge the config options of this node with the config options of the parent node.
        
        :param parent_config: The config options from the parent node.
        """
        # First, merge our current parent object with ourself.
        deepmerge.always_merger.merge(parent_config, self.config)
        
        # Add our tag to the tag hierarchy.
        if not self.pseudo and self.TAG is not None:
            try:
                parent_config['TAG_HIERARCHY'].append(self.TAG)
                
            except KeyError:
                parent_config['TAG_HIERARCHY'] = [self.TAG]
                
    
    def configure(self, config):
        """
        Convert (or attempt to) a config dict to an appropriate configurable object.
        
        :raises Silico_exception: If the CLASS of the configurable is not set or cannot be found.
        :param config: The config dict.
        :returns: A loaded Configurable object.
        """
        config['TYPE'] = self.TYPE
        # These options have no meaning anymore.
        config.pop('SUB_TYPE')
        config.pop('TAG')
        config.pop('NEXT')
        
        # Try and get the configurable class.
        try:
            cls = self.type_class.from_class_handle(config['CLASS'])
        except ValueError:
            raise Silico_exception("Error loading configurable of type '{}' from file '{}'; CLASS '{}' is not recognised".format(self.TYPE, self.file_name, config['CLASS'])) from None
        except KeyError:
            #raise Silico_exception("Error loading configurable of type '{}' from file '{}'; no CLASS set".format(self.TYPE, config_path)) from None
            # If no class set, use the top level class.
            cls = self.type_class
        
        configurable = cls(self.file_name, False, **config) 
        configurable.configure_auto_name()
        configurable.validate()
        return configurable
    
    def size(self):
        """
        The recursive total number of child leaf nodes under this loader. For loaders than only represent one real configurable, size() == 1.
        """
        raise NotImplementedError()
    
    def link(self, loaders):
        """
        Link this partial loader to a number of other loaders.
        
        This implementation does nothing.
        """
        pass
    

class Partial_loader(Configurable_loader):
    """
    A configurable that is made up of multiple linked config objects.
    """
    
    def __init__(self, file_name, TYPE, config, pseudo = False):
        """
        :param file_name: The file from which this loader was parsed.
        :param config: The config options particular to this level of the partial configurable. They will also be set for all child (linked) objects, but may be overwritten.
        :param pseudo: Whether this partial is a pseudo configurable. Pseudo configs act as placeholders only and won't appear as separate options in list etc (instead, the configurables in NEXT will).
        """
        super().__init__(file_name, TYPE, config, pseudo = pseudo)
        
    def size(self):
        """
        The (recursive) total number of child leaf nodes. 
        """
        return  sum(child.size() for child in self.NEXT)
        
    def resolve(self, index, *, parent_offset = 0, parent_config = None):
        """
        Get a child configurable based on a unique index.
        
        Each configurable in the tree generated by this partial configurable can be identified by a unique index, corresponding to the position of that configurable if the tree was flattened in order.
        
        :raises IndexError: If index is out of range.
        :param index: The index to fetch.
        :param parent_offset: The index of this partial configurable in the parent configurable, used when called recursively.
        :param parent_config: The partially resolved dict from the parent configurable, used when called recursively.
        :returns: The resolved Configurable object.
        """
        parent_config = {} if parent_config is None else parent_config
        
        # First, merge our current parent object with ourself.
        self.merge_with_parent(parent_config)
                
        # We need to decide which child object from NEXT to continue to based on the range of possible indexes of each child.
        child_offset = 0
        for child in self.NEXT:
            child_size = child.size()
            if index < parent_offset + child_offset + child_size:
                # The config we want is in this child.
                return child.resolve(index, parent_offset = parent_offset + child_offset, parent_config = parent_config)
            
            # Add this child's size to the offset.
            child_offset += child_size
            
        # If we get this far, index is out of range.
        raise IndexError("Configurable index '{}' is out of range".format(index))

    def link(self, loaders):
        """
        Link this partial loader to a number of other loaders.
        
        This method will resolve the NEXT configs at this node.
        """
        for tag in self.config['NEXT']:
            # Look for a loader with this tag and add.
            # TODO: Should watch out for infinite recursion here.
            matching = [loader for loader in loaders if loader.TAG == tag]
            
            # Panic if we couldn't find any.
            if len(matching) == 0:
                raise Silico_exception("Error loading configurable of type '{}' from file '{}'; TAG '{}' could not be found".format(self.TYPE, self.file_name, tag))
            
            # Add all matching to our NEXT.
            for match in matching:
                # Unset the default TOP value.
                match.TOP = False
                self.NEXT.append(match)
            

class Configurable_list(Partial_loader):
    """
    A list of other configurables.
    
    Unlike partial configurables, configurable lists do not pass on or store any options themselves.
    """
    
    def __init__(self, configs, TYPE = None):
        """
        """
        super().__init__(file_name = None, TYPE = TYPE, config = {}, pseudo = False)        
        
        self.NEXT = configs
        
    
        


class Single_loader(Configurable_loader):
    """
    A configurable loader that represents a single configurable.
    """
    
    def __init__(self, file_name, TYPE, config):
        """
        Constructor for Single_loader objects.
        """
        super().__init__(file_name, TYPE, config)
        
    def size(self):
        """
        The size of this loader.
        
        Always returns 1.
        """
        return 1
    
    def resolve(self, index, *, parent_offset = 0, parent_config = None):
        """
        Get the configurable object that this loader represents.
        """
        parent_config = {} if parent_config is None else parent_config
        
        # First, merge our current parent object with ourself.
        self.merge_with_parent(parent_config)
        
        return self.configure(parent_config)
        
        
        
        
# class Pseudo_configurable():
#     """
#     A configurable that is a placeholder for a number of other configurables.
#     """