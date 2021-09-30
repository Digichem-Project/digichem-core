#
# Classes that represent one or more configurables. Loaders construct real configurable objects from one or more dict objects on demand.
# Loaders are essentially a speed hack that means we don't have to construct 1000s of configurable objects each time we start up.
#

# General imports.
import deepmerge

# Silico imports.
from silico.config.configurable import Configurable
from silico.exception.configurable import Configurable_loader_exception,\
    Tag_path_length_error, Unresolvable_tag_path_error
from silico.exception.base import Silico_exception
from pkg_resources._vendor.pyparsing import ident

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
        :param TYPE: The TYPE of the configurables we represent, a string identifying a class that is a parent to all the configurables we represent. Typically, Destination, Program, Calculation or Basis_set.
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
        
        
    def __iter__(self):
        """
        Iteration magic method
        """
        pos = 1
        while True:
            try:
                yield self.resolve(pos, validate = False)
                
            except IndexError:
                # We're all done.
                return
            
            pos +=1        
    
    def get_concrete_children(self, _concrete_children = None, _current_path = None):
        """
        Get a list of the children of this loader that are concrete.
        
        Concrete here means that the child would appear as a node in a list of options, ie is not pseudo or similar.
        This method is recursive, any direct children that are not concrete (they are pseudo) will also have get_concrete_children() called.
        
        All parameters to this method are only used when called recursively; they should not normally be specified by the user.
        
        :param _concrete_children: The current list of concrete children.
        :param _current_path: A list of pseudo loaders that have already been traversed up to this point.
        :returns: A list of paths to concrete children. Each 'path' is itself a list of loaders which if traversed will lead to the concrete child.
                  Only the last element in the list will be a non-pseudo loader, while all others will be a pseudo loader.
                  If a direct child of this loader is concrete, then the 'path' will be a list containing a single element (which will be that direct child).
                  
        """
        concrete_children = [] if _concrete_children is None else _concrete_children
        if _current_path is None:
            current_path = []
            
        else:
            # We make a new list from the existing current path because it will be modified differently each time this method is called recursively.
            current_path = list(_current_path)
            
            # Add ourself to current path.
            current_path.append(self)
        
        # Iterate through each of our immediate children.
        for child in self.NEXT:
            if child.pseudo:
                # This child is pseudo (ie, not concrete), so we want its children.
                child.get_concrete_children(concrete_children, current_path)
                
            else:
                # This child is concrete, add it to the list.
                current_path.append(child)
                # Add the path to our total list.
                concrete_children.append(current_path)
                
        return concrete_children
    
    def resolve(self, *args, **kwargs):
        raise NotImplementedError()
    
    def resolve_path(self, *args, **kwargs):
        raise NotImplementedError()
    
    def path_by_index(self, *args, **kwargs):
        raise NotImplementedError()
    
    def path_by_tags(self, *args, **kwargs):
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
                
        # Add our filename to the file hierarchy (unless it's None).
        if self.file_name is not None:
            try:
                parent_config['FILE_HIERARCHY'].append((self.TAG, self.file_name))
                
            except KeyError:
                parent_config['FILE_HIERARCHY'] = [(self.TAG, self.file_name)]                
    
    def configure(self, config, validate = True):
        """
        Convert (or attempt to) a config dict to an appropriate configurable object.
        
        :raises Silico_exception: If the class_name of the configurable is not set or cannot be found.
        :param config: The config dict.
        :param validate: Whether to validate the configured object.
        :returns: A loaded Configurable object.
        """
        config['TYPE'] = self.TYPE
        # These options have no meaning anymore.
        config.pop('SUB_TYPE', None)
        config.pop('TAG', None)
        config.pop('NEXT', None)
        config.pop('TOP', None)
        file_hierarchy = config.pop('FILE_HIERARCHY', None)
        
        # Try and get the configurable class.
        try:
            cls = self.type_class.from_class_handle(config['class_name'])
        except ValueError:
            raise Configurable_loader_exception(config, self.TYPE, self.file_name, "class_name '{}' is not recognised".format(config['class_name'])) from None
        except KeyError:
            #raise Silico_exception("Error loading configurable of type '{}' from file '{}'; no class_name set".format(self.TYPE, config_path)) from None
            # If no class set, use the top level class.
            cls = self.type_class
        
        configurable = cls(file_hierarchy, False, **config) 
        configurable.configure_auto_name()
        if validate:
            configurable.validate()
        configurable.finalize()
        return configurable
    
    def size(self):
        """
        The recursive total number of child leaf nodes under this loader. For loaders that only represent one real configurable, size() == 1.
        """
        raise NotImplementedError()
    
    def link(self, loaders):
        """
        Link this partial loader to a number of other loaders.
        
        This implementation does nothing.
        """
        pass
    
    def find(self, tag):
        """
        Search through our child loaders (iteratively) for the first one with a given TAG.
        
        :param tag: The TAG to search for.
        """
        # A list of found configurables.
        found = self.search_for_tag(tag)
        
        # Panic if we've got nothing.
        if len(found) == 0:
            raise Silico_exception("Could not find a configurable with TAG '{}' that is a child of '{}'".format(tag, self.TAG))
        
        # We want the match that is closest to us (fewest path steps away).
        shortest = min((len(loader_path) for loader_path in found))
        
        # Prune all paths that are longer than our min.
        found = [loader_path for loader_path in found if len(loader_path) == shortest]
                
        # If we have more than one match, panic.
        if len(found) > 1:
            raise Unresolvable_tag_path_error(tag, found)
        
        else:
            return found[0]
        
    def search_for_tag(self, tag):
        """
        """
        # A list of found configurables.
        found = []
        
        if self.TAG == tag:
            # It's us!
            found.append([self])
            
        else:
            # None at this level, we need to ask our children if they match the tag.
            for child in self.NEXT:
                
                # Add the loaders our child found to our list, adding ourself to the start.
                for loader_list in child.search_for_tag(tag):
                    loader_list.insert(0, self)
                    found.append(loader_list)
            
        return found
    

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
    
    def resolve(self, identifier, validate = True):
        """
        Get one of the configurables that are represented by this loader.
        
        The configurable can be identified either with a unique index (1 - inf) or a unique list of tag names.
        
        :raises TypeError: If identifier is not an integer, list or tuple.
        :param identifier: The identifier to resolve.
        :param validate: Whether to call validate() on the final resolved configurable.
        :returns: A resolved configurable object.
        """
        # First, build or loader list.
        if isinstance(identifier, int):
            # Identifier is an index.
            path = self.path_by_index(identifier)
        
        elif isinstance(identifier, list) or isinstance(identifier, tuple):
            # Tag list.
            path = self.path_by_tags(identifier)
            
        else:
            # Unrecognised identifier
            raise TypeError("identifier should be either an integer or a list-like/tuple-like")
        
        # Now resolve our path.
        return self.resolve_path(path, validate = validate)
    
    def resolve_path(self, path, parent_config = None, validate = True):
        """
        Resolve a loader path, returning a single combined configurable object.
        
        :param path: A list of configurables to resolve.
        :param parent_config: The currently constructed dictionary of resolved options.
        :param validate: Whether to call validate() on the final resolved configurable.
        """
        if parent_config is None:
            parent_config = {}
        
        # First, merge our current parent object with ourself.
        self.merge_with_parent(parent_config)
        
        # Now continue down the loader path, removing the first item (which is us).
        return path[1].resolve_path(path[1:], parent_config = parent_config, validate = validate)
        
    def path_by_index(self, index, *, parent_offset = 0, path = None):
        """
        Build a list of loaders based on a unique index.
        
        Each configurable in the tree generated by this partial configurable can be identified by a unique index, corresponding to the position of that configurable if the tree was flattened in order.
        
        :raises IndexError: If index is out of range.
        :param index: The index to fetch.
        :param parent_offset: The index of this partial configurable in the parent configurable, used when called recursively.
        :param path: The currently built configurable path, used when called recursively.
        :returns: The loader path (a list).
        """
        path = [] if path is None else path
        
        # First, add ourself to the path.
        path.append(self)
                
        # We need to decide which child object from NEXT to continue to based on the range of possible indexes of each child.
        child_offset = 0
        for child in self.NEXT:
            child_size = child.size()
            if index < parent_offset + child_offset + child_size:
                # The config we want is in this child.
                return child.path_by_index(index, parent_offset = parent_offset + child_offset, path = path)
            
            # Add this child's size to the offset.
            child_offset += child_size
            
        # If we get this far, index is out of range.
        raise IndexError("Configurable index '{}' is out of range".format(index))
    
    def path_by_tags(self, tag_list, *, path = None):
        """
        Build a list of loaders based on a list of TAG names.
        
        :param tag_list: A list of ordered tags indicating which configurables to get.
        :param path: The currently built configurable path, used when called recursively.
        :returns: The list of loaders.
        """
        path = [] if path is None else path
            
        # First, check we've actually got a tag to search for.
        try:
            tag = tag_list[0]
        
        except IndexError:
            raise Tag_path_length_error(tag_list) from None

        # Search through our children for the next tag.
        # This function returns a list of loaders that lead to the tag we're looking for (like a path).
        # The first item is ourself, the last is the next child we'll call resolve_by_tags() on.
        next_children = self.find(tag)
        
        # Add our children to the list, except the last (which will add itself).
        path.extend(next_children[:-1])
                
        # Shorten our tag list.
        new_tag_list = tag_list[1:]
        
        # Now continue in the last child.
        try:
            return next_children[-1].path_by_tags(new_tag_list, path = path)
        
        except Tag_path_length_error:
            raise Tag_path_length_error(tag_list) from None
        

    def link(self, loaders):
        """
        Link this partial loader to a number of other loaders.
        
        This method will resolve the NEXT configs at this node.
        """
        try:
            for tag in self.config['NEXT']:
                # Look for a loader with this tag and add.
                # TODO: Should watch out for infinite recursion here.
                matching = [loader for loader in loaders if loader.TAG == tag]
                
                # Panic if we couldn't find any.
                if len(matching) == 0:
                    raise Configurable_loader_exception(self.config, self.TYPE, self.file_name, "next TAG '{}' could not be found".format(tag))
                
                # Add all matching to our NEXT.
                for match in matching:
                    # Unset the default TOP value.
                    match.TOP = False
                    self.NEXT.append(match)
                    
        except KeyError:
            if 'NEXT' not in self.config:
                # Missing required NEXT
                raise Configurable_loader_exception(self.config, self.TYPE, self.file_name, "missing required option NEXT") from None
                
            else:
                # Something else went wrong.
                raise


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
        
        # A flat version of all our child loaders.
        self.loaders = configs


class Single_loader(Configurable_loader):
    """
    A configurable loader that represents a single configurable, or the end of a chain of loaders that leads to a configurable.
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
    
    def resolve_path(self, path, parent_config = None, validate = True):
        """
        Resolve a loader path, returning a single combined configurable object.
        
        :param path: A list of configurables to resolve.
        :param parent_config: The currently constructed dictionary of resolved options.
        :param validate: Whether to call validate() on the final resolved configurable.
        """
        if parent_config is None:
            parent_config = {}
        
        # First, merge our current parent object with ourself.
        self.merge_with_parent(parent_config)
        
        return self.configure(parent_config, validate = validate)
    
    def path_by_index(self, index, *, parent_offset = 0, path = None):
        """
        Build a list of loaders based on a unique index.
        
        Each configurable in the tree generated by this partial configurable can be identified by a unique index, corresponding to the position of that configurable if the tree was flattened in order.
        
        :raises IndexError: If index is out of range.
        :param index: The index to fetch.
        :param parent_offset: The index of this partial configurable in the parent configurable, used when called recursively.
        :param path: The currently built configurable path, used when called recursively.
        :returns: The loader path (a list).
        """
        if path is None:
            path = []
        
        # Add ourself to the path.
        path.append(self)
        
        return path
    
    def path_by_tags(self, tag_list, *, path = None):
        """
        Build a list of loaders based on a list of TAG names.
        
        :param tag_list: A list of ordered tags indicating which configurables to get.
        :param path: The currently built configurable path, used when called recursively.
        :returns: The resolved Configurable object.
        """
        if path is None:
            path = []
        
        # Add ourself to the path.
        path.append(self)
        
        return path
    
class Update_loader(Single_loader):
    """
    A configurable loader that updates/overwrites another loader.
    """
    
    def resolve(self, *args, **kwargs):
        raise NotImplementedError("resolve() has no meaning for Update_loader objects")
    
    def resolve_path(self, *args, **kwargs):
        raise NotImplementedError("resolve_path() has no meaning for Update_loader objects")
    
    def path_by_index(self, *args, **kwargs):
        raise NotImplementedError("path_by_index() has no meaning for Update_loader objects")
    
    def path_by_tags(self, *args, **kwargs):
        raise NotImplementedError("path_by_tags() has no meaning for Update_loader objects")
    
    def update(self, loaders):
        """
        Update a number of other loaders with this loader.
        """
        # Find matching.
        try:
            matching = [loader for loader in loaders if loader.TAG == self.config['TAG']]
            
            if len(matching) == 0:
                # No matching, panic.
                raise Configurable_loader_exception(self.config, self.TYPE, self.file_name, "update TAG '{}' could not be found".format(self.config['TAG']))
            
            for match in matching:
                # Merge.
                deepmerge.always_merger.merge(match.config, self.config)
            
        except KeyError:
            if 'TAG' not in self.config:
                raise Configurable_loader_exception(self.config, self.TYPE, self.file_name, "missing required option TAG; don't know what to update") from None
            
            else:
                raise
        