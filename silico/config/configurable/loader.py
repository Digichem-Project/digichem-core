"""
Classes that represent one or more configurables.

Loaders construct real Configurable_class_target objects from one or more dict objects on demand.
Loaders are essentially a speed hack that means we don't have to construct 1000s of objects each time we start up.
"""

# General imports.
import deepmerge
import yaml
import copy

# Silico imports.
from silico.exception.configurable import Configurable_loader_exception,\
    Short_tag_path_error, Unresolvable_tag_path_error, Long_tag_path_error
from silico.exception.base import Silico_exception
from silico.config.configurable.base import Configurable_class_target
from silico.misc.base import is_iter
from silico.config.configurable.util import setopt

# Make methods available. We do this because our loaders are going to eventually ask for one of these classes.
# TODO: Importing all this here feels weird, perhaps this file should be moved to the submit package?
import silico.submit.destination
import silico.submit.program
import silico.submit.calculation


class Configurable_loader():
    """
    Abstract top-level class for configurable loaders.
    """
    
    # Whether this loader is a partial loader (a partial loader is any that is not a Single Loader).
    partial = True
    
    @property
    def TAG(self):
        """
        An identifying tag.
        """
        return self.config.get('TAG', None)
    
    @property
    def ALIAS(self):
        """
        An identifying alias.
        """
        return self.config.get('ALIAS', self.config.get('TAG', None))
    
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
        self.type_class = Configurable_class_target.from_class_handle(self.TYPE)
        
        # A list of the next loaders in the chain. For single loaders, len(self.NEXT) == 0.
        self.NEXT = []
        
        self._TOP = True
        self.pseudo = pseudo
        
        # Check our tag name is valid (does not contain /).
        if self.TAG is not None and  "/" in self.TAG:
            raise Configurable_loader_exception(self.config, self.TYPE, self.file_name, "the '/' character is not allowed in TAG names")
        
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
            
    @property
    def sub_node_paths(self):
        """
        A property which maps to the paths leading to the child nodes of this loader.
        For partial loaders this will be the NEXT attribute, while for single loaders (where NEXT is always empty) it will be the CHILDREN attribute.
        """
        return [[child] for child in self.NEXT]
    
    def get_concrete_children(self, show_hidden = False, _concrete_children = None, _current_path = None):
        """
        Get a list of the children of this loader that are concrete.
        
        Concrete here means that the child would appear as a node in a list of options, ie is not pseudo or similar.
        This method is recursive, any direct children that are not concrete (they are pseudo) will also have get_concrete_children() called.
        
        All parameters to this method are only used when called recursively; they should not normally be specified by the user.
        
        :param show_hidden: Whether to include hidden loaders.
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
            # We need to make a new list.
            current_path = list(_current_path)
        
        # Iterate through each of our immediate children.
        for child_path in self.sub_node_paths:
            # If any of the items in the child path are hidden, and we've been asked to ignore hidden items, do so.
            if not show_hidden:
                skip = False
                for child_path_item in child_path:
                    if child_path_item.config.get("meta", {}).get("hidden", False):
                        skip = True
                        
                if skip:
                    continue
            
            # Make a new path so we don't add to the same one for each child.
            new_path = list(current_path)
            
            # Add the child path to the new_path.    
            new_path.extend(child_path)
            
            if child_path[-1].pseudo:
                # This child is pseudo (ie, not concrete), so we want its children.
                
                # Continue in the child.
                child_path[-1].get_concrete_children(show_hidden = show_hidden, _concrete_children = concrete_children, _current_path = new_path)
                
            else:
                # Add the path to our total list.
                concrete_children.append(new_path)
                
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
        deepmerge.always_merger.merge(parent_config, copy.deepcopy(self.config))
                
        # Add ourself to the loader path.
        try:
            parent_config['loader_path'].append(self)
            
        except KeyError:
            parent_config['loader_path'] = [self]
    
    def configure(self, config, validate = True):
        """
        Convert (or attempt to) a config dict to an appropriate configurable object.
        
        :raises Silico_exception: If the class_name of the configurable is not set or cannot be found.
        :param config: The config dict.
        :param validate: Whether to validate the configured object.
        :returns: A loaded Configurable object.
        """
        config['meta']['TYPE'] = self.TYPE
        # These options have no meaning anymore.
        config.pop('SUB_TYPE', None)
        config.pop('TAG', None)
        config.pop('ALIAS', None)
        config.pop('NEXT', None)
        config.pop('PARENTS', None)
        config.pop('PREVIOUS', None)
        config.pop('TOP', None)
        loader_path = config.pop('loader_path', None)
        
        # Try and get the configurable class.
        try:
            cls = self.type_class.from_class_handle(config['meta']['class_name'])
        except ValueError:
            raise Configurable_loader_exception(config, self.TYPE, self.file_name, "meta:class_name '{}' is not recognised".format(config['meta']['class_name']))
        except KeyError:
            raise Configurable_loader_exception(config, self.TYPE, self.file_name, "no meta:class_name set") from None
            # If no class set, use the top level class.
            # IMPORTANT: It's not clear why this might be necessary so it has been disabled for now.
            # If this breaks something it will be reinstated.
            cls = self.type_class

        configurable = cls(loader_path, validate_now = validate, **config)
        configurable.finalize()
        return configurable
    
    def size(self):
        """
        The recursive total number of child leaf nodes under this loader. For loaders that only represent one real configurable, size() == 1.
        """
        raise NotImplementedError()
    
    def link(self, loaders, children):
        """
        Link this partial loader to a number of other loaders.
        
        This implementation does nothing.
        
        :param loaders: A flat list of the other loaders of the same type that have been parsed.
        :param children: The top loader (probably a configurable list) for the child type for this loader.
        """
        pass
    
    def find(self, tag):
        """
        Find the first child loader that has a given TAG.
        
        :param tag: The TAG to search for.
        :returns: A list of loaders leading to the one with the given tag.
        """
        # A list of found configurables.
        loader_lists = self.search_by_tag(tag)
        
        # Panic if we've got nothing.
        if len(loader_lists) == 0:
            if self.TAG is not None:
                message = "Could not find a definition with TAG '{}' that is a child of '{}'".format(tag, self.TAG)
            
            else:
                message = "Could not find a definition with TAG '{}'".format(tag)
            raise Silico_exception(message)
        
        # We want the match that is closest to us (fewest path steps away).
        shortest = min((len(loader_path) for loader_path in loader_lists))
        
        # Prune all paths that are longer than our min.
        loader_lists = [loader_path for loader_path in loader_lists if len(loader_path) == shortest]
                
        # If we have more than one match, panic.
        if len(loader_lists) > 1:
            raise Unresolvable_tag_path_error(tag, loader_lists)
        
        else:
            return loader_lists[0]
    
    def search_by_tag(self, tag):
        """
        Search through our child loaders (iteratively) for the first that matches a given tag.
        
        :param tag: The TAG to searcg for.
        :returns: A list of loader lists (a list leading to the loader with the given tag).
        """
        # A list of found configurables.
        loader_lists = []
        
        if self.TAG == tag:
            # It's us!
            loader_lists.append([self])
            
        else:
            # None at this level, we need to ask our children if they match the tag.
            for child in self.NEXT:
                # Add the loaders our child found to our list, adding ourself to the start.
                for loader_list in child.search_by_tag(tag):
                    loader_list.insert(0, self)
                    loader_lists.append(loader_list)
            
        return loader_lists
    

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
    
    def split_identifier_string(self, identifier, check_length = True):
        """
        Split a string identifying a number of method components (destination/program/calculation) into its constituent parts.
        
        :param identifier: The identifier to split, a string like.
        :param check_length: If True, raise an exception if an incorrect number of components are in identifier.
        :returns: The split parts (a list).
        """
        try:
            #tokens = list(ID_splitter(identifier))
            tokens = identifier.split("/")
            
        except Exception as e:
            raise ValueError("Could not split identifier string '{}'".format(identifier)) from e
        
        if check_length and len(tokens) != 3:
            raise ValueError("The identifier string '{}' contains {} components but must contain exactly 3 components".format(identifier, len(tokens)))
        
        new_tokens = []
        for token in tokens:
            new_tokens.append(yaml.safe_load(token))
#             if is_int(token):
#                 new_tokens.append(int(token))
#                 
#             else:
#                 new_tokens.append(token)
        
        return new_tokens
    
    def resolve_method(self, *identifiers, validate = True):
        """
        Resolve a set of identifiers into a tuple representing a method.
        
        :param identifiers: A number of identifiers to resolve.
        :param validate: Whether to call validate() on each of the final resolved configurables.
        """
        parts = []
        last = None
        next_top = self
        
        for identifier in identifiers:
            # First, resolve our identifier.
            configurable = next_top.resolve(identifier, validate = validate)
            path = configurable.loader_list
            
            # Check our new loader is a possible child of our last.
            if last is not None:
                if not last.valid_child_path(path):
                    raise Silico_exception("configurable '{}' is not a valid child of '{}'".format(
                        configurable.description,
                        parts[-1].description
                    ))

            
            next_top = path[-1].top_child
            last = path[-1]
            parts.append(configurable)
            
        return tuple(parts)
    
    def resolve_method_string(self, identifier, validate = True):
        """
        Resolve a string which identifies a complete method (consisting of a destination, a program and a calculation).
        
        Each part of the method is separated by a forward slash (/), and each part can identify the relevant configurable by index (number) or tag list.
        
        :param identifier: The identifier string.
        :param validate: Whether to call validate() on each of the final resolved configurables.
        """
        parts = self.split_identifier_string(identifier, check_length = True)
        return self.resolve_method(*parts, validate = validate)
    
    def resolve(self, identifier, validate = True):
        """
        Get one of the configurables that are represented by this loader.
        
        The configurable can be identified either with a unique index (1 - inf) or a unique list of tag names.
        
        :raises TypeError: If identifier is not an integer, str, list or tuple.
        :param identifier: The identifier to resolve.
        :param validate: Whether to call validate() on the final resolved configurable.
        :returns: A resolved configurable object and the path from which that object was resolved.
        """
        # First, build or loader list.
        if isinstance(identifier, int):
            # Identifier is an index.
            path = self.path_by_index(identifier)
        
        elif isinstance(identifier, list) or isinstance(identifier, tuple):
            # Tag list.
            path = self.path_by_tags(identifier, allow_incomplete = False)
        
        elif isinstance(identifier, str):
            # Single tag.
            path = self.path_by_tags([identifier], allow_incomplete = False)
        
        else:
            # Unrecognised identifier
            raise TypeError("identifier must be either an integer, str or a list-like/tuple-like, not '{}'".format(type(identifier)))
        
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
        try:
            return path[1].resolve_path(path[1:], parent_config = parent_config, validate = validate)
        
        except IndexError:
            # We ran out of parts of our path before reaching a single loader, give up.
            raise Short_tag_path_error([configurable.TAG for configurable in path])
        
    def index_of_path(self, path, *, parent_offset = 0):
        """
        Get the index of the configurable identified by a unique path.
        
        This method only makes sense if it is called from the top-most loader.
        
        :param path: A loader path, a list of configurable loaders.
        :param parent_offset: The current index total from previous iterations, used when called recursively. This should not normally be given by the user.
        :returns: The index.
        """
        # First, we need the index of the next path segment from our list of NEXT children.
        try:
            index = self.NEXT.index(path[1])
        
        except IndexError:
            # Ran out of path segments (or couldn't find the given segment?)
            raise Short_tag_path_error([configurable.TAG for configurable in path])
            
        # Next, we need the total of all the indexes we skipped (that are before our given index.
        skipped_total = sum([loader.size() for loader in self.NEXT[:index]])
        
        # Add this to our total parent offset from previous iterations.
        parent_offset += skipped_total
        
        # Continue in the next child.
        return path[1].index_of_path(path[1:], parent_offset = parent_offset)
           
        
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
            if (index -1) < parent_offset + child_offset + child_size:
                # The config we want is in this child.
                return child.path_by_index(index, parent_offset = parent_offset + child_offset, path = path)
            
            # Add this child's size to the offset.
            child_offset += child_size
            
        # If we get this far, index is out of range.
        raise IndexError("Configurable index '{}' is out of range".format(index))
    
    def path_by_tags(self, tag_list, *, path = None, allow_incomplete = True):
        """
        Build a list of loaders based on a list of TAG names.
        
        :param tag_list: A list of ordered tags indicating which configurables to get.
        :param path: The currently built configurable path, used when called recursively.
        :param allow_incomplete: Whether to return incomplete paths. Incomplete paths are those that do not end in a single loader. If False, a Short_tag_path_error exception will be raised if the given tag_list is not long enough. 
        :returns: The list of loaders.
        """
        path = [] if path is None else path
            
        # First, check we've actually got a tag to search for.
        try:
            tag = tag_list[0]
        
        except IndexError:
            # We ran out of path segments, if we're allowed to return incomplete paths do so.
            if allow_incomplete:
                # Add ourself as the final element.
                path.append(self)
                return path
            
            # We're not allowed incomplete paths, panic.
            raise Short_tag_path_error(tag_list) from None

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
        
        except Short_tag_path_error:
            raise Short_tag_path_error(tag_list) from None
        

    def link(self, loaders, children = None):
        """
        Link this loader to a number of other loaders.
        
        This method will resolve the NEXT configs at this node.
        
        :param loaders: A flat list of the other loaders of the same type that have been parsed.
        :param children: The top loader (probably a configurable list) for the child type for this loader.
        """
        try:
            for tag in self.config['NEXT']:
                # Look for a loader with this tag and add.
                # TODO: Should watch out for infinite recursion here.
                matching = [loader for loader in loaders if loader.TAG == tag]
                
                # Panic if we couldn't find any.
                if len(matching) == 0:
                    # Check we weren't accidentally given a list (because lists are allowed in NEXT for single loaders, but not for partial loaders).
                    if isinstance(tag, list):
                        raise Configurable_loader_exception(self.config, self.TYPE, self.file_name, "cannot use a NEXT tag '{}' that is a list here".format(tag))
                    
                    else:
                        raise Configurable_loader_exception(self.config, self.TYPE, self.file_name, "NEXT tag '{}' could not be found".format(tag))
                
                # Add all matching to our NEXT.
                for match in matching:
                    # Unset the default TOP value.
                    match.TOP = False
                    self.NEXT.append(match)
                    
        except (TypeError, KeyError):
            if 'NEXT' not in self.config:
                # Missing required NEXT
                raise Configurable_loader_exception(self.config, self.TYPE, self.file_name, "missing required option NEXT") from None
            
            elif not is_iter(self.config['NEXT']):
                raise Configurable_loader_exception(self.config, self.TYPE, self.file_name, "option NEXT is not iterable") from None
            
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
        super().__init__(file_name = None, TYPE = TYPE, config = {}, pseudo = True)
        
        self.NEXT = configs
        
        # A flat version of all our child loaders.
        self.loaders = configs


class Single_loader(Configurable_loader):
    """
    A configurable loader that represents a single configurable, or the end of a chain of loaders that leads to a configurable.
    """
    
    partial = False
    
    def __init__(self, file_name, TYPE, config):
        """
        Constructor for Single_loader objects.
        """
        super().__init__(file_name, TYPE, config)
        
        # A list of loaders that are children of this loader.
        # Child loaders will have a different TYPE to the current loader.
        # For example, loaders for destinations will have programs as child loaders, and programs will have calculations as child loaders.
        self.CHILDREN = []
        # The top loader (almost certainly a configurable list) of our child loaders.
        self.top_child = None
        
    def valid_child_path(self, possible_child_path):
        """
        Determine whether another configurable could be a valid child of this configurable.
        
        :param possible_child_path: A loader path (a list) that might be a child of this single loader.
        :returns: True or False.
        """
        for child_path in self.CHILDREN:
            if possible_child_path[:len(child_path)] == child_path:
                return True
        
        return False
        
    def size(self):
        """
        The size of this loader.
        
        Always returns 1.
        """
        return 1
    
    def get_children(self):
        """
        """
    
    @property
    def sub_node_paths(self):
        """
        A property which maps to the child nodes of this loader.
        For partial loaders this will be the NEXT attribute, while for single loaders (where NEXT is always empty) it will be the CHILDREN attribute,
        """
        return self.CHILDREN
    
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
    
    def index_of_path(self, path, *, parent_offset = 0):
        """
        Get the index of the configurable identified by a unique path.
        
        This method only makes sense if it is called from the top-most loader.
        
        :param path: A loader path, a list of configurable loaders.
        :param parent_offset: The current index total from previous iterations, used when called recursively. This should not normally be given by the user.
        :returns: The index.
        """
        return parent_offset + 1
    
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
        
        if len(tag_list) > 0:
            # Something's gone wrong, this is the last loader in the chain but we've got more segments to process?
            raise Long_tag_path_error(path, tag_list)
        
        return path
    
    def link(self, loaders, children = None):
        """
        Link this loader to a number of other loaders.
        
        :param loaders: A flat list of the other loaders of the same type that have been parsed.
        :param children: The top loader (probably a configurable list) for the child type for this loader.
        """
        # TODO: We need to check that our children are actually of a valid type for us, eg disallow Gaussian calcs to be children of Turbomole programs.
        if children is not None:
            # Save our top.
            self.top_child = children
            
            # Go through our list of next children.
            # Unlike for partial loaders, here NEXT refers to configurables of a different TYPE (they are our child configurables).
            for tag_list in self.config.get('NEXT', []):
                # If only a single tag has been given, wrap it in a list.
                if isinstance(tag_list, str):
                    tag_list = [tag_list]
                
                # Get a path leading to the child.
                child_path = children.path_by_tags(tag_list)
                
                # Add to our list.
                self.CHILDREN.append(child_path)
                
        elif len(self.config.get('NEXT', [])):
            # Panic, we have some loaders listed in NEXT but we have no children to search through.
            raise Configurable_loader_exception(self.config, self.TYPE, self.file_name, "cannot find NEXT loaders; there are no children for TYPE '{}'".format(self.TYPE))
    
    
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
        