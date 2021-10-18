# Nodes (logical data storage) for the calculation browser.

# General imports.
import urwid
from silico.interface.urwid.browser.widget import Loader_leaf_widget,\
    Loader_parent_widget, Loader_top_widget


class Loader_node_mixin():
    """
    Mixin class for Treenodes that represent configurable loaders.
    """
        
    @property
    def loader_type(self):
        """
        Get the TYPE of the configurable loaders we represent.
        """
        return self.get_value()[-1].TYPE
    
    def resolve(self):
        """
        Resolve the configurable represented by this node.
        """
        # First, build our loader path.
        path = self.build_loader_path()[-1]
        
        # Now resolve (starting from the first loader in the path.
        return path[0].resolve_path(path)
        
    def build_loader_path(self, loader_paths = None):
        """
        Build a list of configurable loaders from this node up to the parent node.
        
        :returns: A list of loader paths (each of a different TYPE), each of which is a list of loaders.
        """
        if loader_paths is None:
            # A list of paths.
            loader_paths = [[]]
            
        # Add our list of loaders to the path.
        new_loader_path = list(self.get_value())
        new_loader_path.extend(loader_paths[0])
        loader_paths[0] = new_loader_path
        
        # Get our parent node.
        parent = self.get_parent()
        
        # If the parent is not None and is of the same TYPE, climb up the chain.
        if parent is not None and parent.loader_type == self.loader_type:
            return parent.build_loader_path(loader_paths)
        
        elif parent is not None and parent.loader_type != self.loader_type:
            # The parent is of a different type, start a new list.
            loader_paths.insert(0, [])
            return parent.build_loader_path(loader_paths)
        
        else:
            return loader_paths


class Loader_leaf_node(Loader_node_mixin, urwid.TreeNode):
    """
    Data storage object for leaf nodes
    """
    
    def load_widget(self):
        return Loader_leaf_widget(self)


class Loader_parent_node(Loader_node_mixin, urwid.ParentNode):
    """
    Class for Parent nodes that contain configurable loaders as children.
    """
    
    def __init__(self, loader_list, *args, **kwargs):
        """
        Constructor for Loader_parent_node objects.
        
        Each node can represent multiple configurable loaders, because pseudo loaders do not appear as their own node.
        Thus the series of loaders is stored as a list, of the form: [pseudo1, pseudo2... , partial]
        The name that the node will appear under is taken from the last loader in the list.
        
        :param loader_list: A list of the Configurable loaders we represent. The last loader in the list should have children.
        """
        urwid.ParentNode.__init__(self, loader_list, *args, **kwargs)
        self.concrete_children = []
    
    def load_widget(self):
        return Loader_parent_widget(self)

    def load_child_node(self, key):
        """
        Get a sub-node identified by a unique key.
        """
        child_loader_list = self.concrete_children[key]
        child_depth = self.get_depth() +1
        
        # Decide what kind of node we need.
        if len(child_loader_list[-1].sub_node_paths) != 0:
            # This node has children.
            child_class = Loader_parent_node
        
        else:
            # This loader has no children, it needs a leaf.
            child_class = Loader_leaf_node
        
        return child_class(child_loader_list, parent=self, key=key, depth=child_depth)
    
    def load_child_keys(self):
        """
        Get child keys, which uniquely identify each child.
        """
        loader_list = self.get_value()
        self.concrete_children = loader_list[-1].get_concrete_children()
        return range(len(self.concrete_children))


class Loader_top_node(Loader_parent_node):
    """
    The top-most node in the tree
    """
    
    def __init__(self, methods, *args, **kwargs):
        super().__init__([methods])
    
    def load_widget(self):
        return Loader_top_widget(self)