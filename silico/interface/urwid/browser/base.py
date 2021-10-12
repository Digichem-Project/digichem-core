import urwid

from silico.interface.urwid.browser.node import Loader_top_node
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.browser.calcbox import Calcbox
from silico.interface.urwid.base import Section



class Method_browser(urwid.TreeListBox):
    """
    Class for browsing computational methods.
    """
    
    def __init__(self, top_node, selector):
        """
        Constructor for Method_browser objects.
        
        :param top_node: The top most node to display.
        :param selector: Method_selector object to add selected methods to.
        """
        super().__init__(urwid.TreeWalker(top_node))
        self.selector = selector

    def keypress(self, size, key):
        if key in ['enter', ' ']:
            # First, get the focus.
            focus_node = self.body.focus
            
            # If the node is a leaf, add it.
            # TODO: There might be a better way to ID leaf nodes...
            if not hasattr(focus_node, 'has_children'):
                # Get the method of the node in focus.
                method = focus_node.build_loader_path()
                self.selector.add_method(method)
            
            else:
                # Do nothing.
                return super().keypress(size, key)
        else:
            # Do nothing.
            return super().keypress(size, key)
    
    

class Method_selector(Tab_pile):
    """
    Class for displaying and selecting computational methods.
    """


    def __init__(self, methods):
        """
        Constructor for calculation browser objects.
        """
        self.topnode = Loader_top_node(methods)
        self.listbox = Method_browser(self.topnode, self)
        self.listbox.offset_rows = 1
        self.calcbox = Calcbox()
        
        # Construct our display widgets.
        super().__init__([
            ('weight', 2, Section(self.listbox, "Methods Browser")),
            ('weight', 1.5, Section(self.calcbox, title = "Selected Methods")),
            #(1, urwid.AttrMap(urwid.Filler(urwid.Padding(self.confirm, 'center', 11)), confirm_attr))
        ])
            
    def add_method(self, method):
        """
        Add the method that currently has focus to our calc box.
        
        :param method: The method to add.
        """        
        # Now add the method to our calcbox.
        self.calcbox.add_method(method)

