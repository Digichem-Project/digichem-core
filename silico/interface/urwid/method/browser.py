import urwid

from silico.interface.urwid.method.node import Loader_top_node
from silico.interface.urwid.tree.base import Flaggable_tree_list_box,\
    Flaggable_tree_walker
from silico.config.configurable.option import Option
from silico.interface.urwid.swap.swappable import Swappable
from silico.interface.urwid.layout import Pane
    

class Method_browser(Flaggable_tree_list_box):
    """
    Class for displaying and selecting computational methods.
    """


    def __init__(self, methods, show_hidden, one_type_only = False):
        """
        Constructor for calculation browser objects.
        """
        self.methods = methods
        self.topnode = Loader_top_node(methods, show_hidden = show_hidden, stop_on_single = one_type_only)
        super().__init__(Flaggable_tree_walker(self.topnode))
        
    def is_selectable(self, node):
        """
        Determine whether a given node is selectable.
        """
        # Normally, only base calculations are selectable.
        # We do this extra check incase there are destinations and/or programs without children.
        # However, if we're only showing one type of method, then we'll allow any node without children to be selected.
        loader_list = node.get_value()
        if self.topnode.stop_on_single:
            return not loader_list[-1].partial
        else:
            return not loader_list[-1].partial and node.loader_type == "calculation"
    
    def select(self, focus_node, focus_widget):
        """
        Select (or deselect) the node currently in focus.
        """
        super().select(focus_node, focus_widget)


class Method_selector(Swappable):
    """
    A tree list box widget used to browse and select methods.
    """
    
    show_hidden = Option(help = "Whether to show hidden methods.", default = False, type = bool)

    def __init__(self, top, methods, one_type_only = False):
        """
        Constructor for Method_selector objects.
        
        :param top: Top-most widget to use for display.
        :param method_library: An object that contains methods, destinations, programs and calculations as attributes.
        """
        self.browser = Method_browser(methods, show_hidden = self.show_hidden, one_type_only = one_type_only)
        self.browser.offset_rows = 1
        
        self.code_widget = urwid.AttrMap(urwid.Edit(multiline = True), "editable")
        
        # Assemble our widgets.
        #body = Tab_pile([Pane(self.browser, title = "Method Browser"), ('weight', 0.3, Pane(urwid.Filler(self.code_widget, valign = "top"), "Method Codes"))])
        body = Pane(self.browser, title = "Method Browser")
        
        super().__init__(top, body)
        
    @property
    def selected(self):
        return self.browser.selected
    
    def reset(self):
        return self.browser.reset()
        
    def on_settings_change(self):
        """
        A method that will be called when settings have been changed.
        """
        # We need to update whether we're showing hidden stuff.
        self.browser.topnode._show_hidden = self.show_hidden
        self.browser.topnode.refresh()
        