import urwid

from silico.interface.urwid.method.node import Loader_top_node
from silico.interface.urwid.tree.base import Flaggable_tree_list_box,\
    Flaggable_tree_walker
from silico.config.configurable.option import Option
from silico.interface.urwid.swap.swappable import Swappable
from silico.interface.urwid.layout import Pane
from silico.interface.urwid.misc import Tab_pile
    

class Method_browser(Flaggable_tree_list_box):
    """
    Class for displaying and selecting computational methods.
    """


    def __init__(self, method_library, show_hidden):
        """
        Constructor for calculation browser objects.
        """
        self.method_library = method_library
        self.topnode = Loader_top_node(method_library.destinations, show_hidden = show_hidden)
        super().__init__(Flaggable_tree_walker(self.topnode))
        
    def is_selectable(self, node):
        """
        Determine whether a given node is selectable.
        """
        # Only base calculations are selectable.
        # We do this extra check incase there are destinations and/or programs without children.
        loader_list = node.get_value()        
        return not loader_list[-1].partial and node.loader_type == "calculation"
    
    def select(self, focus_node, focus_widget):
        """
        Select (or deselect) the node currently in focus.
        """
        super().select(focus_node, focus_widget)
#         
#         if focus_widget.flagged:
#             # This node is now selected, add to our other widget.
#             focus_widget.


class Method_selector(Swappable):
    """
    A tree list box widget used to browse and select methods.
    """
    
    show_hidden = Option(help = "Whether to show hidden methods.", default = False, type = bool)

    def __init__(self, top, methods):
        """
        Constructor for Method_selector objects.
        
        :param top: Top-most widget to use for display.
        :param methods: The methods that can be selected.
        """
        self.browser = Method_browser(methods, show_hidden = self.show_hidden)
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
        