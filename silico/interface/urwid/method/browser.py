import urwid

from silico.interface.urwid.method.node import Loader_top_node
from silico.interface.urwid.tree.base import Flaggable_tree_list_box
from silico.config.configurable.option import Option
from silico.interface.urwid.top import View
    

class Method_browser(Flaggable_tree_list_box):
    """
    Class for displaying and selecting computational methods.
    """


    def __init__(self, methods, show_hidden):
        """
        Constructor for calculation browser objects.
        """
        self.topnode = Loader_top_node(methods, show_hidden = show_hidden)
        super().__init__(urwid.TreeWalker(self.topnode))
        
    def is_selectable(self, node):
        """
        Determine whether a given node is selectable.
        """
        # Only base calculations are selectable.
        # We do this extra check incase there are destinations and/or programs without children.
        loader_list = node.get_value()        
        return not loader_list[-1].partial and node.loader_type == "calculation"


class Method_selector(View):
    """
    A tree list box widget used to browse and select methods.
    """
    
    show_hidden = Option(help = "Whether to show hidden methods.", default = False, type = bool)

    def __init__(self, methods):
        """
        Constructor for Method_selector objects.
        
        :param methods: The methods that can be selected.
        """
        self.browser = Method_browser(methods, show_hidden = self.show_hidden)
        self.browser.offset_rows = 1
        super().__init__(self.browser, title = "Method Browser")
        
    def on_settings_change(self):
        """
        A method that will be called when settings have been changed.
        """
        # We need to update whether we're showing hidden stuff.
        self.browser.topnode._show_hidden = self.show_hidden
        self.browser.topnode.refresh()
        # TODO: This is in case a node was removed which was in focus, there is probably a better solution.
        self.browser.body.set_focus(self.browser.topnode)
        