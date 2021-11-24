import urwid

from silico.interface.urwid.method.node import Loader_top_node
from silico.interface.urwid.tree.base import Flaggable_tree_list_box,\
    Flaggable_tree_walker
from silico.config.configurable.option import Option
from silico.interface.urwid.tree.selector import Enhanced_tree_selector
    

class Method_browser(Flaggable_tree_list_box):
    """
    Class for displaying and selecting computational methods.
    """


    def __init__(self, methods, show_hidden):
        """
        Constructor for calculation browser objects.
        """
        self.methods = methods
        self.topnode = Loader_top_node(methods, show_hidden = show_hidden)
        super().__init__(Flaggable_tree_walker(self.topnode))
        
    def is_selectable(self, node):
        """
        Determine whether a given node is selectable.
        """
        # Only base calculations are selectable.
        # We do this extra check incase there are destinations and/or programs without children.
        loader_list = node.get_value()        
        return not loader_list[-1].partial and node.loader_type == "calculation"


class Method_selector(Enhanced_tree_selector):
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
        browser = Method_browser(methods, show_hidden = self.show_hidden)
        browser.offset_rows = 1
        
        manual_widget = urwid.Edit(("body", "Codes: "))
        
        super().__init__(top, browser, manual_widget, title = "Method Browser", manual_widget_title = "Manual Method Codes")
        
    def on_settings_change(self):
        """
        A method that will be called when settings have been changed.
        """
        # We need to update whether we're showing hidden stuff.
        self.browser.topnode._show_hidden = self.show_hidden
        self.browser.topnode.refresh()
        