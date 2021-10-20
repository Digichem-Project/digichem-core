import urwid

from silico.interface.urwid.browser.node import Loader_top_node
from silico.interface.urwid.tree.base import Flaggable_tree_list_box
from silico.config.configurable.option import Option
from silico.interface.urwid.top import View
    

class Method_browser(Flaggable_tree_list_box):
    """
    Class for displaying and selecting computational methods.
    """


    def __init__(self, methods):
        """
        Constructor for calculation browser objects.
        """
        topnode = Loader_top_node(methods)
        super().__init__(urwid.TreeWalker(topnode))    


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
        self.browser = Method_browser(methods)
        self.browser.offset_rows = 1
        
        super().__init__(self.browser, "Method Browser")
        