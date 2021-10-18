import urwid

from silico.interface.urwid.browser.node import Loader_top_node
from silico.interface.urwid.base import Section
from silico.interface.urwid.tree.base import Flaggable_tree_list_box
    

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


class Method_selector(Section):
    """
    A tree list box widget used to browse and select methods.
    """

    def __init__(self, methods):
        """
        Constructor for File_selector objects.
        
        :param starting_dir: The starting directory that will be shown expanded.
        """
        self.browser = Method_browser(methods)
        self.browser.offset_rows = 1
        
        super().__init__(self.browser, "Method Browser")