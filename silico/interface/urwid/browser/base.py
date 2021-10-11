import urwid

from silico.interface.urwid.browser.node import Loader_top_node
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.browser.calcbox import Calcbox, Calcbox_item



class Calculation_browser():
    """
    Class for displaying an urwid compatible calculation browser.
    """


    def __init__(self, methods, palette):
        """
        Constructor for calculation browser objects.
        """
        self.palette = palette
        self.topnode = Loader_top_node(methods)
        self.listbox = urwid.TreeListBox(urwid.TreeWalker(self.topnode))
        self.listbox.offset_rows = 1
        self.calcbox = Calcbox()
        
        # Construct our display widgets.
        self.view = Tab_pile([
            ('weight', 2, urwid.LineBox(self.listbox, title = "Methods Browser", title_align = "left", title_attr = "bold")),
            ('weight', 1, urwid.LineBox(self.calcbox, title = "Selected Methods", title_align = "left", title_attr = "bold")),
            #(1, urwid.AttrMap(urwid.Filler(urwid.Padding(self.confirm, 'center', 11)), confirm_attr))
        ])

    def main(self):
        """Run the program."""

        self.loop = urwid.MainLoop(urwid.AttrMap(self.view, 'body'), self.palette,
            unhandled_input=self.unhandled_input)
        self.loop.run()

    def unhandled_input(self, key):
        if key in ('q','Q'):
            raise urwid.ExitMainLoop()
        
        if key in ['enter', ' ']:
            # First, get the focus.
            focus_node = self.listbox.body.focus
            
            # If the node is a leaf, add it.
            # TODO: There might be a better way to ID leaf nodes...
            if not hasattr(focus_node, 'has_children'):
                self.add_method(focus_node)
            
            else:
                # Expand this parent node instead.
                pass
            
    def add_method(self, focus):
        """
        Add the method that currently has focus to our calc box.
        """
        # Get the method based on the focus node.
        method = focus.build_loader_path()
        
        # Now add the method to our calcbox.
        calcbox_item = Calcbox_item(method, self.calcbox)
        self.calcbox.body.append(calcbox_item)
        self.calcbox.update_positions()

