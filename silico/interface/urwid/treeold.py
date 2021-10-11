import shlex
from logging import getLogger

import urwid
from urwidtrees.widgets import TreeBox
from urwidtrees.decoration import CollapsibleArrowTree
from urwidtrees.tree import SimpleTree

from silico.interface.urwid.misc import Tab_pile
import silico
from silico.interface.urwid.edit import Configurable_editor
from silico.interface.urwid.base import Top, Window
from silico.interface.urwid.dialogue import Confirm
from silico.exception.base import Silico_exception


class Tree_node(urwid.WidgetWrap):
    """
    Selectable text widget used for the individual nodes in the tree.
    """
    
    def __init__(self, name, *, ID, attr = "boldnode", focus_attr = "focus"):
        """
        Constructor for Tree_node objects.
        
        :param name: The name (the body) of this node.
        :param ID: The unique ID of this node (can be None).
        :param attr: Attribute to use for this node.
        :param node_attr: Alternative attribute that is used for nodes where ID == None.
        :param focus_attr: Attribute that is used when this node has focus.
        """
        t = urwid.Text(("{}) ".format(ID) if ID is not None else "") + name)    
        w = urwid.AttrMap(t, attr, focus_attr)
        super().__init__(w)
        self.ID = ID
        
        # This attribute contains a list of child nodes, but note this is not actually used by urwid trees itself, it is just a helper when we construct our list of nodes.
        self.name = name
        self._children = []
        
    def to_tuple(self):
        """
        Convert this Tree_node object (possibly with children) to a tuple of (node, children).
        
        :return: A 2 membered tuple of (node, children) where each child in children is a similar tuple, or children is None if there are no children.
        """
        #children = [self.from_node(child) for child in node.children]
        children = [child.to_tuple() for child in self._children] if len(self._children) > 0 else None
        #return (self(node.name, *args, ID = node.ID, **kwargs), children if len(children) > 0 else None)
        return (self, children)
    
    @classmethod
    def from_configurable_lists(self, configurable_lists, top_name = "Calculations"):
        """
        Create a Tree_node object (with children) from a list of Configurable_list objects.
        
        :param configurable_lists: A 3-membered list/tuple of Configurable_list objects (destinations, programs, calculations).
        :param top_name: The name of the top-most node.
        :return: A 2 membered tuple of (node, children).
        """
        # First, get a top-level node.
        base_node = self(top_name, ID = None)
        
        # Add.
        for config in configurable_lists[0]:
            self.from_configurable(base_node, config, configurable_lists)
            
        # Return
        return base_node.to_tuple()
    
    @classmethod
    def from_configurable(self, base_node, configurable, configurables, **kwargs):
        """
        Create a Tree_node object from a Configurable object.
        
        You are probably looking for from_configurable_lists().
        
        :param base_node: The node under which this new node will be added.
        :param configurable: The configurable to construct from.
        :param configurables: A list/tuple of Configurable_list objects.
        :param **kwargs: Passed to the new node's constructor.
        """
        # First, turn our config into a Node.
        for group_name in configurable.GROUP:
            if len(base_node._children) == 0 or base_node._children[-1].name != group_name:
                # Either the base node has no children, or the last child is not our group.
                # Add a new group.
                base_node._children.append(self(group_name, ID = None))
                
            # Set this new node as our base.
            base_node = base_node._children[-1]
            
        # Now add the config to the (current) base node.
        node = Configurable_node(configurable, **kwargs)
        base_node._children.append(node)
        
        # Add our children, using the new node as base.
        if len(configurables) > 1:
            for child in configurable.get_children(configurables[1]):
                self.from_configurable(node, child, configurables[1:])
            
        # Done, return Node for convenience (although it probably isn't very useful).
        return node
    
    def selectable(self):
        return True

    def keypress(self, size, key):
        return key
    
class Category_node(Tree_node):
    """
    Variation of a tree node used to represent a category.
    """
    
class Configurable_node(Tree_node):
    """
    Variation of a Tree_node that accepts a Configurable object as argument.
    """
    
    def __init__(self, configurable, *, attr = "body", focus_attr = "focus"):
        """
        Constructor for Configurable_node objects.
        
        :param configurable: The configurable object to construct from.
        :param attr: Attribute to use for this node.
        :param node_attr: Alternative attribute that is used for nodes where ID == None.
        :param focus_attr: Attribute that is used when this node has focus.
        """
        name = configurable.NAME if configurable.GROUP_NAME is None else configurable.GROUP_NAME
        self.configurable = configurable
        
        # Decide on which attribute to use.
        #if self.configurable.hidden:
        #    attr = "hiddennode"
        if self.configurable.warning is not None:
            attr = "warningnode"
        
        # Add exclamation if we have a warning.
        #if self.configurable.warning is not None:
        #    name += " (!)"
        #if self.configurable.hidden:
        #    name += " (hidden)"
        
        try:            
            name += " ({})".format(configurable.status)
        except AttributeError:
            # This is ok.
            pass
        except Exception:
            getLogger(silico.logger_name).debug("Failed to retrieve status for '{}'".format(configurable.NAME), exc_info = True)
        
        super().__init__(name, ID = configurable.ID, attr = attr, focus_attr = focus_attr)
    
class Calculation_tree(TreeBox):
    """
    Modified TreeBox.
    """
    
    def __init__(self, destinations, programs, calculations, *, show_hidden = False, include_top = False, calcbox, top, **kwargs):
        """
        :param calcbox: Editable widget that will be populated with selected options.
        :param top: The outer containing Top widget.
        """
        # Save our configurables.
        self.destinations = destinations
        self.programs = programs
        self.calculations = calculations
        
        # Get a converted node.
        self.show_hidden = show_hidden
        self.include_top = include_top
        
        # Get our top node.
        if self.show_hidden:
            top_node = Tree_node.from_configurable_lists((self.destinations, self.programs, self.calculations))
        else:
            top_node = Tree_node.from_configurable_lists((self.destinations.visible(), self.programs.visible(), self.calculations.visible()))
        
        # Now update our forest accordingly.
        forest = (top_node,) if self.include_top else top_node[1]
        
        if forest is None:
            # We have nothing to display, which breaks urwid trees quite badly. Give up.
            raise Silico_exception("Unable to render browser; there is nothing to display.")
        
        super().__init__(CollapsibleArrowTree(
                SimpleTree(forest),
                is_collapsed = lambda x: True
            ), **kwargs)
        
        self.calcbox = calcbox
        self.top = top
    
    def keypress(self, size, key):
        """
        Handler for keypress events.
        """
        if key in ['enter', ' ']:
            self.select()
        elif key in ['m', 'e']:
            self.edit()
        elif key in ['backspace']:
            super().keypress(size, '-')
        elif key in ['delete']:
            # Delete the last calculation.
            self.calcbox.calculations = self.calcbox.calculations[:-1]
        else:
            # Not for us to deal with.
            return super().keypress(size, key)
        
    def select(self):
        """
        Handler function when an item is selected.
        """
        # Toggle expanded state.
        # First, get the focused widget.
        focuspos = self.get_focus()[1]
        
        # Check we can expand (and selected is not a leaf with no children).
        if not self._tree._tree.is_leaf(focuspos):
            # Determine whether widget is currently expanded.
            if self._tree.is_collapsed(focuspos):
                #super().keypress(size, '+')
                self.expand_focussed()
            else:
                #super().keypress(size, '-')
                self.collapse_focussed()
        else:
            # Leaf node, add selected.
            self.calcbox.edit_text += (" " if self.calcbox.edit_text != "" else "") + self.get_calc_string(focuspos) + " "
    
    def expand_focussed(self, *args, _no_prompt = False, **kwargs):
        """
        Overridden method to prompt for confirmation in certain situations.
        """
        focus_widget = self._tree[self.get_focus()[1]].base_widget
        
        # If we have a configurable node with a warning, prompt.
        if not _no_prompt and hasattr(focus_widget, 'configurable') and getattr(focus_widget.configurable, 'warning', None) is not None:
            self.top.swap(Confirm('WARNING', focus_widget.configurable.warning, lambda: self.expand_focussed(_no_prompt = True), self.top))
        else:
            super().expand_focussed(*args, **kwargs)
                
        
    def edit(self):
        """
        Method called when the user uses shift-enter on a calculation.
        """
        #focus_widget = self.get_focus()[0].base_widget
        focus_widget = self._tree[self.get_focus()[1]].base_widget
        
        #print(type(focus_widget).__name__)
        
        # Only continue if we have a Configurable_node.
        if hasattr(focus_widget, 'configurable'):
            # Get an editor widget.
            editor = Configurable_editor(focus_widget.configurable, top = self.top)
            self.top.swap(editor.window())
            
    def get_calc_string(self, pos):
        """
        Get a m/p/c calculation string from a selected child widget.
        """
        parent = self._tree.parent_position(pos)
        if parent is not None:
            parentstr = self.get_calc_string(parent)
        else:
            parentstr = None
            
        selfID = self._tree[pos].base_widget.ID
        if selfID is not None:
            if parentstr is not None:
                return parentstr + "/" + str(selfID)
            else:
                return str(selfID)
        else:
            if parentstr is not None:
                return parentstr
            else:
                return None
            
class Calcbox(urwid.Edit):
    """
    Modified edit widget that offers methods for handling calculation strings.
    """
    
    @property
    def calculations(self):
        """
        Get a list of the calculation strings currently entered into this calcbox.
        
        Note that modifications to the returned list are not automatically reflected in this edit widget.
        """
        return shlex.split(self.edit_text)
    
    @calculations.setter
    def calculations(self, value):
        """
        Set the text of this edit widget from a list of calculation strings.
        """
        self.edit_text = "  ".join(value)
        
        

class Calculation_browser(Tab_pile):
    """
    Urwid box widget that allows selecting calculations.
    """
    
    def __init__(self, destinations, programs, calculations, *, show_hidden = False, confirm_text = "Confirm", calcbox_attr = "body", calcbox_inner_attr = "editable", confirm_attr = "good_button", top, **kwargs):
        """
        Constructor for Calculation_browser widgets.
        """
        self.destinations = destinations
        self.programs = programs
        self.calculations = calculations
        self.show_hidden = show_hidden
        self.top = top
        self.browser_args = kwargs
        
        
        # Whether the tree exited normally (ie, by using the confirm button).
        self.confirmed = False
        
        # Each browser is made up of 3 main child widgets; the browser itself, a text field (calcbox) and the confirm button.
        self.calcbox = Calcbox((calcbox_attr, 'Selected calculations: '))
        self.confirm = urwid.Button(confirm_text, self.stop)
        
        # Call our parent constructor to make the pile.
        super().__init__([
            #urwid.LineBox(self.browser),
            urwid.LineBox(urwid.Filler(urwid.Edit())),
            (4, urwid.LineBox(urwid.AttrMap(urwid.Filler(self.calcbox), calcbox_inner_attr))),
            (1, urwid.AttrMap(urwid.Filler(urwid.Padding(self.confirm, 'center', 11)), confirm_attr))
        ])
        
        self.update_browser()
    
    def keypress(self, size, key):
        """
        Handle keypress events.
        """
        if key in ['s']:
            self.show_hidden = not self.show_hidden
            self.update_browser()
        else:
            return super().keypress(size, key)
    
    def update_browser(self):
        """
        Update the browser window.
        
        Sadly, urwidtrees does not respond well to in-place changes. As such, the entire treebox is regenerated (and replaced) by this method.
        """
        # First, get the configurables we're going to display. This depends on whether we are showing hidden items or not.
        if self.show_hidden:
            destinations = self.destinations
            programs = self.programs
            calculations = self.calculations
        else:
            destinations = self.destinations.visible()
            programs = self.programs.visible()
            calculations = self.calculations.visible()
            
        self.browser = Calculation_tree(
            destinations,
            programs,
            calculations,
            show_hidden = self.show_hidden,
            calcbox = self.calcbox,
            top = self.top,
            **self.browser_args
        )
        
        # Now update our pile.
        self.contents[0] = (urwid.LineBox(self.browser), ('weight', 1))
    
    def stop(self, *args, **kwargs):
        """
        This method is called when the user clicks the confirm button.
        """
        self.confirmed = True
        raise urwid.ExitMainLoop() 
        
    @classmethod
    def run(self, destinations, programs, calculations, palette, include_top = False):
        """
        Interactively run urwid, using a Calculation_browser as the top-most widget.
        
        :param node_tuple: The top-most node used to populate the choices in the calculation browser as a 2 membered tuple of (node, children).
        :param include_top: Whether to include the top-most node in the browser. If False, the children will instead be used as top nodes.
        :return: The selected calculations, as a string.
        """
        
        top = Top(urwid.SolidFill())
        browser = Calculation_browser(destinations, programs, calculations, include_top = include_top, top = top)
        top.original_widget = Window(
            urwid.AttrMap(browser, 'body'),
            title = 'Silico {} Calculation Browser'.format(silico.version),
            help = 'ENTER: select   m: modify   s: show    DELETE: delete   E: expand all   C: contract all   Esc: quit'
        )
    
        #enclose in a frame
        urwid.MainLoop(
            top,
            palette
        ).run()  # go
    
        return browser.calcbox.edit_text if browser.confirmed else ""
    
   
    
