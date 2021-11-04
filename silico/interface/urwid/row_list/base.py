import urwid
from silico.interface.urwid.misc import Blank
from silico.interface.urwid.dialogue import Confirm_dialogue

class Row_control():
    """
    A button to control a Row_item object.
    """

class Row_widget(urwid.WidgetPlaceholder):
    """
    A widget used for displaying a Row_item object.
    """
    
    attr_name = "body"
    focus_attr_name = "editable"
    
    def __init__(self, row_item):
        """
        Constructor for Row_widget objects.
        
        :param item: The Row_item we represent.
        """
        self._inner = None
        self._widget = None
        self._controls = None
        self.row_item = row_item
        
        # A widget we'll use to keep track of our position.
        self.position_widget = urwid.Text("***")
        
        # Widgets to give us spacing between ourself and the next row item.
        #self.divider = urwid.Divider()
        self.divider = Blank()
        self.divider_placeholder = urwid.WidgetPlaceholder(self.divider)
        
        super().__init__(urwid.AttrMap(self.get_widget(), self.attr_name, self.focus_attr_name))
        
    @property
    def movable(self):
        """
        Whether we are allowed to rearrange our position.
        """
        return self.row_item.movable
        
    def get_widget(self, reload = False):
        """
        Get the inner widget we'll use for display.
        """
        if self._widget is None or reload:
            self._widget = self.load_widget()
            
        return self._widget
    
    def load_widget(self):
        """
        Load our inner widget that we'll use for display.
        """
        row_items = []
        
        # If we are movable, add a position widget.
        if self.movable:
            row_items.append(('pack', urwid.AttrMap(self.position_widget, "bold")))
            
        # Add our body.
        row_items.append(self.get_inner())
        
        # Add our controls.
        controls = self.get_controls()
        if len(controls) > 0:
            control_box = urwid.GridFlow(controls, 5, 0, 0, "right")
            row_items.append((10, control_box))
        
        return urwid.Pile([
            urwid.Columns(row_items, dividechars = 2),
            self.divider_placeholder
        ])
    
    def get_inner(self, reload = False):
        """
        Get the widget we'll use to display our main body.
        """
        if self._inner is None or reload:
            self._inner = self.load_inner()
            
        return self._inner
    
    def load_inner(self):
        """
        Load the widget we'll use to display our main body.
        """
        return urwid.Text(self.row_item._value)
        
    
    def get_controls(self, reload = False):
        """
        Get the buttons we use for control.
        """
        if self._controls is None or reload:
            self._controls = self.load_controls()
            
        return self._controls
        
    def load_controls(self):
        """
        Load the buttons that can be used to control this row.
        
        This implementation returns up and down controls (if we are movable) and a remove control.
        
        :returns: A list of control buttons.
        """
        up_button = urwid.Button("↑", lambda button: self.row_item.move_up())
        down_button = urwid.Button("↓", lambda button: self.row_item.move_down())
        delete_button = urwid.Button("r", lambda button: self.row_item.remove())
        
        controls = []
        if self.movable:
            controls.append(urwid.AttrMap(up_button, "button--small", "button--small--focus"))
            controls.append(urwid.AttrMap(down_button, "button--small", "button--small--focus"))
            
        controls.append(urwid.AttrMap(delete_button, "button--bad--small", "button--small--focus"))
        return controls
        
    def update_position(self):
        """
        Update our position widget with our current position in the parent calcbox.
        """
        # Get our index.
        try:
            position = self.row_item.index +1
            position_str = str(position)
            
        except IndexError:
            position_str = "***"
        
        self.position_widget.set_text(position_str)
        
        # If we're in the last position, hide our divider.
        if position == len(self.row_item.row_list.body):
            self.divider_placeholder.original_widget = Blank()
            
        else:
            self.divider_placeholder.original_widget = self.divider


class Row_pointer_widget(Row_widget):
    """
    A widget used for displaying the row pointer.
    """
    
    def load_inner(self):
        """
        Load the widget we'll use to display our main body.
        """
        return urwid.Button(self.row_item._value, lambda button: self.row_item.add_func())
    
    def load_controls(self):
        """
        Load the buttons that can be used to control this row.
        
        This implementation returns up and down controls (if we are movable) and a remove control.
        
        :returns: A list of control buttons.
        """    
        controls = super().load_controls()
        # remove the delete button.
        controls = controls[:-1]
        return controls


class Row_item():
    """
    A single row in a Row_list.
    """
    
    def __init__(self, value, row_list, movable = True):
        """
        
        :param value: The contents (data) of this row.
        :param row_list: The Row_list object we belong to.
        """
        self.movable = movable
        self._value = value
        self.row_list = row_list
        self._widget = None
        
    def get_value(self):
        """
        Get the updated value of this row.
        """
        return self._value
        
    def load_widget(self):
        """
        Load a widget that represents this Row_item().
        
        :returns: The loaded widget.
        """
        return Row_widget(self)
        
    def get_widget(self, reload = False):
        """
        Get the widget that represents us.
        """
        if self._widget is None or reload:
            self._widget = self.load_widget()
            
        return self._widget
        
    @property
    def index(self):
        """
        Return the current index of this row.
        """
        return self.row_list.body.index(self)
    
    def move_up(self):
        """
        Move this row up one.
        """
        return self.move('up', 1)
    
    def move_down(self):
        """
        Move this row down one.
        """
        return self.move('down', 1)
        
    def move(self, position, amount):
        """
        Move this row item either up or down.
        
        :param position: One of either 'up' or 'down'.
        :param amount: The amount to move up or down.
        """
        # Get our current index.
        cur_index = self.index
        
        if position == "up":
            new_index = cur_index -amount
            
        elif position == "down":
            new_index = cur_index +amount
            
        else:
            raise ValueError("position must be one of either 'up' or 'down'")
        
        # Don't allow negatives
        if new_index < 0:
            new_index = 0
        
        # We need to insert and pop to move.
        self.row_list.body.insert(new_index, self.row_list.body.pop(cur_index))
                
        # Set ourself as the focus.
        self.row_list.body.set_focus(new_index if new_index < len(self.row_list.body) else len(self.row_list.body) -1)
        
    def remove(self):
        """
        Remove this Row_item from its parent.
        """
        index = self.index
        self.row_list.body.pop(index)


class Row_pointer(Row_item):
    """
    A dummy row item that indicates where the next item will be inserted.
    """
    
    def __init__(self,  row_list, add_func, movable = True, widget_text = "Add new here"):
        """
        Constructor for Row_pointer objects.
        
        :param value: 
        """
        Row_item.__init__(self, widget_text, row_list, movable)
        # A function we'll call when our button is pressed.
        self.add_func = add_func
        
    def load_widget(self):
        """
        Load a widget that represents this Row_item().
        
        :returns: The loaded widget.
        """
        return Row_pointer_widget(self)


class Row_walker(urwid.SimpleFocusListWalker):
    """
    ListWalker-compatible class for displaying Row_item objects.

    Positions are Row_item objects, while actual display is handled by Row_widgets.
    """

    def __init__(self, data):
        """
        Constructor for Row_walker objects.
        
        :param data: List of data (Row_items) to display.
        """
        super().__init__(data)
        
    
    def _modified(self):
        """
        Function called whenever our list is changed.
        """
        super()._modified()
        self.update_positions()
        
    def update_positions(self):
        """
        Tell each child to look up their positions again.
        """
        for row in self:
            row.get_widget().update_position()
        
    @property
    def focus_row(self):
        """
        Return the Row_item currently in focus.
        """
        return self[self.focus]
    
    @focus_row.setter
    def focus_row(self, value):
        """
        Set the focus to a Row_item.
        """
        self.set_focus(self.index(value))
    
    def get_focus(self):
        """
        Get the current focus.
        
        :returns: A tuple (widget, position) of the current focus.
        """
        try:
            focus = self.focus
            return self[focus].get_widget(), focus
        except (IndexError, KeyError, TypeError):
            return None, None
        
    def get_next(self, position):
        """
        Get the next.
        
        :returns: A tuple (widget, position) of position.
        """
        try:
            position = self.next_position(position)
            return self[position].get_widget(), position
        except (IndexError, KeyError):
            return None, None

    def get_prev(self, position):
        """
        Get the previous.
        
        :returns: A tuple (widget, position) of position.
        """
        try:
            position = self.prev_position(position)
            return self[position].get_widget(), position
        except (IndexError, KeyError):
            return None, None


class Row_list(urwid.ListBox):
    """
    A custom ListBox that displays data in (movable) rows.
    """
    
    def __init__(self, add_func, rearrangeable = True, initial = None):
        """
        
        :param add_func: A function to call when a new item is to be added.
        :param rearrangeable: Whether the rows of this Row_list can be rearranged.
        :param initial: An optional list of initial values to populate the list with.
        """
        initial = [] if initial is None else initial
        
        self.rearrangeable = rearrangeable
        
        # Convert our initial values into actual rows.
        initial = [self.load_row(value) for value in initial]
        
        # Keep track of our pointer
        self.pointer = Row_pointer(self, add_func, rearrangeable)
        # Add to our list of starting items.
        initial.append(self.pointer)
        
        # Construct.
        super().__init__(Row_walker(initial))
        
        # Set our initial focus.
        self.set_focus(0)
        
        self.body.update_positions()
            
    def load_row(self, value):
        """
        Return an object that represents one of our rows.
        
        :param value: The value of the new row.
        :returns: The Row_item object.
        """
        return Row_item(value, self, self.rearrangeable)
            
    def add_row(self, value):
        """
        Add a new row to this Row_list.
        
        :param value: The value of the new row.
        """
        row_item = self.load_row(value)
        
        # Find the current position of our pointer.
        index = self.body.index(self.pointer)
        
        self.body.insert(index, row_item)
        
        # Set the pointer as the current focus.
        self.set_focus(index +1)
        
        
class Row_browser(Row_list):
    """
    A higher level Row_list object that can show a browser to add nodes.
    """
    
    def __init__(self, selector, top, rearrangeable=True, initial = None):
        """
        Constructor for Method_list objects.
        
        :param selector: The selector widget to show when the user wants to add stuff.
        :param top: The topmost widget to use for display.
        """
        self.top = top
        
        self.selector = selector
        
        super().__init__(self.swap_to_browser, rearrangeable=rearrangeable, initial=initial)
        
    def swap_to_browser(self):
        """
        Function called to open the calculation browser.
        """
        self.top.swap_into_window(self.selector, cancel_callback = self.cancel, submit_callback = self.submit)
        
    def submit(self):
        """
        Add selected methods to our list.
        """
        errors = []
        # Go through the nodes the user selected.
        for selected_node in self.selector.browser.selected_nodes:
            value = None
            try:
                value = self.value_from_node(selected_node)
                self.add_row(value)
            
            except Exception as error:
                # We couldn't load the node for some reason, show an error and ignore.
                errors.append((value if value is not None else selected_node.get_value(), error))
                
        # Clear the selected nodes.
        self.selector.browser.reset()
        
        # If we got any errors, show them to the user.
        if len(errors) > 0:
            error_title, error_text = self.get_error_text(errors)
            
            self.top.swap(Confirm_dialogue(error_title, self.top, error_text, error = True, submit_callback = 2))
            # Returning False stops us swapping back before we show the dialogue.
            return False
        
    def get_values(self):
        """
        Get a list of values currently represented by this row list.
        """
        return [item.get_value() for item in self.body if item != self.pointer]
        
    def value_from_node(self, node):
        """
        Get a value to add from a selected node.
        """
        return node.get_value()
        
    def get_error_text(self, errors):
        """
        Get the text to display when an error occurs adding to our List.
        
        :param: A tuple of (title, message)
        """
        return ("Error", "The following {} value(s) could not be loaded and will be ignored:\n\n".format(len(errors)) + "\n\n".join(["{}".format(error) for value, error in errors]))
    
    def cancel(self):
        """
        Reset our method selector.
        """
        self.selector.browser.reset()    

    