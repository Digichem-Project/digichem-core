# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.row_list.base import Row_item, Row_widget,\
    Row_browser
from silico.interface.urwid.method.browser import Method_selector
from silico.interface.urwid.submit.edit import Method_editor
import shlex
import logging
import silico
from silico.config.configurable.loader import Partial_loader


class Method_widget(Row_widget):
    """
    Widget used to display a method in a Method_list.
    """
    
    def __init__(self, row_item):
        Row_widget.__init__(self, row_item)
        self.method_editor = Method_editor(self.row_item.destination, self.row_item.program, self.row_item.calculation)

    def get_code(self):
        """
        Get the unique code of the method we represent.
        """
        return "[{}/{}/{}]".format(
            self.row_item.destination.index(), 
            self.row_item.program.index(),
            self.row_item.calculation.index()
        )
    
    def get_text(self):
        """
        Get the text that we will display.
        """
        return "{} : {} : {}".format(
            self.row_item.destination.name,
            self.row_item.program.name,
            self.row_item.calculation.name,
        )
    
    def load_inner(self):
        """
        Load the inner widget we'll use for display.
        """
        return urwid.Columns([
            ("pack", urwid.Text(self.get_code())),
            urwid.Text(self.get_text())
        ], dividechars = 1)
        
    def load_controls(self):
        """
        Load the buttons that can be used to control this row.
        
        This implementation returns up and down controls (if we are movable) and a remove control.
        
        :returns: A list of control buttons.
        """
        up_button = urwid.Button("↑", lambda button: self.row_item.move_up())
        down_button = urwid.Button("↓", lambda button: self.row_item.move_down())
        modify_button = urwid.Button("m", lambda button: self.switch_to_editer())
        delete_button = urwid.Button("r", lambda button: self.row_item.remove())
        
        controls = []
        if self.movable:
            controls.append(urwid.AttrMap(up_button, "button--small", "button--small--focus"))
            controls.append(urwid.AttrMap(down_button, "button--small", "button--small--focus"))
        
        controls.append(urwid.AttrMap(modify_button, "button--small", "button--small--focus"))
        controls.append(urwid.AttrMap(delete_button, "button--bad--small", "button--small--focus"))
        return controls
    
    def switch_to_editer(self):
        """
        Open a window where settings for this method can be changed.
        """
        self.row_item.row_list.top.swap_into_window(self.method_editor, cancel_callback = self.method_editor.discard, submit_callback = lambda: print("oooh"))


class Method_item(Row_item):
    """
    Logical storage for methods.
    """
    
    def __init__(self, method, row_list, movable=True):
        """
        Constructor for Method_item objects.
        
        :param method: The method we represent, a tuple of configurables (destination, program, calculation).
        :param row_list: The parent Row_list object.
        :param movable: Whether we are movable.
        """
        Row_item.__init__(self, method, row_list, movable=movable)
        # Also resolve and store each of our method parts.
        # Each method consists of 3 parts: 1) destination, 2) program, 3) calculation.
        #destination_path, program_path, calculation_path = method
        #self.destination = destination_path[0].resolve_path(destination_path)
        #self.program = program_path[0].resolve_path(program_path)
        #self.calculation = calculation_path[0].resolve_path(calculation_path)
        self.destination, self.program, self.calculation = method
        
    def load_widget(self):
        """
        Load the widget we'll use for display.
        """
        return Method_widget(self)
    
    def get_value(self):
        """
        Get the current value of this coordinate item object.
        """
        return (self.destination, self.program, self.calculation)


class Method_list(Row_browser):
    """
    A widget for displaying a list of loaded methods.
    """
    
    def __init__(self, top, methods, rearrangeable = True, initial_methods = None):
        """
        Constructor for Method_list objects.
        
        :param top: The topmost widget to use for display.
        :param methods: A Configurable_list or similar that contains available methods.
        :param rearrangeable: Whether the order of chosen methods can be rearranged.
        :param initial_methods: A list of methods (tuples) to initially populate with.
        """
        # Add our starting files.
        initial_methods = [] if initial_methods is None else initial_methods
        
        super().__init__(Method_selector(methods), top, rearrangeable = rearrangeable, initial = initial_methods)
        
    def add_to_list(self):
        """
        Add the items currently selected in our browser to our list.
        """
        # Our parent will take care of adding items from the tree browser.
        super().add_to_list()
        
        # Get logger for errors.
        logger = logging.getLogger(silico.logger_name)
        
        # Now add the methods identified by their code.
        method_strings = shlex.split(self.selector.manual_widget.get_edit_text())
        
        # Parse each string.
        for method_string in method_strings:
            try:
                self.add_row(self.selector.browser.methods.resolve_method_string(method_string))
                
            except Exception:
                logger.warning("Failed to add method from code '{}':".format(method_string), exc_info = True)
                
        # Reset our selector.
        self.selector.reset()
            
    
    def value_from_node(self, node):
        """
        Get a value to add from a selected node.
        """
        destination_path, program_path, calculation_path = node.build_loader_path()
        
        destination = destination_path[0].resolve_path(destination_path)
        program = program_path[0].resolve_path(program_path)
        calculation = calculation_path[0].resolve_path(calculation_path)
        
        return (destination, program, calculation)
    
    def load_row(self, value):
        return Method_item(value, self, self.rearrangeable)

