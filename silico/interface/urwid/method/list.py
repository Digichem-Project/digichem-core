# General imports.
import urwid
import shlex

# Silico imports.
from silico.interface.urwid.row_list.base import Row_item, Row_widget,\
    Row_browser, Row_pointer, Row_pointer_widget
from silico.interface.urwid.method.browser import Method_selector
from silico.interface.urwid.method.edit import Method_editor
import silico.logging
from silico.interface.urwid.dialogue import Edit_dialogue
from silico.interface.urwid.file.browser import File_selector
from silico.submit.base import parse_method_from_file


class Method_widget(Row_widget):
    """
    Widget used to display a method in a Method_list.
    """
    
    def __init__(self, row_item):
        Row_widget.__init__(self, row_item)
        self.method_editor = Method_editor(self.row_item.row_list.top, self.row_item.destination, self.row_item.program, self.row_item.calculation)

    def get_code(self):
        """
        Get the unique code of the method we represent.
        """
        # Carefully get the code of each part of our method, being aware that methods not loaded from the library don't have an index.
        try:
            dest_code = self.row_item.destination.index()
        
        except IndexError:
            dest_code = "N.A"
            
        try:
            prog_code = self.row_item.program.index()
        
        except IndexError:
            prog_code = "N.A"
            
        try:
            calc_code = self.row_item.calculation.index()
        
        except IndexError:
            calc_code = "N.A"
        
                
        return "[{}/{}/{}]".format(dest_code, prog_code, calc_code)
    
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
        self.row_item.row_list.top.swap_into_window(self.method_editor, cancel_callback = self.method_editor.discard, submit_callback = self.method_editor.confirm_callback)


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


class Method_pointer_widget(Row_pointer_widget):
    """
    Specific list pointer widget used by method lists; allows specifying methods via either a browser or directly by a code.
    """
    
    def __init__(self, row_item):
        self.code_popup = Edit_dialogue("Method codes", row_item.row_list.top, submit_callback = self.add_from_edit_code)
        self.file_selector = File_selector(row_item.row_list.top, pane_title = "Select method file(s)", can_choose_folders = False, can_choose_multiple = True)
        super().__init__(row_item)
        
    def add_from_edit_code(self):
        """
        Add the values from our edit widget to our parent list.
        """
        # This reference is long, perhaps there is better access?
        methods = self.row_item.row_list.selector.browser.method_library.methods
        
        # Add each of the methods identified by each code.
        for code in self.code_popup.edit.edit_text.split():
            try:
                self.row_item.row_list.add_row(methods.resolve_method_string(code))
            
            except Exception:
                # Something went wrong, most probably the given code is no good.
                silico.logging.get_logger().error("Failed to resolve method identified by code '{}'".format(code), exc_info = True)
        
        # Clear our edit widget.
        self.code_popup.edit.edit_text = ""
        
    def add_from_files(self):
        """
        Add the values from our file selector to our parent list.
        """
        for method_file in self.file_selector.browser.selected:
            try:
                self.row_item.row_list.add_row(parse_method_from_file(method_file, self.row_item.row_list.selector.browser.method_library))
            
            except Exception:
                # Something went wrong, most probably the given code is no good.
                silico.logging.get_logger().error("Failed to parse method file '{}'".format(method_file), exc_info = True)
                
        self.file_selector.browser.reset()
    
    def load_inner(self):
        """
        Load the widget we'll use to display our main body.
        """
        return urwid.Columns([
            urwid.Button("Browse library", lambda button: self.row_item.browser_callback()),
            urwid.Button("Add from code", lambda button: self.row_item.row_list.top.popup(self.code_popup)),
            urwid.Button("Add from files", lambda button: self.row_item.row_list.top.swap_into_window(self.file_selector, submit_callback = self.add_from_files))
        ], dividechars = 1)


class Method_pointer(Row_pointer):
    def load_widget(self):
        return Method_pointer_widget(self)


class Method_list(Row_browser):
    """
    A widget for displaying a list of loaded methods.
    """
    
    def __init__(self, top, method_library, rearrangeable = True, initial_methods = None):
        """
        Constructor for Method_list objects.
        
        :param top: The topmost widget to use for display.
        :param method_library: An object that contains methods, destinations, programs and calculations as attributes.
        :param rearrangeable: Whether the order of chosen methods can be rearranged.
        :param initial_methods: A list of methods (tuples) to initially populate with.
        """
        # Add our starting files.
        initial_methods = [] if initial_methods is None else initial_methods
        
        super().__init__(Method_selector(top, method_library), top, rearrangeable = rearrangeable, initial = initial_methods)
        
    def get_row_pointer(self, rearrangeable):
        return Method_pointer(self, self.swap_to_browser, rearrangeable)
    
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

