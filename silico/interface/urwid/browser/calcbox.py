import urwid
from urwid import ListBox
from urwid.listbox import SimpleFocusListWalker

class Calcbox_item(urwid.AttrMap):
    """
    """
    
    def __init__(self, calcbox):
        """
        Constructor for Calcbox_item objects.
        
        :param method: The method we represent.
        :param calcbox: The calcbox we are part of.
        """
        self.calcbox = calcbox
        
        # Setup our controls.
        # Not all calcbox items will use all these controls, but we set them up anyway.
        self.up_button = urwid.Button("↑", lambda button: self.move_up())
        self.down_button = urwid.Button("↓", lambda button: self.move_down())
        self.modify_button = urwid.Button("m")
        self.delete_button = urwid.Button("r", lambda button: self.remove())
        
        # Keep track of our position widget (because we'll need to update it).
        self.position_widget = urwid.Text("***")
        
        self.view = self.get_view()
        
        super().__init__(self.view, "body", "editable")
        
    def get_view(self):
        """
        Remove this calcbox item from the parent calcbox.
        
        Should be implemented in inheriting classes.
        """
        raise NotImplementedError()

        
    @property
    def position(self):
        """
        Our position in the calcbox.
        """        
        return self.calcbox.body.index(self) +1
    
    def move_up(self):
        """
        """
        return self.move('up')
    
    def move_down(self):
        """
        """
        return self.move('down')
        
    def move(self, position):
        """
        Move this calcbox item either up or down.
        
        :param position: One of either 'up' or 'down'.
        """
        # Get our current index.
        cur_index = self.calcbox.body.index(self)
        
        if position == "up":
            new_index = cur_index -1
            
        elif position == "down":
            new_index = cur_index +1
            
        else:
            raise ValueError("position must be one of either 'up' or 'down'")
        
        # Don't allow negatives
        if new_index < 0:
            new_index = 0
        
        # We need to insert and pop to move.
        self.calcbox.body.insert(new_index, self.calcbox.body.pop(cur_index))
        
        # Update.
        self.calcbox.update_positions()
        
        # Set ourself as the focus.
        self.calcbox.body.set_focus(new_index if new_index < len(self.calcbox.body) else len(self.calcbox.body) -1)
        
    def update_position(self):
        """
        Update our position widget with our current position in the parent calcbox.
        """
        # Get our index.
        try:
            position = str(self.position)
            
        except IndexError:
            position = "***"
        
        self.position_widget.set_text(position)


class Calcbox_pointer(Calcbox_item):
    """
    A calcbox item that indicates the position at which the next method will be appended.
    """
    
    def get_view(self):
        """
        Get the inner widget that we'll use for display.
        """
        return urwid.Pile([
            urwid.Columns([
                ('pack', urwid.AttrMap(self.position_widget, "bold")),
                urwid.Text(self.get_text()),
                (10, urwid.GridFlow([
                    urwid.AttrMap(self.up_button, "normalButton"),
                    urwid.AttrMap(self.down_button, "normalButton")
                ], 5, 0, 0, "center"))
            ], dividechars = 2),
            urwid.Divider()
        ])

    def get_text(self):
        """
        Get the text that we will display.
        """
        return "<<< Next method inserted here >>>"

class Calcbox_method(Calcbox_item):
    """
    An item in a calbox that represents a calculation method.
    """
    
    def __init__(self, method, calcbox):
        """
        Constructor for Calcbox_item objects.
        
        :param method: The method we represent.
        :param calcbox: The calcbox we are part of.
        """
        self.method = method
        
        # Also resolve and store each of our method parts.
        # Each method consists of 3 parts: 1) destination, 2) program, 3) calculation.
        destination_path, program_path, calculation_path = self.method
        self.destination = destination_path[0].resolve_path(destination_path)
        self.program = program_path[0].resolve_path(program_path)
        self.calculation = calculation_path[0].resolve_path(calculation_path)
        
        super().__init__(calcbox)
        
    def get_view(self):
        """
        Get the inner widget that we'll use for display.
        """
        return urwid.Pile([
            urwid.Columns([
                ('pack', urwid.AttrMap(self.position_widget, "bold")),
                ('pack', urwid.AttrMap(urwid.Text(("bold", self.get_code())), "bold")),
                urwid.Text(self.get_text()),
                (10, urwid.GridFlow([
                    urwid.AttrMap(self.up_button, "normalButton"),
                    urwid.AttrMap(self.down_button, "normalButton"),
                    urwid.AttrMap(self.modify_button, "normalButton"),
                    urwid.AttrMap(self.delete_button, "badButton--small")
                ], 5, 0, 0, "center"))
            ], dividechars = 2),
            urwid.Divider()
        ])
    
    def remove(self):
        """
        Remove this calcbox item from the parent calcbox.
        """
        index = self.calcbox.body.index(self)
        self.calcbox.body.pop(index)
        self.calcbox.update_positions()
    
    def get_code(self):
        """
        Get the unique code of the method we represent.
        """
        return "[{}/{}/{}]".format(
            self.destination.index(), 
            self.program.index(),
            self.calculation.index()
        )
        
    def get_text(self):
        """
        Get the text that we will display.
        """
        return "{} : {} : {}".format(
            self.destination.name,
            self.program.name,
            self.calculation.name,
        )

class Calcbox(ListBox):
    """
    A widget that displays a list of selected calculations.
    """
    
    def __init__(self):
        """
        Constructor for Calcbox objects.
        """
        # Keep track of our pointer object.
        self.pointer = Calcbox_pointer(self)
        
        # We can access our list walker at the self.body attribute.
        super().__init__(SimpleFocusListWalker([self.pointer]))
        
        # Set our initial focus.
        self.set_focus(0)
        
        self.update_positions()
        
    def update_positions(self):
        """
        Tell each child to look up their positions again.
        """
        for child in self.body:
            child.update_position()
            
    def add_method(self, method):
        """
        Add a new method to the calcbox.
        """
        calcbox_item = Calcbox_method(method, self)
        
        # Find the current position of our pointer.
        index = self.body.index(self.pointer)
        
        self.body.insert(index, calcbox_item)
        self.update_positions()
        
    