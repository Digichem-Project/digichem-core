import urwid
from urwid import ListBox
from urwid.listbox import SimpleListWalker, SimpleFocusListWalker

class Calcbox_item(urwid.AttrMap):
    """
    """
    
    def __init__(self, method, calcbox):
        """
        Constructor for Calcbox_item objects.
        
        :param method: The method we represent.
        :param calcbox: The calcbox we are part of.
        """
        self.calcbox = calcbox
        self.method = method
        
        # Also resolve and store each of our method parts.
        # Each method consists of 3 parts: 1) destination, 2) program, 3) calculation.
        destination_path, program_path, calculation_path = self.method
        self.destination = destination_path[0].resolve_path(destination_path)
        self.program = program_path[0].resolve_path(program_path)
        self.calculation = calculation_path[0].resolve_path(calculation_path)
        
        # Setup our controls.
        self.up_button = urwid.Button("↑", lambda button: self.move_up())
        self.down_button = urwid.Button("↓", lambda button: self.move_down())
        self.modify_button = urwid.Button("m")
        self.delete_button = urwid.Button("r", lambda button: self.remove())
        
        # Keep track of our position widget (because we'll beed to update it).
        self.position_widget
        
        self.view = urwid.Pile([
            urwid.Columns([
                ('pack', urwid.AttrMap(urwid.Text(str("***")), "bold")),
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
        
        super().__init__(self.view, "body", "editable")

        
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
    
    
    def remove(self):
        """
        Remove this calcbox item from the parent calcbox.
        """
        index = self.calcbox.body.index(self)
        self.calcbox.body.pop(index)
        self.calcbox.update_positions()
        
    def update_position(self):
        """
        """
        # Get our index.
        try:
            position = str(self.position)
            
        except IndexError:
            position = "***"
        
        self.view.contents[0][0][0].base_widget.set_text(position)
    
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
        # We can access our list walker at the self.body attribute.
        super().__init__(SimpleFocusListWalker([]))
        
    def update_positions(self):
        """
        Tell each child to look up their positions again.
        """
        for child in self.body:
            child.update_position()
        
    