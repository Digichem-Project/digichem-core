# General import.
import urwid
from urwid.font import HalfBlock7x7Font

# Silico import.
import silico
from silico.interface.urwid.window import Window
from silico.interface.urwid.section import Section
from silico.interface.urwid.top import View
import io


class Silico_window(Window):
    """
    Main interactive interface for silico.
    """
    
    def __init__(self, programs = None, initial = None):
        """
        Constructor for Interface objects.
        
        :param programs: A list of programs that we can switch between.
        :param initial: A program to show initially.
        """
        # Save our sub programs.
        self.programs = programs
        
        # Setup our window.
        super().__init__(title = "Silico {}".format(silico.version))
        
        body = urwid.Pile([
            ('pack', urwid.Padding(urwid.BigText("Silico", HalfBlock7x7Font()), align = "center", width = "clip")),
            ('pack', self.get_menu())
        ])
        #self.top.original_widget = urwid.Filler(body)
        #self.top.set_top(urwid.Filler(body))
        self.top.swap(urwid.Filler(body))
        
        # If we've got an initial program, swap to it now,
        if initial is not None:
            self.top.swap_into_window(initial.get_interface(self))
        
    def get_menu(self):
        """
        Get the main menu from which the user can select an option.
        """
        return Section(urwid.Pile(self.get_options()), "Silico Main Menu")
        
    def get_options(self):
        """
        Get the list of options (buttons) that we will display in our main menu.
        
        :returns: A list of widgets.
        """
        program_buttons = [self.program_button(program) for program in self.programs]
        # Add our exit button.
        program_buttons.append(self.get_option_widget("Quit", lambda button: self.top.back()))
        
        return program_buttons
        
    def get_option_widget(self, name, callback):
        """
        Get a widget to display an option from the main menu.
        
        :param name: The name of the opion/button.
        :param callback: A function to call when the button is selected.
        """
        # width = "pack" doesn't seem to work as expected for buttons?
        return urwid.AttrMap(urwid.Padding(urwid.Button(name, callback), width = "pack"), "button--small", "button--small--focus")
    
    def program_button(self, program):
        """
        Get an option that will populate our main menu.
        
        :param name: The name of the option.
        :param widget: A widget to switch to when this option is chosen.
        :returns: A widget.
        """
        return self.get_option_widget(program.command, lambda button: self.top.swap_into_window(program.get_interface(self))) 
    
    def run_loop(self, palette):
        """
        Run an urwid loop using this Window as the top-most frame.
        
        :param palette: The palette to display with.
        """
        super().run_loop(palette)
        

class Program_view(View):
    """
    A widget for running a sub-program.
    """
    
    def __init__(self, window, program):
        """
        Constructor for Program_view objects.
        
        :param window: A Silico_window widget used for display.
        :param program: A program object to run.
        """
        self.window = window
        self.program = program
        super().__init__(self.get_body(), title = program.name, border = False)
        
    def get_body(self):
        """
        Get the widget used to display the body of this program.
        
        :returns: An urwid widget to display.
        """
        raise NotImplementedError("Implement in subclass")
    
    def setup(self):
        """
        Setup our program to ready it to be run.
        """
        # This default implementation does nothing.
    
    def submit(self):
        """
        Submit this program widget, running the program it wraps.
        """
        self.setup()
        retval = self.program.main()
        self.post(retval)
        
        return retval
        
    def post(self):
        """
        Method called once our main program has finished running.
        """
        # This default implementation does nothing.


class Output_catcher(io.StringIO):
    """
    A class for catching output sent to stdout and stderr and redirecting it to an urwid widget.
    """
    
    def __init__(self, top, error = False):
        """
        """
        self.top = top
        self.error = error
        super().__init__()
        
    def write(self, *args, **kwargs):
        retval = super().write(*args, **kwargs)
        
        # If the last char in our buffer is a newline, we can write to out widget.
        curvalue = self.getvalue()
        if curvalue[-1] == "\n":
            # Update our top.
            self.top.output(self.getvalue(), self.error)
            
            # Then empty our buffer (could be a bad idea?)
            self.seek(0)
            self.truncate(0)
        
        return retval
        
        