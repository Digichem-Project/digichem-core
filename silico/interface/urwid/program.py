import urwid
from urwid.font import HalfBlock7x7Font
import io

import silico.log
from silico.interface.urwid.setedit.configurable import make_paginated_configurable_browser
from silico.interface.urwid.layout import Window, Pane, Sub_pane
from silico.program.status import Status_program
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.dialogue import Confirm_dialogue
from silico.author import get_authorship_string
from silico.interface.urwid.method.builder import Method_builder_menu
from silico.interface.urwid.wrapper import Cancel
from silico.interface.urwid.swap.swappable import Swappable
from silico.interface.urwid.swap.top import Top


class Silico_window(Window):
    """
    Main interactive interface for silico.
    """
    
    def __init__(self, program):
        """
        Constructor for Interface objects.
        
        :param program: The interactive program object we represent.
        """
        self.program = program
        
        # A dictionary of program interfaces that we can display.
        # Each key is an instance of each sub program.
        self._interfaces = {}
        
        # The top-most widget, a funky object that allows swapping back and forth between which widget is currently shown.
        top = Top(None)
        
        # A widget which can be swapped to in order to change the main silico settings.        
        self.settings_pane = Sub_pane(make_paginated_configurable_browser(program.config, top, on_change_callback = self.update_settings, page_selector_title = "Settings Type"), "Main Silico Settings")
        
        # A widget for creating new methods.
        self.method_builder = Method_builder_menu(top, program.config)
        
        body = Tab_pile([
            ('pack', urwid.Padding(urwid.BigText("Silico", HalfBlock7x7Font()), align = "center", width = "clip")),
            ('pack', urwid.Padding(urwid.AttrMap(urwid.Button("Silico development team", lambda button: self.popup_authorship_window()), "body", "editable"), align = "center", width = 27)),#19
            ('pack', urwid.Padding(self.get_menu(), align = "center", width = 45))
        ], focus_item = 2)
        
        top.swap(urwid.Filler(body))
        
        # Setup our window.
        super().__init__(top, title = "Silico {}".format(silico.version), help = "Arrow Keys: Navigate  TAB: Next: SHIFT-TAB: Previous  ENTER: Select  ESC: Back")
        
        # If we've got an initial program, swap to it now,
        if self.program.initial is not None:
            self.top.swap_into_window(self.get_interface(self.program.initial))
            
    def update_settings(self, confirm):
        """
        Method called when main silico options are changed.
        
        :param confirm: Whether the settings were updated because they were changed (True) or because they were rolled back (False).
        """
        # Update logging.
        silico.log.set_logging_level(self.program.config['logging']['log_level'], self.program.config['logging']['verbose'])
        
        # Save modified settings to file.
        if confirm:
            self.program.config.save()
            
            
    def get_interface(self, prog_obj):
        """"
        Get an instance of a sub program's interface.
        
        This method will first create a new instance of the program class before caching it so future calls will use the same instance.
        Likewise the returned interface widget will also be saved to a cache, so subsequent calls will return the same interface.
        
        :param prog_obj: The program class to get.
        :returns: A program instance.
        """
        if prog_obj not in self._interfaces:
            self._interfaces[prog_obj] =  prog_obj.get_interface(self)
            
        return self._interfaces[prog_obj]
        
    def get_menu(self):
        """
        Get the main menu from which the user can select an option.
        """
        return Pane(urwid.Pile(self.get_options()), "Main Menu")
        
    def get_options(self):
        """
        Get the list of options (buttons) that we will display in our main menu.
        
        :returns: A list of widgets.
        """
        program_buttons = [self.program_button(prog_cls) for prog_cls in self.program.program_classes if prog_cls is not Status_program]
        # Add a status button.
        # TODO: Should add a smarter interface for the status program.
        program_buttons.append(self.get_option_widget("Status: Check queue status", lambda button: self.program.get_program(Status_program).main()))
        # Add a method builder button.
        program_buttons.append(self.get_option_widget("Method builder", lambda button: self.swap_to_method_builder()))
        # Add a settings button.
        program_buttons.append(self.get_option_widget("Settings", lambda button: self.swap_to_main_settings()))
        # And an exit button.
        program_buttons.append(self.get_option_widget("Quit", lambda button: self.top.back()))
        
        return program_buttons
    
    def popup_authorship_window(self):
        """
        """
        self.top.popup(Confirm_dialogue("The Silico Development Team", get_authorship_string(), self.top), width = ('relative', 60), height = ('relative', 60))
        
    def swap_to_method_builder(self):
        """
        Swap to the widget that can be used to create new methods.
        """
        self.top.swap(Cancel(self.method_builder, self.top, lambda: None))
    
    def swap_to_main_settings(self):
        """
        Swap to the settings editor widget.
        """
        # TODO: Calling refresh should really happen elsewhere (mostly because it's easy to forget if its the job of the swapping widget to refresh).
        self.settings_pane.base_widget.refresh()
        self.top.swap_into_window(self.settings_pane, cancel_callback = self.settings_pane.base_widget.cancel_callback, submit_callback = self.settings_pane.base_widget.confirm_callback)
        
    def get_option_widget(self, name, callback):
        """
        Get a widget to display an option from the main menu.
        
        :param name: The name of the opion/button.
        :param callback: A function to call when the button is selected.
        """
        # width = "pack" doesn't seem to work as expected for buttons?
        return urwid.AttrMap(urwid.Padding(urwid.Button(name, callback), width = "pack"), "button--small", "button--small--focus")
    
    def program_button(self, prog_cls):
        """
        Get an option that will populate our main menu.
        
        
        :returns: A widget.
        """
        return self.get_option_widget("{}: {}".format(prog_cls.command.capitalize(), prog_cls.help), lambda button: self.top.swap_into_window(self.get_interface(self.program.get_program(prog_cls)))) 
    
    def run_loop(self, palette):
        """
        Run an urwid loop using this Window as the top-most frame.
        
        :param palette: The palette to display with.
        """
        super().run_loop(palette)


class Output_catcher(io.StringIO):
    """
    A class for catching output sent to stdout and stderr and redirecting it to an urwid widget.
    """
    
    def __init__(self, window, error = False):
        """
        """
        self.window = window
        self.error = error
        super().__init__()
        
    def write(self, s):
        retval = super().write(s)
        
        # If the last char in our buffer is a newline, we can write to our widget.
        curvalue = self.getvalue()
        if len(curvalue) > 0 and curvalue[-1] == "\n":
            # Update our top.
            self.window.top.output(curvalue, self.error)
            
            # Update our widgets (also could be a bad idea?)
            self.window.loop.draw_screen()
            
            # Then empty our buffer (could be a bad idea?)
            self.seek(0)
            self.truncate(0)
        
        return retval


class Program_view(Swappable):
    """
    A widget for running a subprogram.
    """
    
    def __init__(self, window, program):
        """
        Constructor for Program_view objects.
        
        :param window: A Silico_window widget used for display.
        :param program: A program object to run.
        """
        self.window = window
        self.program = program
        super().__init__(self.window.top, Sub_pane(self.get_body(), title = program.meta['name']))
        
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
        
        :returns: True if program execution was successful, false otherwise.
        """
        try:
            self.setup()
            self.program.main()
            self.post()
            #return True
            
        except Exception:
            logger = silico.log.get_logger()
            logger.exception("Sub-program {} stopped with error".format(self.program.command))
            #return False
        
        except KeyboardInterrupt:
            # The user wanted to stop.
            logger = silico.log.get_logger()
            logger.error("Sub-program {} interrupted by user (ctrl-c)".format(self.program.command))
            #return False
        
        # Return false to prevent swapping back.
        return False
        
    def post(self):
        """
        Method called once our main program has finished running.
        """
        # This default implementation does nothing.
        