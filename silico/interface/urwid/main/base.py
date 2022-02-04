# General import.
import urwid
from urwid.font import HalfBlock7x7Font
import io

# Silico import.
import silico
import silico.logging
from silico.interface.urwid.setedit.configurable import make_paginated_configurable_browser
from silico.interface.urwid.layout import Window, Pane, Sub_pane
from silico.program.status import Status_program
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.dialogue import Confirm_dialogue
from silico.authors import get_authorship_string


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
        
        # Setup our window.
        super().__init__(title = "Silico {}".format(silico.version), help = "Arrow Keys: Navigate  TAB: Next: SHIFT-TAB: Previous  ENTER: Select  ESC: Back")
        
        # A widget which can be swapped to in order to change the main silico settings.        
        self.settings_pane = Sub_pane(make_paginated_configurable_browser(self.program.config, self.top, on_change_callback = self.update_settings, page_selector_title = "Settings Type"), "Main Silico Settings")
        
        body = Tab_pile([
            ('pack', urwid.Padding(urwid.BigText("Silico", HalfBlock7x7Font()), align = "center", width = "clip")),
            ('pack', urwid.Padding(urwid.AttrMap(urwid.Button("Silico development team", lambda button: self.popup_authorship_window()), "body", "editable"), align = "center", width = 27)),#19
            ('pack', urwid.Padding(self.get_menu(), align = "center", width = 45))
        ], focus_item = 2)
        
        self.top.swap(urwid.Filler(body))
        
        # If we've got an initial program, swap to it now,
        if self.program.initial is not None:
            self.top.swap_into_window(self.get_interface(self.program.initial))
            
    def update_settings(self):
        """
        Method called when main silico options are changed.
        """
        # Update logging.
        silico.logging.set_logging_level(self.program.config['logging']['log_level'], self.program.config['logging']['verbose'])
        
        # Save modified settings to file.
        self.program.config.save()
        
        # Clear our cache of programs (so when they're called again they'll be recreated and take the updated options.
        #self.program._programs = {}
            
            
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
        # Add a settings button.
        program_buttons.append(self.get_option_widget("Settings: Permanently change settings", lambda button: self.swap_to_main_settings()))
        # And an exit button.
        program_buttons.append(self.get_option_widget("Quit", lambda button: self.top.back()))
        
        return program_buttons
    
    def popup_authorship_window(self):
        """
        """
        self.top.popup(Confirm_dialogue("The Silico Development Team", get_authorship_string(), self.top), width = ('relative', 60), height = ('relative', 60))
    
    def swap_to_main_settings(self):
        """
        Swap to the settings editor widget.
        """
        # TODO: Calling refresh should really happen elsewhere (mostly because it's easy to forget if its the job of the swapping widget to refresh).
        self.settings_pane.base_widget.refresh()
        self.top.swap_into_window(self.settings_pane, cancel_callback = self.settings_pane.base_widget.discard, submit_callback = self.settings_pane.base_widget.confirm_callback)
        
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
        if curvalue[-1] == "\n":
            # Update our top.
            self.window.top.output(curvalue, self.error)
            
            # Update our widgets (also could be a bad idea?)
            self.window.loop.draw_screen()
            
            # Then empty our buffer (could be a bad idea?)
            self.seek(0)
            self.truncate(0)
        
        return retval
        
        