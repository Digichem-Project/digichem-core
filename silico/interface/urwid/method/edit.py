#from silico.interface.urwid.pages import Pages
from silico.interface.urwid.setedit.configurable import Configurable_browser
from silico.interface.urwid.layout import Pane
from silico.interface.urwid.setedit.base import Paginated_settings_browser


class Method_editor(Paginated_settings_browser):
    """
    A widget for editing different parts of a method.
    """
    
    def __init__(self, top, destination, program, calculation):
        """
        """
        calculation_pane = Pane(Configurable_browser.from_configurable(top, calculation, can_reset = False), self.get_title(calculation))
        program_pane = Pane(Configurable_browser.from_configurable(top, program, can_reset = False), self.get_title(program))
        destination_pane = Pane(Configurable_browser.from_configurable(top, destination, can_reset = False), self.get_title(destination))
        # We show in reverse order, because typically the calculation is what the user wants to change anyway.
        super().__init__({
            "Calculation": calculation_pane,
            "Program": program_pane,
            "Destination": destination_pane
        })
        
    def get_title(self, method_target):
        """
        Get the title to display for a given method_target.
        """
        try:
            return "Editing {} with ID: {}".format(method_target.TYPE, method_target.index())
        
        except IndexError:
            # No index.
            return "Editing {}".format(method_target.TYPE)