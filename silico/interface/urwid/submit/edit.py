from silico.interface.urwid.pages import Pages
from silico.interface.urwid.setedit.base import Settings_editor
from silico.interface.urwid.setedit.configurable import Configurable_browser


class Method_editor(Pages):
    """
    A widget for editing different parts of a method.
    """
    
    def __init__(self, destination, program, calculation):
        """
        """
        self.calculation_editor = Settings_editor(Configurable_browser(calculation), self.get_title(calculation))
        self.program_editor = Settings_editor(Configurable_browser(program), self.get_title(program))
        self.destination_editor = Settings_editor(Configurable_browser(destination), self.get_title(destination))
        # We show in reverse order, because typically the calculation is what the user wants to change anyway.
        super().__init__([
            ("Calculation", self.calculation_editor),
            ("Program", self.program_editor),
            ("Destination", self.destination_editor)
        ])
        
    def discard(self):
        """
        Reset any changed values back to their defaults.
        """
        self.calculation_editor.browser.discard()
        self.program_editor.browser.discard()
        self.destination_editor.browser.discard()
        
    def save(self):
        """
        Save any changes made.
        """
        self.calculation_editor.browser.save()
        self.program_editor.browser.save()
        self.destination_editor.browser.save()
        
    def get_title(self, method_target):
        """
        Get the title to display for a given method_target.
        """
        return "Editing {} with ID: {}".format(method_target.TYPE, method_target.index())