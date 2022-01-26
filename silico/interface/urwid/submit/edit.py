from silico.interface.urwid.pages import Pages
from silico.interface.urwid.setedit.configurable import Configurable_browser
from silico.interface.urwid.layout import Pane


class Method_editor(Pages):
    """
    A widget for editing different parts of a method.
    """
    
    def __init__(self, top, destination, program, calculation):
        """
        """
        self.calculation_pane = Pane(Configurable_browser(top, calculation), self.get_title(calculation))
        self.program_pane = Pane(Configurable_browser(top, program), self.get_title(program))
        self.destination_pane = Pane(Configurable_browser(top, destination), self.get_title(destination))
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
        self.calculation_editor.base_widget.discard()
        self.program_editor.base_widget.discard()
        self.destination_editor.base_widget.discard()
        
    def save(self):
        """
        Save any changes made.
        """
        self.calculation_editor.base_widget.save()
        self.program_editor.base_widget.save()
        self.destination_editor.base_widget.save()
        
    def get_title(self, method_target):
        """
        Get the title to display for a given method_target.
        """
        return "Editing {} with ID: {}".format(method_target.TYPE, method_target.index())