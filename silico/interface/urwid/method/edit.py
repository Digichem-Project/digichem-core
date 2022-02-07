# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.setedit.configurable import Configurable_browser
from silico.interface.urwid.layout import Pane
from silico.interface.urwid.setedit.base import Paginated_settings_browser,\
    Setedit_editor_mixin
from silico.interface.urwid.edit.popup import Choices_edit
from silico.interface.urwid.misc import Tab_pile
from silico.interface.urwid.dialogue import Confirm_or_cancel_dialogue
from silico.submit.calculation import Calculation_target
from silico.submit.program import Program_target
from silico.submit.destination import Destination_target


# TODO: Probably need a better name?
class Method_page(Tab_pile, Setedit_editor_mixin):
    """
    """
    
    def __init__(self, top, method_target, parent_class):
        """
        """
        self.top = top
        self.parent_class = parent_class
        
        # A widget which allows us to change the class of the configurable we're editing.
        self.class_widget = Choices_edit(top, parent_class.known_handles(), initial = method_target.class_name, title = "Class name", change_callback = self.popup_class_confirmation_dialogue)
        
        # The actual browser where we can change settings.
        self.browser = Configurable_browser.from_configurable(top, method_target, can_reset = False)
        self.browser_pane = Pane(self.browser, self.get_title(method_target))
        
        # Keep track of the method target we're editing, so we can rollback to it if we need.
        self.last_method_target = method_target
        
        super().__init__([("pack", Pane(self.class_widget, "Class")), self.browser_pane])
        
    def popup_class_confirmation_dialogue(self):
        """
        Show a dialogue that prompts the user if they're sure before changing the method target class.
        """
        # If the chosen class is the same as our current class, do nothing.
        if self.class_widget.value == self.browser.configurable.class_name:
            return
        
        def cancel_callback():
            self.class_widget.value = self.browser.configurable.class_name
        
        # Otherwise, show a confirmation box.
        self.top.popup(
            Confirm_or_cancel_dialogue("Change class?", "Are you certain that you wish to change class? This will discard all current options.",
                self.top,
                error = True,
                submit_callback = self.update_class,
                cancel_callback = cancel_callback 
            )
        )
        
    def update_class(self):
        """
        Create a new empty method_target from the current value of the class picker widget.
        """
        new_target = self.parent_class.from_class_handle(self.class_widget.value)(validate_now = False)
        self.update_current_method_target(new_target)
        
    def update_current_method_target(self, new_target):
        """
        Change the method target that is currently being edited.
        
        :param new_target: The method_target to set.
        """
        self.browser = Configurable_browser.from_configurable(self.top, new_target, can_reset = False)
        self.browser_pane.inner_body = self.browser
        self.class_widget.value = self.browser.configurable.class_name
        
    def get_title(self, method_target):
        """
        Get the title to display for a given method_target.
        """
        try:
            return "Editing {} with ID: {}".format(method_target.TYPE, method_target.index())
        
        except IndexError:
            # No index.
            return "Editing {}".format(method_target.TYPE)
        
    def refresh(self):
        self.browser.refresh()
        
    def save(self, validate=True):
        self.browser.save(validate)
        # Save was successful, update the rollback point.
        self.last_method_target = self.browser.configurable
        
    def discard(self):
        # If the current method target that's being edited is not the same as the last one that was saved, we'll rollback the entire configurable.
        if self.last_method_target is not self.browser.configurable:
            # Reset.
            self.update_current_method_target(self.last_method_target)
            
        else:
            self.browser.discard()
        
    def validate(self):
        self.browser.validate()


class Method_editor(Paginated_settings_browser):
    """
    A widget for editing different parts of a method.
    """
    
    def __init__(self, top, destination, program, calculation, on_change_callback = None):
        """
        """
        self.top = top
        calculation_pane = Method_page(top, calculation, Calculation_target)
        program_pane = Method_page(top, program, Program_target)
        destination_pane = Method_page(top, destination, Destination_target)
        # We show in reverse order, because typically the calculation is what the user wants to change anyway.
        super().__init__({
            "Calculation": calculation_pane,
            "Program": program_pane,
            "Destination": destination_pane
        }, on_change_callback = on_change_callback)
        
    @property
    def method(self):
        return (
            self.pages['Destination'].browser.configurable,
            self.pages['Program'].browser.configurable,
            self.pages['Calculation'].browser.configurable,
        )
        
        