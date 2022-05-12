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


class Method_target_editor(Tab_pile, Setedit_editor_mixin):
    """
    A widget that edits part of a method (a calculation, program or destination)
    """
    
    def __init__(self, top, parent_class, method_target = None):
        """
        Constructor for Method_page objects.
        
        :param method_target: The method target (one part of a method: a calculation, program or destination) to edit.
        :param parent_class: The parent class of method_target. This class should inherit from Dynamic_parent and is used to change the class of method_target that is being edited.
        """
        self.top = top
        self.parent_class = parent_class
        
        # A widget which allows us to change the class of the configurable we're editing.
        self.class_widget = Choices_edit(top, parent_class.known_handles(), initial = method_target.class_name if method_target is not None else None, title = "Class name", change_callback = self.popup_class_confirmation_dialogue)
        
        # The actual browser where we can change settings.
        # If we do not have a configurable class to start with, this will be an empty placeholder for now.
        if method_target is None:
            self.browser = None
            self.browser_swapper = urwid.WidgetPlaceholder(urwid.Filler(urwid.Text("")))
            
        else:
            self.browser = Configurable_browser.from_configurable(top, method_target, can_reset = False)
            self.browser_swapper = urwid.WidgetPlaceholder(
                Pane(self.browser, method_target.TYPE)
            )
        
        
        # Keep track of the method target we're editing, so we can rollback to it if we need.
        self.last_method_target = method_target
        
        super().__init__(self.get_body())
        
    @property
    def configurable(self):
        """
        The configurable that is being edited.
        
        Note that thid property can return None if no initial method_target was given and the user has not yet chosen a configurable class.
        """
        try:
            return self.browser.configurable
        
        except AttributeError:
            if self.browser is None:
                return None
            
            else:
                raise
        
    def get_body(self):
        return [
            ("pack", Pane(self.class_widget, "Select class")),
            self.browser_swapper
        ]
        
    def popup_class_confirmation_dialogue(self):
        """
        Show a dialogue that prompts the user if they're sure before changing the method target class.
        """
        if self.browser is None:
            self.update_class()
            return
        
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
        new_target = self.parent_class.from_class_handle(self.class_widget.value)(class_name = self.class_widget.value, validate_now = False)
        self.update_current_method_target(new_target)
        
    def update_current_method_target(self, new_target):
        """
        Change the method target that is currently being edited.
        
        :param new_target: The method_target to set.
        """
        self.browser = Configurable_browser.from_configurable(self.top, new_target, can_reset = False)
        self.browser_swapper.original_widget = Pane(self.browser, new_target.TYPE)
        self.class_widget.value = self.browser.configurable.class_name
        
    def refresh(self):
        # If no class has yet been selected, do nothing.
        if self.class_widget.value is None:
            return
        
        self.browser.refresh()
        
    def save(self, validate=True):
        """
        Save (or attempt to) any changes made since the last made.
        """
        # If no class has yet been selected, do nothing.
        if self.class_widget.value is None:
            return
        
        self.browser.save(validate)
        # Save was successful, update the rollback point.
        self.last_method_target = self.browser.configurable
        # Call finalize on the target, so child objects created from its template will see the new changes.
        # TODO: Perhaps there should be an option to call this from a Configurable_browser object?
        self.browser.configurable.finalize()
        
    def discard(self, save = True):
        """
        Discard any changes made since the last save.
        
        :param save: Whether to save (and validate) the rolled back changes. Normally this is safe, but if this browser is only editing part of the options of a configurable and some of those other options have not been reset (or otherwise contain invalid options), this method will fail (because all the options have to be validated, not just those edited by this browser).
        """
        # If no class has yet been selected, do nothing.
        if self.class_widget.value is None:
            return
        
        # If the current method target that's being edited is not the same as the last one that was saved, we'll rollback the entire configurable.
        if self.last_method_target is not self.browser.configurable:
            # Reset.
            self.update_current_method_target(self.last_method_target)
            
        else:
            self.browser.discard(save = save)
        
    def validate_setedits(self):
        # If no class has yet been selected, do nothing.
        if self.class_widget.value is None:
            return
        
        self.browser.validate_setedits()


class Method_editor(Paginated_settings_browser):
    """
    A widget for editing different parts of a method.
    """
    
    def __init__(self, top, destination, program, calculation, on_change_callback = None):
        """
        """
        self.top = top
        calculation_page = self.get_page(top, Calculation_target, calculation)
        program_page = self.get_page(top, Program_target, program)
        destination_page = self.get_page(top, Destination_target, destination)
        # We show in reverse order, because typically the calculation is what the user wants to change anyway.
        super().__init__({
            "Calculation": calculation_page,
            "Program": program_page,
            "Destination": destination_page
        }, on_change_callback = on_change_callback)
        
    def get_page(self, *args, **kwargs):
        return Method_target_editor(*args, **kwargs)
        
    @property
    def method(self):
        return (
            self.pages['Destination'].configurable,
            self.pages['Program'].configurable,
            self.pages['Calculation'].configurable,
        )
        
        