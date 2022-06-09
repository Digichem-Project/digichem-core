# General imports.
import urwid
import yaml

# Silico imports.
from silico.interface.urwid.method.edit import Method_target_editor
from silico.interface.urwid.layout import Pane
from silico.interface.urwid.misc import Tab_pile
from silico.submit.base import Method_target, parse_method_from_file
from silico.interface.urwid.popup import Method_target_picker , Output_edit
from silico.interface.urwid.setedit.edit import Setedit_editor_mixin,\
    Paginated_settings_browser
from silico.submit.destination.base import Destination_target
from silico.submit.calculation.base import Calculation_target
from silico.submit.program.base import Program_target
from silico.interface.urwid.swap.swappable import Swappable
import silico.logging
from silico.interface.urwid.file.browser import File_selector
from silico.interface.urwid.dialogue import Confirm_or_cancel_dialogue


class Method_target_builder(Tab_pile, Setedit_editor_mixin):
    """
    """
    
    def __init__(self, top, parent_class, method_library, method_type, initial = None, initially_from_library = None):
        """
        """
        self.top = top
        self.parent_class = parent_class
        self.method_library = method_library
        self.method_type = method_type
        
        # Whether to start by showing the method editor or library picker.
        if initially_from_library is None:
            initially_from_library = not isinstance(initial, Method_target) and initial is not None
        
        # A checkbox which will change the main body contents depending on its state.
        self.from_library_widget = urwid.CheckBox("Select an existing {} from the library or create new?".format(method_type), initially_from_library, on_state_change = self.switch_target_editor)
        
        # An editor which can be used to change settings.
        self.editor = Method_target_editor(self.top, self.parent_class, initial if not initially_from_library else None)
        
        # A widget which can pick an existing method from the library.
        if initially_from_library:
            # Get the configurable referred to by our ID.
            try:
                initial = method_library.resolve(initial) if initial is not None else None
            
            except Exception:
                silico.logging.get_logger().error("Failed to resolve {} identified by ID '{}'".format(method_type, initial), exc_info = True)
                initial = None
            
        self.library_picker_widget = Method_target_picker(top, method_library, initial = initial if initially_from_library else None)
        self.library_picker_pane = Pane(self.library_picker_widget, "Choose from library")
        
        # The widget whose contents will change depending on the value of the from_library widget.
        self.body = urwid.WidgetPlaceholder(urwid.Filler(self.library_picker_pane, valign = "top") if initially_from_library else self.editor)
        
        # Some attributes to enable rollback.
        self.previous_from_library = initially_from_library
        self.previous_library_configurable = self.library_picker_widget.value
        
        super().__init__([
            ("pack", Pane(self.from_library_widget, "From library")),
            self.body
        ])
        
    
    @property
    def from_library(self):
        """
        Whether the method_target being represented by this editor is picked from the method library.
        """
        return self.from_library_widget.state
        
    def switch_target_editor(self, check_box, new_state):
        """
        Switch between the two types of editor that this widget supports based on the current value of the checkbox widget.
        """
        if not new_state:
            # Not ticked, manual widget.
            self.body.original_widget = self.editor
        
        else:
            # Ticked, library widget.
            self.body.original_widget = urwid.Filler(self.library_picker_pane, valign = "top")
            
    def dump(self):
        """
        Get the current value of the method_target that is being edited as a serializable object (probably a str, dict or list).
        """
        try:
            # First, if we represent a value picked from the library, return the tag path leading to it.
            if self.from_library:
                tag_hierarchy = self.library_picker_widget.value.tag_hierarchy
                if len(tag_hierarchy) == 1:
                    return tag_hierarchy[0]
                
                else:
                    return tag_hierarchy
            
            # We have an actual configurable.
            return self.editor.configurable.dump()
        
        except AttributeError:
            if (self.from_library and self.library_picker_widget.value is None) or (not self.from_library and self.editor.configurable is None):
                # No configurable has been chosen yet.
                #raise Exception("No {} has been chosen".format(self.method_type)) from None
                silico.logging.get_logger().warning("No {} has been chosen; this part of the method will be left blank".format(self.method_type))
                return None
            
            else:
                raise
    
    def refresh(self):
        """
        Refresh each of the pages of options.
        """
        # Not sure it makes sense to refresh here?
        raise NotImplementedError("refresh has no meaning for '{}'".format(type(self).__name__))
    
    def save(self, validate = True):
        """
        Save the changes made to this part of the method.
        """
        if self.from_library:
            # Check the chosen value is valid.
            if validate and self.library_picker_widget.value is not None:
                # This method will throw an exception if the tag path cannot be resolved.
                self.method_library.path_by_tags(self.library_picker_widget.value.tag_hierarchy)
            
        else:
            self.editor.save(validate)
            
        # Update rollback attributes.
        self.previous_from_library = self.from_library
        self.previous_library_configurable = self.library_picker_widget.value
            
    def validate_setedits(self):
        """
        Validate each of the pages of options (without saving changes first).
        """
        if self.from_library:
            # Check the chosen value is valid.
            # This method will throw an exception if the tag path cannot be resolved.
            self.method_library.path_by_tags(self.library_picker_widget.value.tag_hierarchy)
            
        else:
            self.editor.validate_setedits()
            
    def discard(self):
        """
        Discard any changes made.
        """
        self.editor.discard()
        self.from_library_widget.set_state(self.previous_from_library)
        self.library_picker_widget.value = self.previous_library_configurable


class Method_builder(Swappable, Setedit_editor_mixin):
    """
    A widget that can be used to build/create a new calculation method.
    """
    
    def __init__(self, top, method_library, destination = None, program = None, calculation = None, file_path = None):
        """
        Constructor for Method_builder objects.
        
        :param top: The top-most widget to use for swapping.
        :param method_library: A loaded method library, containing destinations, programs and calculations as attributes.
        :param destination: Optional destination (either an ID or an object) to edit by default.
        :param destination: Optional program (either an ID or an object) to edit by default.
        :param destination: Optional calculation (either an ID or an object) to edit by default.
        :param file_path: A default file path to show.
        """
        calculation_page = Method_target_builder(top, Calculation_target, method_library = method_library.calculations, method_type = "calculation", initial = calculation)
        program_page = Method_target_builder(top, Program_target, method_library = method_library.programs, method_type = "program", initial = program)
        destination_page = Method_target_builder(top, Destination_target, method_library = method_library.destinations, method_type = "destination", initial = destination)
        
        # We show in reverse order, because typically the calculation is what the user wants to change anyway.
        self.editor = Paginated_settings_browser({
            "Calculation": calculation_page,
            "Program": program_page,
            "Destination": destination_page
        })
        
        # The save location.
        self.output_widget = Output_edit(top, file_path, default_file_name = "method.sim")
        
        super().__init__(top, Tab_pile([
            self.editor,
            ('pack', Pane(self.output_widget, "Save location"))
        ]))
        
    def submit_callback(self):
        """
        Save the method to file.
        """
        try:
            self.write()
            silico.logging.get_logger().info("Saved method to file '{}' successfully".format(self.output_widget.value))
        
        except Exception:
            silico.logging.get_logger().error("Failed to write method to file", exc_info = True)
            return False
            
    def write(self):
        """
        Write the method to file.
        """
        # Panic if there's no file to save to.
        if self.output_widget.value is None:
            raise Exception("No output location chosen")
            
        
        # First, save changes.
        self.save()
        
        # Build a method dict.
        method = {
            'destination': self.editor.page['Destination'].dump(),
            'program': self.editor.page['Program'].dump(),
            'calculation': self.editor.page['Calculation'].dump()    
        }
        
        # Open the file pointed at by the output widget.
        try:
            with open(self.output_widget.value, "wt") as method_file:
                yaml.dump(method, method_file)
        
        except Exception as e:
            raise Exception("Failed to open method file for writing") from e
            
    def cancel_callback(self):
        return self.discard()
    
    def refresh(self):
        """
        Refresh each of the pages of options.
        """
        return self.editor.refresh()
            
    def save(self, validate = True):
        """
        Save changes made to the method.
        """
        return self.editor.save(validate)

    def validate_setedits(self):
        """
        Validate each of the pages of options (without saving changes first).
        """
        return self.editor.validate_setedits()

    def discard(self):
        """
        Discard any changes made.
        """
        return self.editor.discard()


class Method_builder_menu(Swappable):
    """
    A widget (menu) which provides access to a method builder.
    """
    
    def __init__(self, top, method_library):
        """
        :param top: The top-most widget to use for swapping.
        :param method_library: A loaded method library, containing destinations, programs and calculations as attributes.
        """
        self.method_library = method_library
        
        # The method builder widget that can be switched to create/edit methods.
        self.builder = None
        
        # A file browser to use to pick existing method files.
        self.file_selector = File_selector(top, pane_title = "Select method file", can_choose_folders = False, can_choose_multiple = False)
        
        self.menu = urwid.Pile([
            ('pack', urwid.AttrMap(urwid.Button("Create a new method", lambda button: self.create_new()), "button--small", "button--small--focus")),
            ('pack', urwid.AttrMap(urwid.Button("Load a method from file", lambda button: self.popup_file_browser()), "button--small", "button--small--focus")),
        ])
        self._done_add_resume_button = False
        
        super().__init__(top, urwid.Filler(urwid.Padding(Pane(self.menu, "Method Builder"), width = 45, align = "center")))
        
    def cancel_callback(self):
        self.top.back()
        
    def popup_file_browser(self):
        """
        Open the file browser widget that can be used to select an existing method file to edit.
        """
        self.top.swap_into_window(self.file_selector, submit_callback = self.load_from_file)
        
    def add_resume_editing_button(self):
        """
        Add a button to the central menu that switches back to the method builder.
        """
        if not self._done_add_resume_button:
            self.menu.contents.insert(0,
                (urwid.AttrMap(urwid.Button("Continue editing", lambda button: self.swap_to_builder()), "button--small", "button--small--focus"), self.menu.options(height_type = "pack")) 
                )
            self.menu.focus_position = 0
        self._done_add_resume_button = True
        
    def swap_to_builder(self):
        """
        Open the method builder.
        """
        self.top.swap_into_window(self.builder, submit_callback = self.builder.submit_callback)
    
    def load_from_file(self):
        """
        Load a method from a file.
        """
        # If no file selected, do nothing.
        if len(self.file_selector.selected) == 0:
            return
        
        method_file_path = self.file_selector.selected[0]
        self.file_selector.reset()
        
        # First, load the method from file (or try to anyway).
        try:
            method = parse_method_from_file(method_file_path, self.method_library, resolve = False)
            
        except Exception:
            silico.logging.get_logger().error("Failed to parse method file '{}'".format(method_file_path), exc_info = True)
            return
        
        # TODO: This swapping is a bit clumsy.
        self.top.back()
        self.change_builder(Method_builder(self.top, self.method_library, method[0], method[1], method[2], method_file_path))
        
        return False
        
    def create_new(self):
        """
        Create a new, blank method.
        """
        self.change_builder(Method_builder(self.top, self.method_library))
        
    def change_builder(self, new_builder):
        """
        Change the method builder, prompting the user if they wish to discard changes if a builder has already been set.
        """
        def change_callback():
            self.builder = new_builder
            self.add_resume_editing_button()
            self.swap_to_builder()
            
        # If a builder has already been set, prompt the user whether they're sure they want to change.
        if self.builder is not None:
            self.top.popup(Confirm_or_cancel_dialogue("Discard Method?", "There is already a method that is open for editing. Opening a new method will discard any unsaved changes to this old method. Are you sure you want to do this?", self.top, error = True, submit_callback = change_callback))
        
        else:
            # No builder currently set, change without prompting.
            change_callback()
        
            
        