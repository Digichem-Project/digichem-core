# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.row_list.base import Row_browser, Row_item,\
    Row_widget
from silico.interface.urwid.result.browser import Result_selector
import silico.logging
from silico.parser.base import parse_multiple_calculations, parse_calculation
from silico.result.alignment.base import Alignment


class Result_widget(Row_widget):
    """
    Widget for displaying Result_item objects.
    """
    
    def __init__(self, *args, **kwargs):
        self.charge_edit = None
        self.mult_edit = None
        super().__init__(*args, **kwargs)

    def load_inner(self):
        """
        Load the widget we'll use to display our main body.
        """
        return urwid.Columns([
            ('weight', 1, urwid.Text(self.row_item._value.metadata.identity_string)),
            # TODO: Careful, what if result has no file?
            ('weight', 2, urwid.Text(str(self.row_item._value.metadata.log_files[0])))
        ])


class Result_item(Row_item):
    """
    An item that appears in a Result_list, represents a loaded coordinates file.
    """
    
    def __init__(self, result, row_list, movable = True,):
        """
        Constructor for Result_item objects.
        
        :param result: A Result_set object.
        :param row_list: Our parent row_list object.
        :param movable: Whether we should show rearranging buttons.
        """
        super().__init__(result, row_list, movable = movable)
        
    def load_widget(self):
        """
        Load the widget we'll use for display.
        """
        return Result_widget(self)


class Result_list(Row_browser):
    """
    A widget for displaying a list of parsed results.
    """
    
    def __init__(self, top, start_dir = None, rearrangeable = True, initial_results = None, subprocess_init = None):
        """
        Constructor for Result_list objects.
        
        :param top: A top-most widget to use for display.
        :param start_dir: The directory to show when first opening the file selector.
        :param rearrangeable: Whether the list of loaded results can be rearranged by the user.
        :param initial_results: A list of parsed results to show by default.
        :param subprocess_init: A function to call to initialise subprocesses when parsing results in parallel.
        """
        top = top
        self.subprocess_init = subprocess_init
        
        # Widget for choosing files, stored at self.selector
        file_selector = Result_selector(top, start_dir, title = "Select Result Files")
        
        # Add our starting files.
        initial_results = [] if initial_results is None else initial_results
        
        super().__init__(file_selector, top = top, rearrangeable = rearrangeable, initial = initial_results)
        
    def add_to_list(self):
        """
        Add selected files to our list.
        """
        # Get the logger item we'll use for communication.
        logger = silico.logging.get_logger()
        logger.info("Parsing calculation results")
        
        # Get the chosen alignment class.
        alignment = Alignment.from_class_handle(self.selector.alignment)
        
        # Parse (in parallel).
        # We currently do not do this in parallel because logging to urwid gets messed up...
        #results = parse_multiple_calculations(*self.selector.selected, alignment_class = alignment, init_func = self.subprocess_init, processes = self.selector.num_CPUs)
        results = [parse_calculation(selected, alignment_class = alignment) for selected in self.selector.selected]
        
        # Add each.
        for result in results:
            self.add_row(result)
                 
        # Clear the selected files.
        self.selector.reset()
         
        logger.info("Finished parsing calculation results")
        
    def load_row(self, value):
        return Result_item(value, self)