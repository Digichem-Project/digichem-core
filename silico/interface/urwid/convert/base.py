# General imports.
import urwid

# Silico imports.
from silico.interface.urwid.program import Program_view
from silico.interface.urwid.edit.popup import File_edit, Output_edit,\
    Choices_edit
from silico.file.input import Silico_coords
from silico.interface.urwid.misc import Tab_pile, IntEditZero
from silico.interface.urwid.layout import Pane


class Convert_interface(Program_view):
    """
    Class for controlling interface to result program.
    """
    
    @classmethod
    def input_formats(self):
        """
        An ordered list of input formats that are supported.
        """
        formats = list(Silico_coords.input_formats().keys())
        formats.sort(key = lambda item: item.upper())
        formats.insert(0, None)
        return formats
    
    @classmethod
    def output_formats(self):
        """
        An ordered list of input formats that are supported.
        """
        formats = list(Silico_coords.output_formats().keys())
        formats.sort(key = lambda item: item.upper())
        formats.insert(0, None)
        return formats
    
    def __init__(self, window, program):
        """
        Constructor for calculation submitter objects.
        
        :param window: The program window we are being rendered in.
        :param program: A program object.
        """
        # Keep track of our individual widgets.
        # The file to convert.
        self.input_widget = File_edit(window.top, program.args.input_file, "Input File", can_choose_folders = False)
        # The format of the input file.
        self.input_format_widget = Choices_edit(window.top, self.input_formats(), program.args.input_format, "Input format")
        
        # Multiplicity and charge.
        self.charge_widget = IntEditZero(("body", "Charge: "), program.args.charge)
        self.multiplicity_widget = IntEditZero(("body", "Multiplicity: "), program.args.multiplicity) 
        
        # Whether to gen3D.
        self.gen3D_widget = urwid.CheckBox(("body", "Gen3D (whether to generate 3D coordinates from 2D files)"), program.args.gen3D)
        # The location to write the formatted results to.
        self.output_widget = Output_edit(window.top, program.args.output_file if program.args.output_file != "-" else None)
        # The output format.
        self.output_format_widget = Choices_edit(window.top, self.output_formats(), program.args.output_format, "Output format")
        
        super().__init__(window, program)
        
    def get_body(self):
        """
        Get the widget used to display the body of this program.
        
        :returns: An urwid widget to display.
        """
        return Tab_pile([
            Pane(urwid.Filler(
                urwid.Pile([
                    urwid.Columns([
                        ('pack', urwid.Text("Input file:")),
                        urwid.AttrMap(self.input_widget, "editable"),
                    ], dividechars = 1),
                    
                    urwid.Columns([
                        ('pack', urwid.Text("Input format:")),
                        urwid.AttrMap(self.input_format_widget, "editable"),
                    ], dividechars = 1),
                ]),
            valign = "top"), "Input Options"),
            
            Pane(urwid.Filler(
                urwid.Pile([
                    urwid.Columns([
                        ('pack', urwid.AttrMap(self.gen3D_widget, "editable")),
                    ], dividechars = 1),
                    
                    urwid.Columns([
                        urwid.AttrMap(self.charge_widget, "editable"),
                        urwid.AttrMap(self.multiplicity_widget, "editable"),
                    ], dividechars = 1),
                    
                    urwid.Columns([
                        ('pack', urwid.Text("Output format:")),
                        urwid.AttrMap(self.output_format_widget, "editable"),
                    ], dividechars = 1),
                    
                    urwid.Columns([
                        ('pack', urwid.Text("Save location:")),
                        urwid.AttrMap(self.output_widget, "editable"),
                    ], dividechars = 1)
                
                ]),
            valign = "top"), "Output Options")
        ])
        
    def submit(self):
        self.program.logger.info("Converting file...")
        return super().submit()
    
    def post(self):
        self.program.logger.info("Done converting file")
        
    def setup(self):
        """
        Setup our program to ready it to be run.
        """
        self.program.args.input_file = self.input_widget.value
        self.program.args.input_format = self.input_format_widget.value
        
        self.program.args.charge = self.charge_widget.value()
        self.program.args.multiplicity = self.multiplicity_widget.value()
        
        self.program.args.output_format = self.output_format_widget.value
        self.program.args.output_file = self.output_widget.value
        self.program.args.gen3D = self.gen3D_widget.get_state()
        
        