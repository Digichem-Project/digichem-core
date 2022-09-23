from pathlib import Path

from silico.input.base import Input_file


class Calculation_directory_input(Input_file):
    """
    Class for using a directory containing a previous calculation as input.
    
    Typically this input type is used for restarting an old calculation.
    If you just want to submit a new calculation using the output coordinates of a previous calculations, use the normal Input_coords class instead.
    """
    
    def __init__(self, calculation_directory, name = None):
        """
        """
        super().__init__(charge = None, multiplicity = None, name = name, file_name = calculation_directory)
        self.calculation_directory = Path(calculation_directory)