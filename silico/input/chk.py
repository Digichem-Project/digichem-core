from silico.input import Input_file


class Chk_input(Input_file):
    """
    An Input_file type for submitting Gaussian calculations from an existing chk file.
    """
    
    def __init__(self, chk_file, name = None):
        """
        Constructor for Chk_input objects.
        
        Note that unlike real coord input objects, chk files are not cached in memory (because they are huge).
        Normally this is fine, but it will cause problems if the chk file is moved before the calculation starts,
        or if the calculation is run on another machine.
        """
        # Charge and mult is saved on the chk file.
        super().__init__(charge = None, multiplicity = None, name = name, file_name = chk_file)
        self.chk_file = chk_file
        