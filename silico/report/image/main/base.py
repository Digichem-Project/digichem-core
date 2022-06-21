class Cube_setup():
    """
    ABC for classes that create cubes for reports.
    """
        
    def setup(self, output_dir, output_name):
        """
        Perform setup.
        
        Calling this method will set cube objects in the parent report.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        # This default implementation does nothing.
        raise NotImplementedError()
        
    def cleanup(self):
        """
        Remove any intermediate files that may have been created by this object.
        """
        pass


class Partial_cube_setup(Cube_setup):
    """
    ABC for classes that setup cube maker objects for generating orbital and spin density images.
    """
    
    def __init__(self, report, *, result, do_orbitals = None, do_spin = None, options, calculation = None):
        """
        Constructor for Orbital_and_spin_cube_setup objects.
        
        These objects process the results from a single sub-calculation at a time (in the case where multiple calculations have been merged).
        
        :param report: The report object we will setup images for.
        :param result: The result corresponding to the (sub) calculation we will make images from.
        :param do_orbitals: Whether to generate orbitals.
        :param do_spin: Whether to generate spin density plots.
        :param options: Dictionary of config options.
        :param calculation: Optional calculation which will be used as a template for new calculations to generate new images.
        """
        self.report = report
        self.result = result
        self.options = options
        self.calculation = calculation