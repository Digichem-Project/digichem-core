

class Image_setup():
    """
    Top-level class for classes that setup images for report objects.
    """
    
    def __init__(self, report, *, metadata, do_orbitals = None, do_spin = None, options, calculation = None):
        """
        Constructor for Image_setup objects.
        
        :param report: The report object we will setup images for.
        :param metadata: The metadata corresponding to the (sub) calculation we will make images from.
        :param do_orbitals: Whether to generate orbitals.
        :param do_spin: Whether to generate spin density plots.
        :param options: Dictionary of config options.
        :param calculation: Optional calculation which will be used as a template for new calculations to generate new images.
        """
        self.report = report
        self.metadata = metadata
        self.options = options
        self.calculation = calculation
        
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