from silico.report import Report

class Gaussian_report(Report):
    """
    A specialised report object for processing Gaussian results.
    """
    
    def __init__(self, results, *, chk_file_path = None, fchk_file_path = None,):
        """
        Constructor for Gaussian reports.
        
        One (or both) of either chk_file_path or fchk_file_path must be given in order for orbital/molecular images to be rendered.
        
        :param results: A result_set object.
        :param chk_file_path: Optional path to a chk_file.
        :param fchk_file_path: Optional path to an fchk_file.
        """
        self.chk_file_path = chk_file_path
        self.fchk_file = fchk_file_path
        
        # Set up our image makers.
        images = {}
        
        
        # First
        
        
        # Continue in parent.
        super().__init__(results)