from silico.submit import Configurable_target
from silico.submit.structure.directory import Calculation_directory
from uuid import uuid4
from silico.exception import Submission_error

class Method_target(Configurable_target):
    """
    Top-level class for classes that implement a method of submitting calculations.
    
    'Methods' define how and where a calculation is run, but don't manage the calculation software itself (see silico.submit.program for those definitions).
    Methods, for example, handle submission to a scheduling software (SLURM, TORQUE etc), running as a daemon, submission to a networked server etc.
    """
    
    CLASS_HANDLE = ("method",)
        
    @property
    def unique_name(self):
        """
        Get a name that is unique for this calculation instance.
        
        Some methods may provide their own get_unique_name() methods; this default implementation returns a random string with a very low collision chance.
        """
        if getattr(self, "_unique_name", None) is None:
            self._unique_name = uuid4().hex
            
        return self._unique_name
        
    @property
    def status(self):
        """
        This method is called to get 'status' about this method.
        
        Status is always a string, but what it contains depends entirely on the Method_target.
        It is used, for example, to report the number of free nodes for SLURM.
        
        This default implementation raises NotImplementedError.
        """
        raise NotImplementedError
    
    
    #############################
    # Class creation mechanism. #
    #############################
    
    class _actual(Configurable_target._actual):
        """
        Inner class for methods.
        """
        
        def __init__(self):
            """
            Constructor for method objects.
            """
            self.program = None
            self.calc_dir = None

        def submit(self):
            """
            Submit this method.
            
            This default implementation creates the required directory structure.
            """
            # Set output directory.
            self.calc_dir = Calculation_directory.from_calculation(self.program.calculation)
            
            # Now create it.
            try:
                self.calc_dir.create_structure(True)
            except Exception:
                raise Submission_error(self, "could not create directory structure; try setting a different output directory ('-O')")
            
    
    