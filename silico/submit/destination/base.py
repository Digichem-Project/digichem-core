from silico.submit import Configurable_target
from silico.submit.structure.directory import Calculation_directory
from uuid import uuid4
from silico.exception import Submission_error
from silico.config.configurable.option import Option

class Destination_target(Configurable_target):
    """
    Top-level class for classes that implement a destination to where calculations can be submitted.
    
    'Destinations' define how and where a calculation is run, but don't manage the calculation software itself (see silico.submit.program for those definitions).
    Destinations, for example, handle submission to a scheduling software (SLURM, TORQUE etc), running as a daemon, submission to a networked server etc.
    """
    
    CLASS_HANDLE = ("destination",)
            
    @property
    def unique_name(self):
        """
        Get a name that is unique for this calculation instance.
        
        Some destinations may provide their own get_unique_name() methods; this default implementation returns a random string with a very low collision chance.
        """
        if getattr(self, "_unique_name", None) is None:
            self._unique_name = uuid4().hex
            
        return self._unique_name
        
    @property
    def status(self):
        """
        This method is called to get 'status' about this destination.
        
        Status is always a string, but what it contains depends entirely on the Destination_target.
        It is used, for example, to report the number of free nodes for SLURM.
        
        This default implementation raises NotImplementedError.
        """
        raise NotImplementedError
    
    
    #############################
    # Class creation mechanism. #
    #############################
    
    class _actual(Configurable_target._actual):
        """
        Inner class for destinations.
        """
        
        def __init__(self):
            """
            Constructor for destination objects.
            """
            self.program = None
            self.calc_dir = None

        def submit(self):
            """
            Submit this destination.
            
            This default implementation creates the required directory structure.
            """
            # Set output directory.
            self.calc_dir = Calculation_directory.from_calculation(self.program.calculation)
            
            # Now create it.
            try:
                self.calc_dir.create_structure(True)
            except Exception:
                raise Submission_error(self, "could not create directory structure; try setting a different output directory ('-O')")
            
    
    