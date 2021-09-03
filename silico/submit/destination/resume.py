from silico.exception.uncatchable import Submission_paused
import dill
#import pickle as dill
from pathlib import Path
from silico.submit.structure.flag import Flag
from silico.submit.destination import Destination_target

class Resumable_destination(Destination_target):
    """
    Mixin class for destination that are resumable.
    
    Resumable here means that program execution stops during the submission process. Submission is then 'resumed' in a new process immediately before the calculation proper begins.
    This mechanism is required by several destinations, for example SLURM (calculation has to occur from the SLURM node) and SSH (calculation occurs on a different machine entirely).
    
    The 'resume' is achieved via pickle, so your class should be picklable.
    """
    
    CLASS_HANDLE = ()
    
    ############################
    # Class creation mechanism #
    ############################
    
    class _actual(Destination_target._actual):
        """
        Inner class for SLURM.
        """
        def __init__(self, *args, **kwargs):
            """
            Constructor for Resumable_destination objects.
            
            :param output: Path to the file that we will write to resume from.
            """
            super().__init__(*args, **kwargs)
            # A flag keeping track of which side of the pickle reload we are.
            self._resumed = False
            
        def pre(self):
            """
            Method called prior to pausing, inheriting classes should define setup here.
            
            This default implementation does nothing.
            """
            pass
        
        def post(self):
            """
            Method called after pausing, inheriting classes should define setup here.
            
            This default implementation does nothing.
            """
            pass
            
        def submit(self):
            """
            Submit to this resumable destination.
            
            In resumable destinations, this method will get called twice (automatically) during the submission process.
            This method can raise Submission_paused exceptions as part of this process,
            you should not go out of you way to catch these exceptions unless you know what you are doing
            (and you should almost certainly be re-raising once you are done).
            """
            # If we have not yet resumed, create our directory structure.
            # We need to do this before the resume because this is where we'll write our SLURM batch file.
            if not self._resumed:
                # First, call our parent (creates directories).
                super().submit()
                
                # Call user defined pre function.
                self.pre()

                # Pause here.
                self.pause()
                
                # Set our pending flag (we will be going into a queue).
                self.calc_dir.set_flag(Flag.PENDING)
                
                # Call user defined post function.
                self.post()
                
                # Use a special 'exception' to prevent normal submission.
                raise Submission_paused()
            
            # If we get this far, then we have resumed and can continue as normal.
            # Delete the pending flag.
            self.calc_dir.del_flag(Flag.PENDING)
            
        @property
        def resume_file_path(self):
            """
            Path to our pickled resume file. 
            """
            return Path(self.calc_dir.input_directory, "silico.resume.pickle")
            
        def pause(self):
            """
            The first part of the 'resume' mechanism, pause() should set-up the class for resuming later. Typically this involves writing to a pickle file.
            """
            with open(self.resume_file_path, "bw") as pickle_file:
                dill.dump(self, pickle_file)
            
        def resume(self):
            """
            The second part of the 'resume' mechanism, resume() is called after the pickle file has been re-loaded. Execution should continue from here.
            """
            self._resumed = True
            
            # Continue submitting.
            self.program.calculation.submit()
        