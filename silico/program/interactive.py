# General imports.
from openbabel import pybel
import sys

# Silico imports.
from silico.program.base import Program
from silico.program.submit import Submit_program
from silico.program.report import Report_program
from silico.program.result import Result_program
from silico.program.status import Status_program
from silico.program.convert import Convert_program
from silico.interface.urwid.main.base import Silico_window, Output_catcher
import silico.base


class Interactive_program(Program):
    """
    The Silico interactive program.
    """
    
    name = "Silico (interactive)"
    command = "interactive"
    description = "run Silico interactively with a console interface"
    help = "Run interactively"
        
    def __init__(self, args, config, logger, initial = None):
        """
        Constructor for the interactive program.
        """
        super().__init__(args, config, logger)
        # Save these for later.
        self.init_args = (args, config, logger)
        
        # A list of program instances.
        # We init this with out initial program, which has already been constructed.
        self.initial = initial
        self._programs = {type(initial): initial}
        
        # All the programs we can switch between.
        self.program_classes = [Submit_program, Report_program, Result_program, Status_program, Convert_program]

    def get_program(self, prog_cls):
        """
        Get an instance of a sub program.
        
        This method will first create a new instance of the program class before caching it so future calls will return the same instance.
        
        :param prog_cls: The program class to get.
        :returns: A program instance.
        """
        if prog_cls not in self._programs:
            self._programs[prog_cls] = prog_cls.load_from_argparse(*self.init_args)
            
        return self._programs[prog_cls]
        
    def main(self):
        """
        Logic for our program.
        """
        window = Silico_window(self)
        
        # Setup stdout redirection.
        with Output_catcher(window, False) as stdout, Output_catcher(window, True) as stderr:
            sys.stdout = stdout
            sys.stderr = stderr
            
            # Also update our logging output, which is handled separately.
            old_log_stream = silico.base.LOGGING_HANDLER.stream
            silico.base.LOGGING_HANDLER.stream = stderr
            
            # Turn of warnings from pybel, they're probably not useful anyway and get dumped randomly in the screen.
            log_level = pybel.ob.obErrorLog.GetOutputLevel()
            pybel.ob.obErrorLog.SetOutputLevel(-1)
            
            try:
                return window.run_loop(self.config.urwid_palette)
                
            finally:
                # Restore stdout.
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
                silico.base.LOGGING_HANDLER.stream = old_log_stream
                
                # And obabel logging.
                pybel.ob.obErrorLog.SetOutputLevel(log_level)
                
        
        
        