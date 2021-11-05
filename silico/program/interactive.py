# General imports.

# Silico imports.
from silico.program.base import Program
from silico.program.submit import Submit_program
from silico.program.report import Report_program
from silico.program.result import Result_program
from silico.program.status import Status_program
from silico.program.convert import Convert_program
from silico.interface.urwid.main.base import Silico_window, Output_catcher
import sys
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
        
        # Our list of sub programs.
        self.programs = {
            prog_cls.command: prog_cls.load_from_argparse(args, config, logger) 
            for prog_cls 
            in [Submit_program, Report_program, Result_program, Status_program, Convert_program]
            if prog_cls != type(initial)
        }
        
        self.initial = initial
        if initial is not None:
            self.programs[initial.command] = initial
        
        
    def main(self):
        """
        Logic for our program.
        """
        window = Silico_window(programs = self.programs.values(), initial = self.initial)
        
        # Setup stdout redirection.
        with Output_catcher(window.top, False) as stdout, Output_catcher(window.top, True) as stderr:
            sys.stdout = stdout
            sys.stderr = stderr
            # Also update our logging output, which is handled separately.
            old_log_stream = silico.base.LOGGING_HANDLER.stream
            silico.base.LOGGING_HANDLER.stream = stderr
            
            try:
                return window.run_loop(self.config.urwid_palette)
                
            finally:
                # Restore stdout.
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
                silico.base.LOGGING_HANDLER.stream = old_log_stream
        
        
        