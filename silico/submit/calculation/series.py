from silico.submit.calculation import Calculation_target
from silico.config.configurable.option import Option
from silico.exception.configurable import Configurable_exception
from silico.submit.base import Memory

class Calculation_series(Calculation_target):
    """
    A pseudo calculation that automatically runs a number of calculations in series.
    """
    
    CLASS_HANDLE = ["Series"]
    
    # Configurable options.
    calculation_IDs = Option(
        "calculations",
        help = "A list of calculations to perform in series",
        required = True,
        type = tuple,
        no_edit = True
    )
    
    memory = Option(help = "Force each calculation in this series to use this amount of memory", default = None, type = Memory, rawtype = str)
    num_CPUs = Option(help = "Force each calculation in this series to use this many CPUs", default = None, type = int)
            
    def expand(self, calculations):
        """
        Expand this calculation target if it represents multiple real calcs.
        
        For Calculation_series, we return the calcs we represent.
        
        :return: A list of ready-to-go calculation targets.
        """
        calcs = []
        
        for tag_path in self.calculation_IDs:
            try:
                # Get the calc we represent.
                calc = calculations.resolve(tag_path)
                
            except Exception as e:
                raise Configurable_exception(self, "Could not expand to real calculation with TAG path '{}'".format(tag_path)) from e
                
            # Overwrite the memory and CPUs if given.
            if self.num_CPUs is not None:
                calc.num_CPUs = self.num_CPUs
        
            if self.memory is not None:
                calc.memory = self.memory
            
            calcs.append(calc)
        
        return calcs