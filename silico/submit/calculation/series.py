from silico.submit.calculation import Calculation_target
from silico.config.configurable.option import Option
from silico.exception.configurable import Configurable_exception

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
        type = tuple
    )
            
    def expand(self, calculations):
        """
        Expand this calculation target if it represents multiple real calcs.
        
        For Calculation_series, we return the calcs we represent.
        
        :return: A list of ready-to-go calculation targets.
        """
        calcs = []
        
        for tag_path in self.calculation_IDs:
            try:
                calcs.append(calculations.resolve(tag_path))
                
            except Exception as e:
                raise Configurable_exception(self, "Could not expand to real calculation with TAG path '{}'".format(tag_path)) from e
        
        return calcs