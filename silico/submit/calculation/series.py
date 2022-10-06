from silico.submit.calculation import Calculation_target
from silico.config.configurable.option import Option
from silico.exception.configurable import Configurable_exception
from silico.submit.memory import Memory


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
        type = None,
        list_type = list,
        # TODO: Need to add some way to actually edit properly; each item can be an int, a string or a list.
        no_edit = True
    )
    
    memory = Option(help = "Force each calculation in this series to use this amount of memory", default = None, type = Memory)
    num_cpu = Option(help = "Force each calculation in this series to use this many CPUs", default = None, type = int)
    _combined_report_name = Option("combined_report_name", help = "The name to use for the folder in which the combined report for this series calculation will be written", default = None)
            
    @property
    def combined_report_name(self):
        return self._combined_report_name if self._combined_report_name is not None else self.meta['name']
    
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
            if self.num_cpu is not None:
                calc.num_cpu = self.num_cpu
        
            if self.performance['memory'] is not None:
                # NOTE: memory is an object, this will set the same object reference for all sub-calcs.
                calc.performance['memory'] = self.performance['memory']
            
            calcs.append(calc)
        
        return calcs