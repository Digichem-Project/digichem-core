from silico.submit.calculation import Calculation_target
from silico.config.configurable.option import Option
from silico.config.configurable.exception import Configurable_exception
from silico.submit.memory import Memory
from silico.config.configurable.options import Options
from silico.config.configurable.identifier import Identifier
from silico.config.configurable.util import getopt, setopt


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
        type = Identifier,
        list_type = list,
    )
    
    performance = Options(help = "Options for controlling functional aspects of the calculation.",
        memory = Option(help = "Force each calculation in this series to use this amount of memory.", default = None, type = Memory),
        num_cpu = Option(help = "Force each calculation in this series to use this many CPUs.", default = None, type = int)
    )
    post_process = Options(name = "post", help = "Options that control post processing of the calculation results.",
        write_summary = Option(help = "Whether to write Silico summary text files to the 'Results' folder at the end of each calculation. Set to False to disable for all calculations.", default = None, choices = (True, False, None)),
        write_report = Option(help = "Whether to write a Silico PDF report to the 'Report' folder at the end of each calculation. Set to False to disable for all calculations.", default = None, choices = (True, False, None)),
        write_combined_report = Option(help = "Whether to write the combined report for this series calculation", default = True, type = bool),
        combined_report_name = Option(help = "The name to use for the folder in which the combined report for this series calculation will be written. If not specified, the name of this series will be used.", default = None),
    )
    
            
    @property
    def combined_report_name(self):
        return self.post_process['combined_report_name'] if self.post_process['combined_report_name'] is not None else self.meta['name']
    
    def expand(self, calculations):
        """
        Expand this calculation target if it represents multiple real calcs.
        
        For Calculation_series, we return the calcs we represent.
        
        :return: A list of ready-to-go calculation targets.
        """
        calcs = []
        
        for identifer in self.calculation_IDs:
            try:
                # Get the calc we represent.
                # TODO: It feels strange that resolve cannot accept an identifier.
                calc = calculations.resolve(identifer.value)
                
            except Exception as e:
                raise Configurable_exception(self, "Could not expand to real calculation with TAG path '{}'".format(identifer.value)) from e
                
            # Overwrite specific properties if given.
            for prop in (("performance", "num_cpu"), ("performance", "memory"), ("post_process", "write_summary"), ("post_process", "write_report")):
                if getopt(self, *prop) is not None:
                    setopt(calc, *prop, getopt(self, *prop))
            
            calcs.append(calc)
        
        return calcs